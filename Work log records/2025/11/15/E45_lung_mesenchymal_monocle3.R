setwd('~/project/03_3D_spatial/02_result/251115_E45_lung_Mesenchymal_monocle3/')
library(Seurat)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(dplyr)
library(patchwork)
library(stringr)
library(rjson)
library(data.table)
library(harmony)
library(randomcoloR)
obj.list = list()
stage_list <- c('E25','E35','E45')
for (i in 1:length(stage_list)){
  obj.list[[i]] <- readRDS(paste0('../../00_data/lung/',stage_list[i],'_new_lung.rds'))
  obj.list[[i]]$Stage <- stage_list[i]
}
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
head(obj)
table(obj$celltype_gut)
#res1 <- data.frame(obj$cell,obj$celltype_gut,obj$Stage)
#write.table(res1,'../../00_data/lung/E25_E35_45_lung_all_cellID.txt',sep = '\t',row.names = F)
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells'))

obj_sub1 <- subset(obj_sub,subset = Stage == "E45")
obj_sub1 <- NormalizeData(obj_sub1,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(obj_sub1,group.by = 'Spatial_snn_res.1',label = T,raster=FALSE,label.size = 5)+NoLegend()
FeaturePlot(obj_sub1,features = "NKX6-1")
all_markers <- FindAllMarkers(object = obj_sub1,group.by = 'Spatial_snn_res.1',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

color_mapping <- export_clusters_to_ply(obj_sub1,
                                        output_dir = "./E45_lung_mesenchymal_phy",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.1")

obj_sub2 <- subset(obj_sub1,subset = Spatial_snn_res.1 != "5")
obj_sub2$type <- "E45"
obj_sub2$type[obj_sub2$Spatial_snn_res.1 %in% c("6","10")] <- 'E45_single_slice'
table(obj_sub2$type)
obj_sub2[["Spatial"]] <- split(obj_sub2[["Spatial"]], f = obj_sub2$type)  ## 对时期间和单张芯片去批次
obj_sub2 <- NormalizeData(obj_sub2)
obj_sub2 <- FindVariableFeatures(obj_sub2)
obj_sub2 <- ScaleData(obj_sub2)
obj_sub2 <- RunPCA(obj_sub2)
obj_sub2 <- FindNeighbors(obj_sub2, dims = 1:30, reduction = "pca")
obj_sub2 <- IntegrateLayers(object = obj_sub2, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj_sub2[["Spatial"]] <- JoinLayers(obj_sub2[["Spatial"]])
obj_sub2 <- FindNeighbors(obj_sub2, reduction = "integrated.harmony", dims = 1:30)
obj_sub2 <- FindClusters(obj_sub2, resolution = 1)
obj_sub2 <- RunUMAP(obj_sub2, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)

DimPlot(obj_sub2,group.by = 'type',label = T,raster=FALSE,label.size = 8)+NoLegend()
FeaturePlot(obj_sub2,features = "FGF10")
color_mapping <- export_clusters_to_ply(obj_sub2,
                                        output_dir = "./E45_lung_mesenchymal_phy",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.1")

#####monocle3
#monocle3
data <- GetAssayData(obj_sub2,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj_sub2@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('cds.umap')
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj_sub2, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('int.umap')
p2

cds <- cluster_cells(cds,cluster_method = 'louvain')
#p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
#p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
#  ggtitle("label by partitionID")
#p = wrap_plots(p1, p2)
#p

## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p


p + geom_vline(xintercept = seq(1,4,0.5)) + geom_hline(yintercept = seq(0.5,4,0.5))
embed <- data.frame(Embeddings(obj_sub2, reduction = "umap"))
embed <- subset(embed, umap_1 > 1 & umap_1 < 4 & umap_2 > 0.5 & umap_2 < 4)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
#pdf('E35_E45_lung_main_bronchi_monocle3_pseudotime_umap.pdf',width = 7, height = 5)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1.5,
           group_label_size = 2,
          trajectory_graph_segment_size = 1.5)+ 
  theme_classic()+
  theme(axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16))
#dev.off()



Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=20)
Track_genes_sig <- Track_genes %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()

Track_genes_sig <- c('PTN','WNT5A')
#pdf('E35_E45_lung_main_bronchi_monocle3_top30Genes.pdf',width = 12,height = 8)
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
#dev.off()

saveRDS(obj_sub2,'E45_lung_mesenchymal_251116.rds')
