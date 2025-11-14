setwd('~/project/03_3D_spatial/02_result/251113_E35E45_lung_Mesenchymal_monocle3/')
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

#E45 单张芯片lung结果
cellID <- read.table('E35_45_lung_mesenchymal_cellID.txt',sep = '\t',header = T)
cellID1 <- cellID[cellID$obj_sub.Stage == "E45",]
cellID2 <- cellID1[cellID1$obj_sub.Spatial_snn_res.1 == "3",]
cellID2$obj_sub.cell
obj_sub@meta.data[cellID2$obj_sub.cell, "Stage"] <- "E45_single_slice"
table(obj_sub$Stage)

obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间和单张芯片去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.harmony", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)

head(obj_sub)
DimPlot(obj_sub,group.by = 'Spatial_snn_res.1',label = T,label.size = 10)
FeaturePlot(obj_sub,features = c('NKX6-1'))
res <- data.frame(obj_sub$cell,obj_sub$Spatial_snn_res.1,obj_sub$Stage)
write.table(res,'E35_45_lung_mesenchymal_cellID2.txt',sep = '\t',row.names = F)

obj_sub$newcelltype <- 'lung'
obj_sub@meta.data[cellID2$obj_sub.cell, "newcelltype"] <- "E45_single_slice"
DimPlot(obj_sub,group.by = 'newcelltype',label = T,label.size = 10)
obj_sub$newcelltype <- NULL

#删除低质量细胞群
obj_sub1 <- subset(obj_sub,subset = Spatial_snn_res.1 != "8")
DimPlot(obj_sub1,group.by = 'Spatial_snn_res.1',label = T,label.size = 10)
obj_sub1[["Spatial"]] <- split(obj_sub1[["Spatial"]], f = obj_sub1$Stage)  ## 对时期间和单张芯片去批次
obj_sub1 <- NormalizeData(obj_sub1)
obj_sub1 <- FindVariableFeatures(obj_sub1)
obj_sub1 <- ScaleData(obj_sub1)
obj_sub1 <- RunPCA(obj_sub1)
obj_sub1 <- FindNeighbors(obj_sub1, dims = 1:30, reduction = "pca")
obj_sub1 <- IntegrateLayers(object = obj_sub1, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj_sub1[["Spatial"]] <- JoinLayers(obj_sub1[["Spatial"]])
obj_sub1 <- FindNeighbors(obj_sub1, reduction = "integrated.harmony", dims = 1:30)
obj_sub1 <- FindClusters(obj_sub1, resolution = 1)
obj_sub1 <- RunUMAP(obj_sub1, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub1,group.by = 'Spatial_snn_res.1',label = T,label.size = 10)
FeaturePlot(obj_sub1,features = c('FGF10'))

all_markers <- FindAllMarkers(object = obj_sub,group.by = 'Spatial_snn_res.1',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)



#monocle3
data <- GetAssayData(obj_sub1,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj_sub1@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_gut") + ggtitle('cds.umap')
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj_sub1, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_gut") + ggtitle('int.umap')
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


p + geom_vline(xintercept = seq(2.5,4.5,0.5)) + geom_hline(yintercept = seq(-3.5,0,0.5))
embed <- data.frame(Embeddings(obj_sub1, reduction = "umap"))
embed <- subset(embed, umap_1 > 2.5 & umap_1 < 4.5 & umap_2 > -3.5 & umap_2 < 0)
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

#pdf('E35_E45_lung_main_bronchi_monocle3_top30Genes.pdf',width = 12,height = 8)
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))
#dev.off()

FeaturePlot(obj_sub,features = c('DUSP6','HAPLN1','FAM237A','PTN','AGR2','TNC'))
head(obj_sub)
DimPlot(obj_sub,group.by = 'Stage',label = T)
FeaturePlot(obj_sub,features = c('SOX9'))

library(ClusterGVis)
genes <- row.names(subset(Track_genes, q_value == 0 & morans_I > 0.25))
pre_pseudotime_matrix <- getFromNamespace("pre_pseudotime_matrix","ClusterGVis")
mat <- pre_pseudotime_matrix(cds_obj = cds,
                             gene_list = genes)
ck <- clusterData(obj = mat,
                  cluster.method = "kmeans",
                  cluster.num = 5)
