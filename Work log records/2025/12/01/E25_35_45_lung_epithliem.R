setwd('../251201_E25_35_45_lung_epithliem_monocle3/')
library(Seurat)
library(ggplot2)
library(dplyr)
library(monocle3)
options(future.globals.maxSize = 100 * 1024^3)
obj.list = list()
stage_list <- c('E25','E35','E45')
for (i in 1:length(stage_list)){
  obj.list[[i]] <- readRDS(paste0('../../00_data/lung/',stage_list[i],'_lung_bronchi.rds'))
  obj.list[[i]]$Stage <- stage_list[i]
}
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
DimPlot(obj,group.by = 'Stage')
obj_sub<-obj
#obj_sub<- UpdateSeuratObject(obj_sub)
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.rpca",k.anchor= 0.25,k.weight=15,
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub,group.by = 'Stage')+ggtitle('k.anchor 0.25 k.weight 15')
FeaturePlot(obj_sub,features = 'NKX2-1')

pdf('E25_35_45_lung_epithliem_umap1.pdf',width = 7,height = 6)
DimPlot(obj_sub,group.by = 'Stage')+
  ggtitle('')+
  labs(x = 'UMAP1',y='UMAP2')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size =20))
dev.off()

#monocle3
data <- GetAssayData(obj_sub,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj_sub@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_gut") + ggtitle('cds.umap')
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_gut") + ggtitle('int.umap')
p2

cds <- cluster_cells(cds,cluster_method = 'louvain')
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p

p + geom_vline(xintercept = seq(-3,-1,0.5)) + geom_hline(yintercept = seq(-3.5,-1.5,0.5))
embed <- data.frame(Embeddings(obj_sub, reduction = "umap"))
embed <- subset(embed, umap_1 > -3 & umap_1 < -1 & umap_2 > -3.5 & umap_2 < -1.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
pdf('E25_35_45_lung_epithliem_monocle3_pseudotime_umap.pdf',width = 8, height = 6)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 1,
          trajectory_graph_segment_size = 1)+ 
  labs(x = 'UMAP1',y='UMAP2')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.title = element_text(size=20))+
  scale_color_viridis_c(name = "Pseudotime", option = "plasma")
dev.off()

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=20)
Track_genes_sig <- Track_genes %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()

pdf('E25_35_45_lung_epithliem_monocle3_top30Genes.pdf',width = 16,height = 12)
plot_genes_in_pseudotime(cds[Track_genes_sig,],
  color_cells_by = "Stage",min_expr=0.5, ncol = 5)+ 
  labs(x='Pseudotime')+
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "italic"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()

FeaturePlot(obj_sub,features = c('DUSP6','HAPLN1','FAM237A','PTN','AGR2','TNC','SHH','NKX2-1','NKX2-8'))


library(ClusterGVis)
library(ComplexHeatmap)
library(dplyr)
library(org.Gg.eg.db)
library(clusterProfiler)

modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 24)
genes <- row.names(subset(modulated_genes, q_value < 0.01 & morans_I > 0.1))
pre_pseudotime_matrix <- getFromNamespace("pre_pseudotime_matrix","ClusterGVis")
mat <- pre_pseudotime_matrix(cds_obj = cds,gene_list = genes)
head(mat[1:5,1:5])
ck <- clusterData(obj = mat,
                  clusterMethod = "kmeans",
                  clusterNum = 6)
ck$cluster.list

enrich_bp <- enrichCluster(object = ck,
                       OrgDb = org.Gg.eg.db,
                       type = "BP",
                       idTrans = T,
                       #fromType = "SYMBOL",
                       #toType = c("ENTREZID"),
                       readable = T,
                       pvalueCutoff = 0.5,
                       topn = 10)
enrich_cc <- enrichCluster(object = ck,
                       OrgDb = org.Gg.eg.db,
                       type = "CC",
                       idTrans = T,
                       #fromType = "SYMBOL",
                       #toType = c("ENTREZID"),
                       readable = T,
                       pvalueCutoff = 0.5,
                       topn = 10)
enrich_mf <- enrichCluster(object = ck,
                       OrgDb = org.Gg.eg.db,
                       type = "MF",
                       idTrans = T,
                       #fromType = "SYMBOL",
                       #toType = c("ENTREZID"),
                       readable = T,
                       pvalueCutoff = 0.5,
                       topn = 10)
enrich <- rbind(enrich_bp,enrich_cc,enrich_mf)


pdf('E25_35_45_lung_epithliem_cluster2.pdf',width = 8,height = 12)
visCluster(
  object = ck,border = F,
  addSampleAnno = F,
  show_row_dend = F,
  #annoTermData = enrich, #软件问题
  #markGenes = sample(rownames(mat),80),
  plotType = "both"
)
dev.off()


library(GO.db)
eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
