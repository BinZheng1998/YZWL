setwd('~/project/03_3D_spatial/04_cell_annotation/new_data_2509/lung')
library(monocle3)
library(Seurat)
library(tidyverse)
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
options(future.globals.maxSize = 100* 1e10)
obj.list = list()
stage_list <- c('E2.5','E3.5','E4.5')
library(Seurat)

for (i in 1:length(stage_list)) {
  obj.list[[i]] <- readRDS(paste0('../../new_data_20250827/', stage_list[i], '_foregut.rds'))
  obj.list[[i]] <- UpdateSeuratObject(obj.list[[i]])
  obj.list[[i]]$Stage <- stage_list[i]
}

obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])


obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                   'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Bud1','Lung Primary Bronchi'))
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                           new.reduction = "integrated.rpca",
                           k.weight = 50,k.anchor = 1,
                           verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))

FeaturePlot(obj_sub,features = c('TNC','WNT9A'))
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


## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method = 'louvain')
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

## 识别轨迹
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p

p + geom_vline(xintercept = seq(-3,0,0.5)) + geom_hline(yintercept = seq(-4.5,-2.5,0.5))
embed <- data.frame(Embeddings(obj_sub, reduction = "umap"))
embed <- subset(embed, umap_1 > -3 & umap_1 < 0 & umap_2 > -4.5 & umap_2 < -2.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1)


Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=20)
Track_genes_sig <- Track_genes %>% top_n(n=20, morans_I) %>%
  pull(gene_short_name) %>% as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,],
                         min_expr=0.5, ncol = 4)


genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                             cell_group=colData(cds)$predicted.id)
head(cds)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
