setwd('~/project/03_3D_spatial/02_result/251203_E35_lung_epithliem_replot/')
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

data <- readRDS("~/project/03_3D_spatial/00_data/lung/E35_lung_bronchi.rds")
DimPlot(data)
obj <- data
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2500,)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj)
obj$lung_celltype <- 'E35 Lung Epithliem'

pdf('E35_lung_epithliem_umap1.pdf',width = 7,height = 6)
DimPlot(obj,group.by = 'lung_celltype',cols = 'green')+
  ggtitle('')+
  labs(x = 'UMAP1',y='UMAP2')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size =20))+
  NoLegend()
dev.off()

pdf('E35_lung_epithliem_genes.pdf',width = 10,height = 8)
FeaturePlot(obj,features = c('NKX2-1','FGFR1','FGFR2','FGFR3','SOX9','ETV4','ETV5','BMP4','SOX2','SHH','ID2','SPRY2'),
coord.fixed = T,order = T) &
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = "italic")) &
  NoLegend()
dev.off()



#monocle3
data <- GetAssayData(obj,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj@meta.data
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
int.embed <- Embeddings(obj, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype_gut") + ggtitle('int.umap')
p2

cds <- cluster_cells(cds,cluster_method = 'louvain')
cds <- learn_graph(cds)
p = plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p

p + geom_vline(xintercept = seq(1,4,0.5)) + geom_hline(yintercept = seq(-3,-0.5,0.5))
embed <- data.frame(Embeddings(obj, reduction = "umap"))
embed <- subset(embed, umap_1 > 1 & umap_1 < 4 & umap_2 > -3 & umap_2 < -0.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
pdf('E35_lung_epithliem_monocle3_pseudotime_umap.pdf',width = 8, height = 6)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 1,
          trajectory_graph_segment_size = 1)+ 
  labs(x = 'UMAP1',y='UMAP2')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.title = element_text(size=20))+
  scale_color_viridis_c(name = "Pseudotime", option = "plasma")
dev.off()

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=20)
write.table(Track_genes,'E35_lung_epithliem_monocle3_res.txt',sep = '\t',row.names = F)
Track_genes_sig <- Track_genes %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()


pdf('E35_lung_epithliem_monocle3_top30Genes.pdf',width = 16,height = 12)
plot_genes_in_pseudotime(cds[Track_genes_sig,],
  color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
  labs(x='Pseudotime')+
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15, face = "italic"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        plot.title = element_text(face = "italic"))+
  guides(color = guide_legend(override.aes = list(size = 5)))
dev.off()
