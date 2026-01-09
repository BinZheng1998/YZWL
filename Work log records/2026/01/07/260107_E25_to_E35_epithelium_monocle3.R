setwd('~/project/03_3D_spatial/02_result/260107_E25_to_E25_epithelium_monocle3/')
library(Seurat)
library(ggunchull)
library(dplyr)
options(future.globals.maxSize = 50 * 1024^3)
df1 <- readRDS('../../00_data/lung/E25_lung_bronchi.rds')
df1$celltype <- 'E25 Lung Epithelium'
df1$Stage <- 'E25'
df2 <- readRDS('../../00_data/lung/E35_lung_epithelium_251222.rds')
df2$celltype <- 'E35 Lung Epithelium'
df2$Stage <- 'E35'
df2 <- NormalizeData(df2,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(df2,group.by = 'Spatial_snn_res.0.5',label = T,raster=FALSE)+NoLegend()
FeaturePlot(df2,features = 'AGR2')
df2$celltype[df2$Spatial_snn_res.0.5 == '0'] <- 'Proximal Epithelium'
df2$celltype[df2$Spatial_snn_res.0.5 == '1'] <- 'Distal Epithelium'

head(df2)

obj <- merge(x = df1, 
            y = df2, 
            add.cell.ids = c("E25_Epi","E35_Epi"), 
            project = "Chicken_Lung_Dev")
DefaultAssay(obj) <- 'Spatial'

obj_sub<-obj
#obj_sub<- UpdateSeuratObject(obj_sub)
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.rpca",k.anchor= 1,k.weight=10,
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)

DimPlot(obj_sub,group.by = 'celltype')+ggtitle('k.anchor 1 k.weight 10')

obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",theta = 0.75,lambda = 1,
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.harmony", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
p1<-DimPlot(obj_sub,group.by = 'celltype',label = T,label.size = 5)+ 
    labs(x = '',y='UMAP2',title = '')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size=15))+
  NoLegend()
p1
FeaturePlot(obj_sub,features = "FABP3")


###monocle3
#monocle3
library(monocle3)
data <- GetAssayData(obj_sub,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('cds.umap')


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('int.umap')


cds <- cluster_cells(cds,cluster_method = 'louvain')
cds <- learn_graph(cds)
p<- plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
p
p + geom_vline(xintercept = seq(-4.5,-2,0.5)) + geom_hline(yintercept = seq(-2,1,0.5))
embed <- data.frame(Embeddings(obj_sub, reduction = "umap"))
embed <- subset(embed, umap_1 > -4.5 & umap_1 < -2 & umap_2 > -2 & umap_2 < 1)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
#pdf('E25_E35_lung_epithelium_monocle3_pseudotime_umap.pdf',width = 7, height = 5)
p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 2,
          trajectory_graph_segment_size = 1)+ 
    #labs(x = 'UMAP1',y='UMAP2')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  NoLegend()
p2
#dev.off()

p1+p2

cds_sub1 <- choose_graph_segments(cds,clear_cds = F,return_list = T)
cds_sub1 <- cds[,cds_sub1$cells]

cds_sub2 <- choose_graph_segments(cds,clear_cds = F,return_list = T)
cds_sub2 <- cds[,cds_sub2$cells]


p3<-plot_cells(cds_sub1, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
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
        axis.title = element_text(size=15))+
  NoLegend()
p3

p4<-plot_cells(cds_sub2, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 1,
          trajectory_graph_segment_size = 1)+ 
  labs(x = 'UMAP1',y='')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size=15))+
  NoLegend()
p4

#pdf('E25_E35_lung_epithelium_monocle3_pseudotime_umap.pdf',width = 10, height = 8)
p1+p2+p3+p4
#dev.off()



Track_genes1 <- graph_test(cds_sub1, neighbor_graph="principal_graph", cores=20)
#write.table(Track_genes,'E45_lung_mesenchymal_monocle3_res.txt',sep = '\t')
Track_genes_sig1 <- Track_genes1 %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig1 <- c('TNC','WNT9A','APOA2','PRTG')
p5<-plot_genes_in_pseudotime(cds_sub1[Track_genes_sig1,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
    labs(x='')+
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        strip.text = element_text(face = "italic", size = 12),
        axis.title = element_text(size=15))+NoLegend()
p5

Track_genes2 <- graph_test(cds_sub2, neighbor_graph="principal_graph", cores=20)
#write.table(Track_genes,'E45_lung_mesenchymal_monocle3_res.txt',sep = '\t')
Track_genes_sig2 <- Track_genes2 %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig2 <- c('AGR2','ANXA8L1','NPM3','HSP90AA1')
p6<-plot_genes_in_pseudotime(cds_sub2[Track_genes_sig2,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
  theme_classic()+
  labs(x='Pseudotime')+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        strip.text = element_text(face = "italic", size = 12),
        axis.title = element_text(size=15))+NoLegend()
p6

p <- p1+p2+p3+p4
p

library(patchwork)
pdf('E25_E35_lung_epithelium_monocle3_pseudotime_umap.pdf',width = 10, height = 14)
p / (p5 / p6) + plot_layout(heights = c(1.75, 1.25))
dev.off()
