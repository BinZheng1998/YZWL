setwd('~/project/03_3D_spatial/02_result/251122_E25_lung_mesenchymal/')
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

data <- readRDS("~/project/03_3D_spatial/04_cell_annotation/new_data_2509/E2.5_foregut.rds")
table(data$Spatial_snn_res.4)
DimPlot(data,group.by = 'Spatial_snn_res.4',label = T,label.size = 8)
FeaturePlot(data,features = 'NKX2-1')
head(data)
df <- readRDS("~/project/03_3D_spatial/00_data/E2.5_SP_new_coord_3d_2511.rds")
metadata <- data.frame(df$cell,df$trans_x,df$trans_y,df$trans_z)
#rownames(metadata)
colnames(metadata) <- c("cell","trans_x","trans_y","trans_z")
data$trans_x <- NULL
data$trans_y <- NULL
data$trans_z <- NULL
data <- AddMetaData(data,metadata = metadata)

all_markers <- FindAllMarkers(object = data,group.by = 'Spatial_snn_res.4',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

source("~/project/03_3D_spatial/03-script/export_clusters_to_ply.R")
#data <- subset(obj,subset = Stage == "E45")
color_mapping <- export_clusters_to_ply(data,
                                        output_dir = "./E25_ply/",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.4")

obj <- subset(data,subset = Spatial_snn_res.4 %in% c("3","12","16","38","35","1","4","6","7","9","27","17","23"))
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = 2000,)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = c(1,2,3,4,5))
obj <- RunUMAP(obj, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)

DimPlot(obj,label = T,group.by = "Spatial_snn_res.3",label.size = 8)
FeaturePlot(obj,features = c("NKX2-1","BARX1","FOXE3","WNT2B"))
all_markers <- FindAllMarkers(object = obj,group.by = 'Spatial_snn_res.3',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

color_mapping <- export_clusters_to_ply(obj,
                                        output_dir = "./E25_ply/2",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.3")


obj1 <- subset(obj,subset = Spatial_snn_res.3 %in% c('11','6','8','13','12','16','0','1','4','7','15','17'))
obj1 <- NormalizeData(obj1)
obj1 <- FindVariableFeatures(obj1,selection.method = "vst",nfeatures = 2000,)
obj1 <- ScaleData(obj1)
obj1 <- RunPCA(obj1)
obj1 <- FindNeighbors(obj1, dims = 1:30, reduction = "pca")
obj1 <- FindClusters(obj1, resolution = c(1,2,3,4,5))
obj1 <- RunUMAP(obj1, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj1,label = T,group.by = "Spatial_snn_res.2",label.size = 8)
FeaturePlot(obj1,features = c('WNT2B'))
table(obj1$Spatial_snn_res.2)
 
color_mapping <- export_clusters_to_ply(obj1,
                                        output_dir = "./E25_ply/test",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.2")


obj2 <- subset(obj1,subset = Spatial_snn_res.2 %in% c('0','5','7'))
saveRDS(obj2,'../../00_data/lung/E25_lung_bronchi.rds')
obj2 <- NormalizeData(obj2)
obj2 <- FindVariableFeatures(obj2,selection.method = "vst",nfeatures = 2000,)
obj2 <- ScaleData(obj2)
obj2 <- RunPCA(obj2)
obj2 <- FindNeighbors(obj2, dims = 1:30, reduction = "pca")
obj2 <- FindClusters(obj2, resolution = c(1,1.2,1.3,1.5,3,4,5))
obj2 <- RunUMAP(obj2, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj2,label = T,group.by = "Spatial_snn_res.1.2",label.size = 8)
FeaturePlot(obj2,features = 'BARX1')
all_markers <- FindAllMarkers(object = obj2,group.by = 'Spatial_snn_res.1.2',
                              only.pos = TRUE,
                              min.pct = 0.15,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
#write.csv(obj2$cell,'../../00_data/lung/E25_lung_bronchi_cellID.csv')

color_mapping <- export_clusters_to_ply(obj2,
                                        output_dir = "./E25_ply/4",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.1.5")

data$lung <- 'lung'
data$lung[data$cell %in% cellid] <- 'E25'
DimPlot(data,group.by = 'lung',label = T,label.size = 8)
DimPlot(data,group.by = 'celltype_gut',label = T,label.size = 5)
FeaturePlot(data,features = 'RIPPLY3')

obj3 <- subset(obj1,subset = Spatial_snn_res.2 %in% c('6','9'))
df <- readRDS('../../00_data/lung/E25_lung_mesenchymal.rds')
library(Seurat)
DimPlot(df)
saveRDS(obj3,'../../00_data/lung/E25_lung_mesenchymal.rds')
#table(obj3$Spatial_snn_res.2)
write.csv(df$cell,'../../00_data/lung/E25_lung_mesenchymal_cellID.csv')

obj3 <- subset(obj1,subset = Spatial_snn_res.2 %in% c('1','2','4','6','8','9'))
obj3 <- NormalizeData(obj3)
obj3 <- FindVariableFeatures(obj3,selection.method = "vst",nfeatures = 2000,)
obj3 <- ScaleData(obj3)
obj3 <- RunPCA(obj3)
obj3 <- FindNeighbors(obj3, dims = 1:30, reduction = "pca")
obj3 <- FindClusters(obj3, resolution = c(1,1.2,1.3,1.5,2,3,4,5))
obj3 <- RunUMAP(obj3, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj3,label = T,group.by = "Spatial_snn_res.1",label.size = 8)
FeaturePlot(obj3,features = c('FGF10','BMP4','WNT2','WNT2B','BARX1'))
FeaturePlot(obj3,features = c('WNT2'))
all_markers <- FindAllMarkers(object = obj3,group.by = 'Spatial_snn_res.1',
                              only.pos = TRUE,
                              min.pct = 0.15,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
color_mapping <- export_clusters_to_ply(obj3,
                                        output_dir = "./E25_ply/5",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.1")

obj3$lung_celltype <- 'MSCs'
obj3$lung_celltype[obj3$Spatial_snn_res.1 == "2"] <- 'WNT2+ MSCs'

DimPlot(obj3,group.by = 'lung_celltype')
pdf('E25_lung_mesenchymal_umap1.pdf',width = 7,height = 6)
DimPlot(obj3,group.by = 'lung_celltype',label = T,label.size = 8)+
  ggtitle('')+
  labs(x = 'UMAP1',y='UMAP2')+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size =20))+
  NoLegend()
dev.off()

all_markers <- FindAllMarkers(object = obj3,group.by = 'lung_celltype',
                              #only.pos = TRUE,
                              min.pct = 0.15,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

FeaturePlot(obj3,features = 'GATA4')
library(scRNAtoolVis)
jjVolcano(diffData = all_markers,flip = T,legend.position = c(-4,2))
pdf('E25_MSCs_volcano.pdf')
markerVolcano(markers = all_markers,
              topn = 7,
              labelCol = ggsci::pal_npg()(9))
dev.off()
