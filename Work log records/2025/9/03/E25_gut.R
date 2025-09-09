setwd(dir = "~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E2.5_new_sp.rds")
#txt_data <- read_excel('../metadata/20250825_SP_E2.5_meta.data.xlsx',sheet = 1)
#txt_data <- txt_data[,-1]
#colnames(txt_data) <- c("Cell_ID", "Tissue", "Celltype")
#missing_cells <- txt_data$Cell_ID[!txt_data$Cell_ID %in% rownames(df@meta.data)]
#if (length(missing_cells) > 0) {
#  warning("以下 Cell ID 在 Seurat 对象中未找到：", paste(missing_cells, collapse = ", "))
#}
#txt_data <- data.frame(Tissue = txt_data$Tissue, Celltype = txt_data$Celltype, 
#                       row.names = txt_data$Cell_ID)
#df@meta.data$Tissue <- txt_data[rownames(df@meta.data), "Tissue"]
#df@meta.data$Celltype <- txt_data[rownames(df@meta.data), "Celltype"]

###
df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
df@meta.data
#data1 <- subset(df,subset = Tissue %in% c("Heart"))
#table(data@meta.data$Celltype)
data <- subset(df,subset = Celltype %in% c("Digestive Gland Progenitors",
                                           "Stomach Mesenchyme","Foregut Primordium",
                                           "Outflow Tract","Hepatocytes",
                                           "Pancreatic Primordium","Periportal Mesenchyme",
                                           "Midgut Primordium"))

data <- readRDS("E25_gut.rds")
head(data@meta.data)
DimPlot(data, reduction = "umap", group.by = "Celltype",raster = F,label = T)+theme(legend.position = 'none')
#DimPlot(df, reduction = "umap", group.by = "Celltype",raster = F,label = T)+theme(legend.position = 'none')
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,4,4.5,5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
p1<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.1",label = T) + ggtitle("Res 0.1")
p2<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.2",label = T) + ggtitle("Res 0.2")
p3<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.5",label = T) + ggtitle("Res 0.5")
p3
p4<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.3",label = T) + ggtitle("Res 3")
p4
DimPlot(data, reduction = "umap", group.by = "Celltype",label = T) 
all_markers <- FindAllMarkers(object = data,group.by = 'Spatial_snn_res.3',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)
#top50_markers <- read.table("../E2.5_gut_res5_top50markers.txt",sep = '\t',header = T)
write.table(top50_markers,file = 'E2.5_gut_res3_top50markers.txt',sep = '\t',row.names = F)

head(data@meta.data)
cellID <- data@meta.data
selected_data <- cellID %>% 
  mutate(Cell_ID = rownames(cellID)) %>%  # 将行名转为新列
  select(Cell_ID, celltype_gut)                # 选择目标列

write.csv(selected_data,'E2.5_cellID.csv',row.names = F)

DotPlot(data,
        group.by ='Spatial_snn_res.1',
        features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                     'GUCA2B','TTR','HAND1','HAND2',
                     'S100A6','ANXA1','CYP2C23a',
                     'CDX1','CDX2'))

FeaturePlot(df,features = c('CDX1','GUCA2B'),reduction = "umap",raster=FALSE)
head(df@meta.data)
DotPlot(df,features = c('GUCA2B','CDX1','CDX2'),group.by ='Celltype')

FeaturePlot(data,features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                              'ANXA1','CYP2C23a','CDX1','CDX2','FGB','FGA','FGG',
                              'GUCA2B','TTR','PDX1','RFX6','HBZ','KRT7','KRT24'))
colnames(data@meta.data)[colnames(data@meta.data) == "Spatial_snn_res.3"] <- 'celltype_gut'
head(data@meta.data)
saveRDS(data,'E25_gut.rds')
