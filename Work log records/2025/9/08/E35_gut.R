setwd(dir = "~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/")
library(Seurat)
library(dplyr)
library(ggplot2)

#df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E3.5_new_sp.rds")
#df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
#df@meta.data
#data <- subset(df,subset = Celltype %in% c("Midgut Mesenchyme","Midgut Epithelium",
#                                           "Hindgut Epithelium","Hindgut Mesenchyme",
#                                           "Foregut Epithelium",
#                                            #"Peritoneal Mesothelium",
 #                                          "Hepatocytes","Hepatic Mesothelium","Hepatic Primordium",
#                                           "Cloacal Epithelium",
 #                                           "Gizzard Mesenchyme","Proventriculus Mesenchyme",
 #                                          "Lung Epithelium","Lung Mesenchyme",
 #                                           "Pancreatic Primordium"
#                                            ))

data <- readRDS("E35_gut.rds")
head(data@meta.data)
DimPlot(data, reduction = "umap", group.by = "Celltype",raster = F,)
DimPlot(df, reduction = "umap", group.by = "Celltype",raster = F,)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(2,3,4,5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
p1<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.1",label = T) + ggtitle("Res 0.1")
p2<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.2",label = T) + ggtitle("Res 0.2")
p3<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.5",label = T) + ggtitle("Res 0.5")
p3
p4<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.5",label = T) + ggtitle("Res 3")
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
write.table(top50_markers,file = 'E3.5_gut_res3_top50markers.txt',sep = '\t',row.names = F)

DotPlot(data,
        group.by ='Spatial_snn_res.1.5',
        features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                     'GUCA2B','TTR','HAND1','HAND2',
                     'S100A6','ANXA1','CYP2C23a',
                     'CDX1','CDX2'))

FeaturePlot(data,features = c('FGG','BARX1','NKX2-1','WNT7B'))

FeaturePlot(data,features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                              'ANXA1','CYP2C23a','CDX1','CDX2','FGB','FGA','FGG',
                              'GUCA2B','TTR','PDX1','RFX6','HBZ','KRT7','KRT24'))
colnames(data@meta.data)[colnames(data@meta.data) == "Spatial_snn_res.3"] <- 'celltype_gut'
colnames(data@meta.data)
saveRDS(data,'E35_gut.rds')

####
data <- readRDS("E35_gut.rds")
head(data@meta.data)
data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Duodenum-Midgut Progenitor Cells",
  
  "1" = "Hindgut Progenitor Cells",
  "8" = "Hindgut Progenitor Cells",
  "29" = "Hindgut Progenitor Cells",
  
  "3" = "Gizzard Progenitor Cells",
  "11" = "Gizzard Progenitor Cells",
  "16" = "Gizzard Progenitor Cells",
  "21" = "Gizzard Progenitor Cells",
  "26" = "Gizzard Progenitor Cells",
  
  "4" = "Proventriculus Progenitor Cells",
  "5" = "Proventriculus Progenitor Cells",
  "18" = "Proventriculus Progenitor Cells",
  "20" = "Proventriculus Progenitor Cells",
  
  "6" = "Duodenum Enterocyte-Neurons Cells",
  "15" = "Duodenum Enterocyte-Neurons Cells",
  
  "7" = "Hepatocyte Progenitor Cells",
  "10" = "Hepatocyte Progenitor Cells",
  "24" = "Hepatocyte Progenitor Cells",
  "32" = "Hepatocyte Progenitor Cells",
  "34" = "Hepatocyte Progenitor Cells",
  "38" = "Hepatocyte Progenitor Cells",
  "39" = "Hepatocyte Progenitor Cells",
  
  "9" = "Gizzard-Duodenum Progenitor Cells",
  "12" = "Gizzard-Duodenum Progenitor Cells",
  
  "13" = "Foregut Mesothelial cells",
  "25" = "Foregut Mesothelial cells",
  "28" = "Foregut Mesothelial cells",
  "36" = "Foregut Mesothelial cells",

  "14" = "Duodenum Progenitor Cells",
  "17" = "Duodenum Progenitor Cells",
  
  "19" = "Hindgut Mucosa Epithelial Cells",
  
  "22" = "Lung Mesenchymal Progenitor Cells",
  "27" = "Lung Mesenchymal Progenitor Cells",
  "23" = "Lung Secondary Bronchi",
  "33" = "Lung Primary Bronchi",
  
  "30" = "Colacal Epithelial Cells",
  
  "35" = "Proventriculus Mucosa Epithelial Cells",
  "37" = "Gizzard Mucosa Epithelial Cells",
  "40" = "Duodenum Mucosa Epithelial Cells",
  
  "41" = "Pancreatic Bud",
  "42" = "Pancreatic Bud"

)
head(data@meta.data,n=50)

data1 <- subset(data,subset = celltype_gut %in% c("Duodenum-Midgut Progenitor Cells",
                                                  "Hindgut Progenitor Cells",
                                                  "Gizzard Progenitor Cells",
                                                  "Proventriculus Progenitor Cells",
                                                  "Duodenum Enterocyte-Neurons Cells",
                                                  "Hepatocyte Progenitor Cells",
                                                  "Gizzard-Duodenum Progenitor Cells",
                                                  "Foregut Mesothelial cells",
                                                  "Duodenum Progenitor Cells",
                                                  "Hindgut Mucosa Epithelial Cells",
                                                  "Lung Mesenchymal Progenitor Cells",
                                                  "Lung Secondary Bronchi",
                                                  "Lung Primary Bronchi",
                                                  "Colacal Epithelial Cells",
                                                  "Proventriculus Mucosa Epithelial Cells",
                                                  "Gizzard Mucosa Epithelial Cells",
                                                  "Duodenum Mucosa Epithelial Cells",
                                                  "Pancreatic Bud"))                                                                            
head(data1@meta.data)
DimPlot(data1, reduction = "umap", group.by = "celltype_gut",label = T) 
saveRDS(data1,'E3.5_foregut.rds')
