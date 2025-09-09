setwd(dir = "~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/")
library(Seurat)
library(dplyr)
library(ggplot2)

df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E4.5_new_sp.rds")
#df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
df@meta.data
data <- subset(df,subset = Celltype %in% c("Gut Smooth Muscle Progenitors","Gastric Epithelium",
                                            "Midgut Smooth Muscle Progenitors","Midgut Epithelium","Midgut Mesenchyme",
                                            #"Gonadal Ridge Mesenchyme","Gonadal Primordium",
                                           "Proventriculus Mesenchyme","Gizzard Mesenchyme",
                                           "Cloacal Epithelium",
                                           "Hepatic Mesenchyme","Hepatic Endoderm","Hepatocytes",
                                            "Foregut Mesenchyme","Foregut Epithelium",
                                           "Hindgut Smooth Muscle Progenitors","Hindgut Epithelium",
                                            "Pancreatic Primordium",
                                           "Lung Epithelium","Lung Mesenchyme"))

data <- readRDS('E45_gut.rds')
head(data@meta.data)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.2,0.3,0.5,0.7,1,1.5,2,3)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
p1<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.1",label = T) + ggtitle("Res 0.1")
p2<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.2",label = T) + ggtitle("Res 0.2")
p3<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.5",label = T) + ggtitle("Res 0.5")
p3
p4<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.3",label = T) + ggtitle("Res 3")
p4
DimPlot(data, reduction = "umap", group.by = "Celltype",label = T) +theme(legend.position = 'none')
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
write.table(all_markers,file = 'E4.5_gut_res3_allmarkers.txt',sep = '\t',row.names = F)
write.table(top50_markers,file = 'E4.5_gut_res3_top50markers.txt',sep = '\t',row.names = F)

head(data@meta.data)
cellID <- data@meta.data
selected_data <- cellID %>% 
  mutate(Cell_ID = rownames(cellID)) %>%  # 将行名转为新列
  select(Cell_ID, celltype_gut)                # 选择目标列

write.csv(selected_data,'E4.5_cellID.csv',row.names = F)
DotPlot(data,
        group.by ='Spatial_snn_res.1.5',
        features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                     'GUCA2B','TTR','HAND1','HAND2',
                     'S100A6','ANXA1','CYP2C23a',
                     'CDX1','CDX2'))

FeaturePlot(data,features = c('ANXA1','GCG','PDX1','INS'))

FeaturePlot(data,features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                              'ANXA1','CYP2C23a','CDX1','CDX2','FGB','FGA','FGG',
                              'GUCA2B','TTR','PDX1','RFX6','HBZ','KRT7','KRT24'))
head(data@meta.data)
colnames(data@meta.data)[colnames(data@meta.data) == "Spatial_snn_res.3"] <- 'celltype_gut'
saveRDS(data,'E45_gut.rds')

#####################
data <- readRDS('E45_gut.rds')
head(data@meta.data)
data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Hepatocyte Progenitor Cells",
  "2" = "Hepatocyte Progenitor Cells",
  "5" = "Hepatocyte Progenitor Cells",
  "31" = "Hepatocyte Progenitor Cells",
  "40" = "Hepatocyte Progenitor Cells",
  
  "1" = "Gizzard Progenitor Cells",
  "21" = "Gizzard Progenitor Cells",
  "23" = "Gizzard Progenitor Cells",
  "41" = "Gizzard Progenitor Cells", 
  
  "3" = "Gizzard-Duodenum Mesenchymal Progenitor Cells",
  
  "4" = "Trachea Mesenchymal Progenitor Cells", 
  
  "6" = "Lung Mesenchymal Progenitor Cells",
  "49" = "Lung Mesenchymal Progenitor Cells",
  
  "7" = "Midgut-Hindgut Progenitor Cells",
  "17" = "Midgut-Hindgut Progenitor Cells",
  "43" = "Midgut-Hindgut Smooth Muscle Progenitors",
  
  "8" = "Foregut Mesothelial cells",
  "12" = "Foregut Mesothelial cells",
  "25" = "Foregut Mesothelial cells",
  "29" = "Foregut Mesothelial cells",
  
  "9" = "Proventriculus Progenitor Cells",
  "15" = "Proventriculus Progenitor Cells",
  "30" = "Proventriculus Progenitor Cells",
  "32" = "Proventriculus Progenitor Cells",
  "42" = "Proventriculus Progenitor Cells",
  
  "10" = "Lung Progenitor Cells",
  "20" = "Lung Secondary Bronchi",
  "33" = "Lung Primary Bronchi",
  
  "11" = "Gizzard Mesenchymal Progenitor Cells",
  "19" = "Gizzard Mesenchymal Progenitor Cells",
  "44" = "Gizzard Mesenchymal Progenitor Cells",
  "50" = "Gizzard Mesenchymal Progenitor Cells",

  "13" = "Duodenum-Midgut Progenitor Cells",

  "14" = "Colacal Epithelial Cells",
  
  "16" = "Hindgut Progenitor Cells",
  "27" = "Hindgut Progenitor Cells",
  "36" = "Hindgut Progenitor Cells",
  "18" = "Hindgut Smooth Muscle Progenitor Cells",
  
  "22" = "Duodenum Mesenchymal Progenitor Cells",
  "35" = "Duodenum Progenitor Cells",
  "45" = "Duodenum Progenitor Cells",
  "46" = "Duodenum Progenitor Cells",
  "47" = "Duodenum Progenitor Cells",
  
  "24" = "Duodenum Enterocyte-Neurons Cells",
  "38" = "Duodenum Enterocyte-Neurons Cells",
  
  "26" = "Hindgut Mucosa Epithelial Cells",
  "28" = "Gizzard Mucosa Epithelial Cells",
  "34" = "Proventriculus Mucosa Epithelial Cells",
  "48" = "Duodenum Mucosa Epithelial Cells",
  
  "37" = "Pancreatic Bud",
  "52" = "Pancreatic Bud",
  
  "39" = "Spleno-Pancreatic Mesenchyme Progenitor Cells",

  "51" = "Hepatopancreatic Ductal Progenitor Cells"
)
head(data@meta.data,n=50)
table(data$celltype_gut)
data1 <- subset(data,subset = celltype_gut %in% c("Hepatocyte Progenitor Cells",
                                                  "Gizzard Progenitor Cells",
                                                  "Gizzard-Duodenum Mesenchymal Progenitor Cells",
                                                  "Trachea Mesenchymal Progenitor Cells",
                                                  "Lung Mesenchymal Progenitor Cells",
                                                  "Midgut-Hindgut Progenitor Cells",
                                                  "Midgut-Hindgut Smooth Muscle Progenitors",
                                                  "Foregut Mesothelial cells",
                                                  "Proventriculus Progenitor Cells",
                                                  "Lung Progenitor Cells",
                                                  "Lung Secondary Bronchi",
                                                  "Lung Primary Bronchi",
                                                  "Gizzard Mesenchymal Progenitor Cells",
                                                  "Duodenum-Midgut Progenitor Cells",
                                                  "Colacal Epithelial Cells",
                                                  "Hindgut Progenitor Cells",
                                                  "Hindgut Smooth Muscle Progenitor Cells",
                                                  "Duodenum Mesenchymal Progenitor Cells",
                                                  "Duodenum Progenitor Cells",
                                                  "Duodenum Enterocyte-Neurons Cells",
                                                  "Hindgut Mucosa Epithelial Cells",
                                                  "Gizzard Mucosa Epithelial Cells",
                                                  "Proventriculus Mucosa Epithelial Cells",
                                                  "Duodenum Mucosa Epithelial Cells",
                                                  "Pancreatic Bud",
                                                  "Spleno-Pancreatic Mesenchyme Progenitor Cells",
                                                  "Hepatopancreatic Ductal Progenitor Cells"

                                                  ))                                                                            
head(data1@meta.data)
DimPlot(data1, reduction = "umap", group.by = "celltype_gut",label = T) +theme(legend.position = 'none')
saveRDS(data1,'E4.5_foregut.rds')
