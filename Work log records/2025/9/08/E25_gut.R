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
                                           "Outflow Tract","Hepatocytes","Hepatic Mesothelium",
                                           "Pancreatic Primordium","Periportal Mesenchyme",
                                           "Midgut Primordium"))

#data <- readRDS("E25_gut.rds")
head(df@meta.data)
table(df@meta.data$Celltype)
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

DotPlot(data,
        group.by ='Spatial_snn_res.1',
        features = c('BARX1','NKX3-2','ISL1','WNT2B','NKX6-1',
                     'GUCA2B','TTR','HAND1','HAND2',
                     'S100A6','ANXA1','CYP2C23a',
                     'CDX1','CDX2'))

FeaturePlot(df,features = c('CDX1','GUCA2B'),reduction = "umap",raster=FALSE)
head(df@meta.data)
DotPlot(data,features = c('BARX1','NKX2-1','CDX2','PDX1'),group.by ='Spatial_snn_res.3')

FeaturePlot(data,features = c('BARX1','NKX2-8','WNT7B','RIPPLY3'))
colnames(data@meta.data)[colnames(data@meta.data) == "Spatial_snn_res.3"] <- 'celltype_gut'
head(data@meta.data)
saveRDS(data,'E25_gut.rds')


######
library(dplyr)
data <- readRDS("E25_gut.rds")
head(data@meta.data)
data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Proventriculus Progenitor Cells",
  "8" = "Proventriculus Progenitor Cells",
  "9" = "Proventriculus Progenitor Cells",
  "14" = "Proventriculus Progenitor Cells",
  
  "1" = "Foregut Progenitor Cells1",
  
  "2" = "Hepatocyte Progenitor Cells",
  "24" = "Hepatocyte Progenitor Cells",
  
  "3" = "Duodenum-Midgut Progenitor Cells",
  "20" = "Duodenum-Midgut Progenitor Cells",
  
  "4" = "Gizzard Progenitor Cells",
  "11" = "Gizzard Progenitor Cells",
  "15" = "Gizzard Progenitor Cells",
  "18" = "Gizzard Progenitor Cells",
  "19" = "Gizzard Progenitor Cells",
  "22" = "Gizzard Progenitor Cells",
  "23" = "Gizzard Progenitor Cells",
  "25" = "Gizzard Progenitor Cells",
  "28" = "Gizzard Progenitor Cells",
  
  "5" = "Lung Bud1",
  "6" = "Lung Bud1",
  "17" = "Lung Bud1",
  "26" = "Lung Bud2",
  
  "7" = "Foregut Progenitor Cells2",
  "31" = "Foregut Progenitor Cells2",
  "32" = "Foregut Progenitor Cells2",
  
  "10" = "Gizzard-Duodenum Progenitor Cells",
  
  "13" = "Foregut Progenitor Cells3",
  "29" = "Foregut Progenitor Cells3",
  
  "16" = "Hepatocyte Mesenchymal Progenitor Cells",
  
  "21" = "Pancreatic Bud",
  
  "27" = "Duodenum Progenitor Cells",
  "34" = "Duodenum Progenitor Cells",
  
  "30" = "Hepatopancreatic Ductal Progenitor Cells",
  
  "33" = "Foregut Mucosa Epithelial Progenitor Cells"
  
)
head(data@meta.data,n=50)

data1 <- subset(data,subset = celltype_gut %in% c("Proventriculus Progenitor Cells",
                                                  "Foregut Progenitor Cells1",
                                                  "Foregut Progenitor Cells2",
                                                  "Foregut Progenitor Cells3",
                                                  "Hepatocyte Progenitor Cells",
                                                  "Duodenum-Midgut Progenitor Cells",
                                                  "Gizzard Progenitor Cells",
                                                  "Lung Bud1",
                                                  "Lung Bud2",
                                                  "Gizzard-Duodenum Progenitor Cells",
                                                  "Hepatocyte Mesenchymal Progenitor Cells",
                                                  "Pancreatic Bud",
                                                  "Duodenum Progenitor Cells",
                                                  "Hepatopancreatic Ductal Progenitor Cells",
                                                  "Foregut Mucosa Epithelial Progenitor Cells"))                                                                            
head(data1@meta.data)
DimPlot(data1, reduction = "umap", group.by = "celltype_gut",label = T) +theme(legend.position = 'none')
saveRDS(data1,'E2.5_foregut.rds')

data <- readRDS("../E4.5_foregut.rds")
head(data@meta.data)
cellID <- data@meta.data
#colnames(cellID) <- make.unique(colnames(cellID))
selected_data <- cellID %>% 
  mutate(Cell_ID = rownames(cellID)) %>%  # 将行名转为新列
  select(Cell_ID, celltype_gut)                # 选择目标列

write.csv(selected_data,'E4.5_foregut_cellID.csv',row.names = F)

table(data@meta.data$celltype_gut)
cellnumber <- data.frame(table(data@meta.data$celltype_gut))

