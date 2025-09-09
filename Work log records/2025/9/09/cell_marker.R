setwd("~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/")
library(Seurat)
library(dplyr)

####E2.5
data <- readRDS("E2.5_foregut.rds")
head(data@meta.data,n=100)
#DefaultAssay(data) <- 'Spatial'
#data <- NormalizeData(data,assay = 'Spatial') %>%
#  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
#  ScaleData() %>%
#  RunPCA(npcs = 30) %>%
#  FindNeighbors(dims=1:30) %>%  # 减少维度
#  FindClusters(resolution=c(2,3)) %>%  # 使用Leiden
#  RunUMAP(dims=1:30)
table(data$celltype_gut)
celltype <- c('Foregut Progenitor Cells1','Foregut Progenitor Cells2','Foregut Progenitor Cells3','Lung Bud1','Lung Bud2',
              'Proventriculus Progenitor Cells','Gizzard Progenitor Cells','Gizzard-Duodenum Progenitor Cells',
              'Hepatocyte Progenitor Cells','Hepatocyte Mesenchymal Progenitor Cells',
              'Foregut Mucosa Epithelial Progenitor Cells',
              'Duodenum Progenitor Cells','Duodenum-Midgut Progenitor Cells',
              'Pancreatic Bud','Hepatopancreatic Ductal Progenitor Cells')
data$celltype_gut <- factor(data$celltype_gut,levels = celltype)
p1 <-DotPlot(data,
        group.by ='celltype_gut',
        features = c('FILIP1','RAB3IP','CREB5','ALDH1A2','AQP1','AQP3',
                     'RIPPLY3','SHH','MME','WNT2','WNT2B','SMOC2','BARX1','NKX3-2',
                     'UTS2B','SLC51B','FGG','FGA','ABCC9','VWDE',
                     'SOX2','CLDN1',
                     'PDX1','RFX6','GUCA2B','SLC7A9','INS','GCG','SOX17','HHEX'))+
  labs(x='',y='')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1
ggsave('E2.5_marker.pdf',p1,dpi = 300,height = 6,width = 14)

#####E3.5
data <- readRDS("E3.5_foregut.rds")
head(data@meta.data,n=100)

table(data$celltype_gut)
celltype <- c('Foregut Mesothelial cells','Lung Primary Bronchi','Lung Secondary Bronchi','Lung Mesenchymal Progenitor Cells',
              'Proventriculus Progenitor Cells','Proventriculus Mucosa Epithelial Cells',
              'Gizzard Progenitor Cells','Gizzard Mucosa Epithelial Cells','Gizzard-Duodenum Progenitor Cells',
              'Hepatocyte Progenitor Cells','Duodenum Progenitor Cells','Duodenum Mucosa Epithelial Cells',
              'Duodenum Enterocyte-Neurons Cells','Duodenum-Midgut Progenitor Cells','Pancreatic Bud',
              'Hindgut Progenitor Cells','Hindgut Mucosa Epithelial Cells','Colacal Epithelial Cells')
data$celltype_gut <- factor(data$celltype_gut,levels = celltype)
p2 <-DotPlot(data,
             group.by ='celltype_gut',
             features = c('TNC','LHX2','NKX2-1','RIPPLY3','ADCY8','MLLT3','ZFPM2','APCDD1L',
                          'WNT2B','SMOC2','CYP2C23a','AGR2','BARX1','NKX3-2','GKN2','NKX6-2',
                          'GPC3','CCDC85A','FGG','FGA','PODXL','AHNAK2','AQP5','TM4SF4',
                          'UTS2B','SLC51B','CD9','LIMCH1','INS','GCG','ANO1','MAB21L1',
                          'SATB2','GUCA2B','HOXD13','EVX1'))+
  labs(x='',y='')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p2
ggsave('E3.5_marker.pdf',p2,dpi = 300,height = 6,width = 14)


#####E4.5
data <- readRDS("E4.5_foregut.rds")
head(data@meta.data,n=100)
table(data$celltype_gut)
celltype <- c('Foregut Mesothelial cells','Lung Primary Bronchi','Lung Secondary Bronchi',
              'Lung Progenitor Cells','Lung Mesenchymal Progenitor Cells','Trachea Mesenchymal Progenitor Cells',
              'Proventriculus Progenitor Cells','Proventriculus Mucosa Epithelial Cells','Gizzard Progenitor Cells',
              'Gizzard Mucosa Epithelial Cells','Gizzard-Duodenum Mesenchymal Progenitor Cells','Gizzard Mesenchymal Progenitor Cells',
              'Hepatocyte Progenitor Cells','Hepatopancreatic Ductal Progenitor Cells','Duodenum Progenitor Cells',
              'Duodenum Mucosa Epithelial Cells','Duodenum Enterocyte-Neurons Cells','Duodenum-Midgut Progenitor Cells',
              'Duodenum Mesenchymal Progenitor Cells','Pancreatic Bud','Midgut-Hindgut Progenitor Cells',
              'Midgut-Hindgut Smooth Muscle Progenitors','Hindgut Progenitor Cells','Hindgut Smooth Muscle Progenitor Cells',
              'Hindgut Mucosa Epithelial Cells','Colacal Epithelial Cells',
              'Spleno-Pancreatic Mesenchyme Progenitor Cells'
              )
data$celltype_gut <- factor(data$celltype_gut,levels = celltype)
p3 <-DotPlot(data,
             group.by ='celltype_gut',
             features = c('TNC','LHX2','NKX2-1','RIPPLY3','MME','WNT4','ADCY8','TBX4','CCDC80','ADCYAP1R1',
                          'JAG1','WNT11','WNT2B','SMOC2','CYP2C23a','AGR2','BARX1','NKX3-2','GKN2','NKX6-2',
                          'GDF6','BMP4','LMO3','CNGA3','FGG','FGA','SOX17','HHEX','CYTL1','PI15',
                          'AQP5','TM4SF4','UTS2B','SLC51B','CD9','ABCC9','LIMCH1','AHNAK2',
                          'INS','GCG','RDH10','HAND1','MYLK','MYL4','HOXD11','HOXA11','ANO1','MAB21L1',
                          'SPINK5','GUCA2B','HOXD13','EVX1','TLX1','GUCY1A2'))+
  labs(x='',y='')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
p3
ggsave('E4.5_marker.pdf',p2,dpi = 300,height =7.5,width = 18)
