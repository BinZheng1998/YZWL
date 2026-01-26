setwd('~/project/03_3D_spatial/02_result/260121_cellchat/')
library(Seurat)
library(CellChat)
library(SeuratData)
library(tidyverse)
packageVersion("CellChat")
options(future.globals.maxSize = 200 * 1024^3)
load('~/project/03_3D_spatial/00_data/database/cellchat/CellChatDB.chicken_jtwen.rda')

df <- readRDS('../../00_data/lung/E45_lung_epithelium_251220.rds')
head(df)
df1 <- readRDS('../../00_data/lung/E45_lung_mesenchymal_251220.rds')
df1$lung_celltype[df1$lung_celltype == "Distal Mesenchymal"] <- "Distal Mesenchyme"
df1$lung_celltype[df1$lung_celltype == "Proximal Mesenchymal"] <- "Proximal Mesenchyme"
df1$celltype <- df1$lung_celltype
head(df1)

sce <- merge(x = df,y = df1, 
                       add.cell.ids = c("Epithelium", "Mesenchyme"), 
                       project = "E45_Lung")

sce <- UpdateSeuratObject(sce)
sce@meta.data$cell <- sce@meta.data$celltype

sce <- sce %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData

Idents(sce) <-"cell"
sce
sce <- JoinLayers(sce,assay = 'Spatial')
data.input = Seurat::GetAssayData(sce, assay = "Spatial", layer = "data")
#data.input = sce@assays$Spatial@layers$data
data.input[1:4,1:4]
meta = sce@meta.data[,"cell",drop=F]# a dataframe with rownames containing cell mata data
colnames(meta) <-"labels"
meta$samples <- 'sample1'
meta$samples <- factor(meta$samples)
meta$labels <- paste0('c',meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by ="labels")
cellchat
levels(cellchat@idents)# show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))# number of cells in each cell group
groupSize

CellChatDB.use <- subsetDB(CellChatDB.chicken)
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat

cellchat <- subsetData(cellchat)# This step is necessary even if using the whole database
future::plan("multisession", workers = 10)# do parallel
# CellChat identifies over-expressed ligands or receptors in one cell group
cellchat <- updateCellChat(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat

cellchat <- computeCommunProb(cellchat,type="triMean")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

df.net <- subsetCommunication(cellchat, thresh = 0.05)
head(df.net)




.

pdf('E45_cellchat_chicken_WNT_CCI_heatmap.pdf',width = 6,height = 5)
pathways.show <-"ncWNT"
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap ="Reds")
dev.off()

pdf('E45_cellchat_chicken_all_CCI_heatmap.pdf',width = 6,height = 5)
netVisual_heatmap(
  cellchat, 
  #signaling = "WNT",
  color.heatmap = "Reds")

dev.off()

pdf('E45_cellchat_chicken_CCI_dotplot.pdf',width = 12,height = 12)
netVisual_bubble(cellchat, 
  sources.use = c("cProximal Mesenchyme","cProximal Epithelium","cDistal Mesenchyme","cDistal Epithelium"), 
  targets.use = c("cProximal Mesenchyme","cProximal Epithelium","cDistal Mesenchyme","cDistal Epithelium"), 
  remove.isolate = FALSE)
dev.off()



saveRDS(cellchat, file = "E45_cellchat_chicken.rds")
