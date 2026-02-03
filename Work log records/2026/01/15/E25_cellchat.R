setwd('~/project/03_3D_spatial/02_result/260114_E25_cellchat/')
library(Seurat)
library(CellChat)
library(SeuratData)#加载seurat数据集
library(tidyverse)
packageVersion("CellChat")
load('CellChatDB.chicken.rda')

df <- readRDS('../../00_data/lung/E25_lung_bronchi.rds')
df$lung_celltype <- 'E25 Lung Epithelium'

df1 <- readRDS('../../00_data/lung/E25_lung_mesenchymal_1216.rds')
head(df1)

sce <- merge(x = df, 
                       y = df1, 
                       add.cell.ids = c("Epithelium", "Mesenchyme"), 
                       project = "E25_Lung")
sce <- UpdateSeuratObject(sce)
sce@meta.data$cell <- sce@meta.data$lung_celltype

sce <- sce %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData

Idents(sce) <-"cell"

data.input = sce@assays$Spatial@data# normalized data matrix
data.input[1:4,1:4]
meta = sce@meta.data[,"cell",drop=F]# a dataframe with rownames containing cell mata data
colnames(meta) <-"labels"
meta$samples <- 'sample1'
meta$samples <- factor(meta$samples)
meta$labels <- paste0('c',meta$labels)
head(meta)
table(meta)

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


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pdf('cellchat_chicken_WNT_CCI_heatmap.pdf',width = 6,height = 5)
pathways.show <-"ncWNT"
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap ="Reds")
dev.off()

pdf('cellchat_chicken_all_CCI_heatmap.pdf',width = 6,height = 5)
netVisual_heatmap(
  cellchat, 
  #signaling = "WNT",
  color.heatmap = "Reds")

dev.off()

pdf('cellchat_chicken_CCI_dotplot.pdf',width = 8,height = 12)
netVisual_bubble(cellchat, 
  sources.use = c("cWNT2+ MSCs","cE25 Lung Epithelium","cMSCs"), 
  targets.use = c("cWNT2+ MSCs","cE25 Lung Epithelium","cMSCs"), 
  remove.isolate = FALSE)
dev.off()

saveRDS(cellchat, file = "cellchat_chicken.rds")
