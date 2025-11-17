setwd('~/project/03_3D_spatial/02_result/251101_lung_hdWGCNA//')
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(patchwork)
library(stringr)
library(rjson)
library(data.table)
library(harmony)
library(randomcoloR)
library(WGCNA)
library(hdWGCNA)
library(scplotter)
enableWGCNAThreads(nThreads = 48)
set.seed(12345)
options(future.globals.maxSize = 100* 1e10)
obj.list = list()
stage_list <- c('E3.5','E4.5')

for (i in 1:length(stage_list)) {
  obj.list[[i]] <- readRDS(paste0('../../04_cell_annotation/new_data_20250827/', stage_list[i], '_foregut.rds'))
  obj.list[[i]] <- UpdateSeuratObject(obj.list[[i]])
  obj.list[[i]]$Stage <- stage_list[i]
}

obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])


obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Trachea Mesenchymal Progenitor Cells',
                                                   'Lung Secondary Bronchi',
                                                   #'Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
#obj_sub <- UpdateSeuratObject(obj_sub)
#head(obj_sub)
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = HarmonyIntegration, orig.reduction = "pca", 
                           new.reduction = "harmony",
                           #k.weight = 50,k.anchor = 1,
                           verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "harmony", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))

saveRDS(object = obj_sub,file = 'E35_45_lung_progenitor_cells.rds')


obj_sub <- readRDS("E35_45_lung_progenitor_cells.rds")

#hdWGCNA
seurat_obj <- SetupForWGCNA(
  obj_sub,
  gene_select = "fraction",
  fraction = 0.05, 
  wgcna_name = "lung_progenitor_cells" 
)

seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype_gut"), 
  reduction = 'umap', 
  k = 25,
  max_shared = 10, 
  ident.group = 'celltype_gut' 
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  #group_name = "INH", 
  #group.by='cell_type', 
  assay = 'Spatial', 
  layer = 'data' 
)

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' 
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)

seurat_obj <- ConstructNetwork(seurat_obj,
                               #tom_name = 'INH'
)
PlotDendrogram(seurat_obj, main='Lung hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)

seurat_obj <- ModuleEigengenes(seurat_obj,group.by.vars="celltype_gut")

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(seurat_obj)
seurat_obj <- ResetModuleNames(seurat_obj,new_name = "Lung-M")

p <- PlotKMEs(seurat_obj, ncol=4)
p

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', 
  order=TRUE 
)
wrap_plots(plot_list, ncol=4)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

FeaturePlot(object = obj_sub,features = c("WNT5A","WNT11","WNT2B","HOXB5"))


###########GO enrichment
library(org.Gg.eg.db)
library(clusterProfiler)
library(GO.db)
lungM2 <- modules[modules$module == "Lung-M2",]
genes <- lungM2$gene_name
eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
id <- as.character(eg[, 2])

ego <- enrichGO(gene = id, OrgDb = "org.Gg.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
ego <- as.data.frame(ego)
ego <- ego[ego$pvalue < 0.05,]

top15_pvalue <- ego %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(15)
head(top15_pvalue)
top15_pvalue$Description <- str_to_title(top15_pvalue$Description)

pdf('Lung_Progenitor_cells_hdWGCNA_M2_GO.pdf',width = 10,height = 7)
ggplot(top15_pvalue, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)),  # 按p-value排序
           fill = RichFactor)) +
  scale_fill_distiller(palette = 'RdPu', direction = 1) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_classic() +
  labs(y = '', x = '-log10(p-value)')+
  theme(axis.text = element_text(size=12))
dev.off()

###KEGG
kegg_res <- enrichKEGG(gene = id,
                       organism = 'gga',
                       keyType = 'kegg',
                       pvalueCutoff = 1,
                       pAdjustMethod = 'BH')
kegg_res1 <- kegg_res@result
