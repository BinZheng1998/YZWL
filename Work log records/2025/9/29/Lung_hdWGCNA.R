setwd('~/project/03_3D_spatial/04_cell_annotation/new_data_2509/lung')
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
enableWGCNAThreads(nThreads = 8)
set.seed(12345)
options(future.globals.maxSize = 100* 1e10)
obj.list = list()
stage_list <- c('E2.5','E3.5','E4.5')
#library(Seurat)

for (i in 1:length(stage_list)) {
  obj.list[[i]] <- readRDS(paste0('../../new_data_20250827/', stage_list[i], '_foregut.rds'))
  obj.list[[i]] <- UpdateSeuratObject(obj.list[[i]])
  obj.list[[i]]$Stage <- stage_list[i]
}

obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])


obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                   'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Bud1','Lung Primary Bronchi'))
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                           new.reduction = "integrated.rpca",
                           k.weight = 50,k.anchor = 1,
                           verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))
saveRDS(obj_sub,'LungBud1.rds')

p <- DimPlot(obj_sub, group.by='celltype_gut', label=TRUE) +
  umap_theme() +NoLegend()+
  ggtitle('Lung')

p

seurat_obj <- SetupForWGCNA(
  obj_sub,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "lung" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype_gut", "sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype_gut' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  #group_name = "INH", # 选择感兴趣的细胞类型分析
  #group.by='cell_type', # 选择细胞类型所在列
  assay = 'Spatial', # using RNA assay
  layer = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table)


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  #tom_name = 'INH'
)
PlotDendrogram(seurat_obj, main='Lung hdWGCNA Dendrogram')

TOM <- GetTOM(seurat_obj)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="celltype_gut"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj
)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Lung-M"
)


p <- PlotKMEs(seurat_obj, ncol=5)
p

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', 
  order=TRUE 
)
wrap_plots(plot_list, ncol=6)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
write.table(modules,'lung_Primary_Bronchi_hdWGCNA.txt',row.names = F)


################################
#提取cellID进行空间构象绘图
head(seurat_obj)
seurat_obj@misc$lung
ME_matrix <- seurat_obj@misc$lung$MEs
module_name <- "Lung-M19"
threshold <- quantile(ME_matrix[, module_name], 0.85) # 或 quantile(ME_matrix[, module_name], 0.75) 等
cell_ids <- rownames(ME_matrix)[ME_matrix[, module_name] > threshold]
seurat_obj$ModuleCells <- ifelse(colnames(seurat_obj) %in% cell_ids, "ModuleCells", "Other")

head(seurat_obj@meta.data)
DimPlot(seurat_obj, group.by = "ModuleCells", label = FALSE, pt.size = 0.5) +
  ggtitle("Cells Associated with Module") +
  umap_theme() +NoLegend()

E35_lung_M19 <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$ModuleCells == "ModuleCells" & seurat_obj@meta.data$Stage == "E3.5"]
write.table(E35_lung_M19,'E35_lung_M19.txt',row.names = F)
E45_lung_M19 <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$ModuleCells == "ModuleCells" & seurat_obj@meta.data$Stage == "E4.5"]
write.table(E45_lung_M19,'E45_lung_M19.txt',row.names = F)
E35_lung <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$Stage == "E3.5"]
write.table(E35_lung,'E35_lung.txt',row.names = F)
E45_lung <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$Stage == "E4.5"]
write.table(E45_lung,'E45_lung.txt',row.names = F)




###########GO enrichment
library(org.Gg.eg.db)
library(clusterProfiler)
library(GO.db)
lungM19 <- modules[modules$module == "Lung-M19",]
genes <- lungM19$gene_name
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

pdf('LungBud1_hdWGCNA_M19_GO.pdf',width = 10,height = 7)
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
