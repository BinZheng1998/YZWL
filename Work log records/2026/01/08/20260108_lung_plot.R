setwd('~/project/03_3D_spatial/02_result/260108_lung_plot/')
library(Seurat)
library(ggunchull)
library(dplyr)
df1 <- readRDS('../../00_data/lung/E25_lung_bronchi.rds')
df1 <- NormalizeData(df1,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(df1,group.by = "Spatial_snn_res.0.5")
df1$celltype <- 'E25 Lung Epithelium'

df2 <- readRDS('../../00_data/lung/E25_lung_mesenchymal_1216.rds')
df2$celltype <- 'E25 Lung Mesenchyme'

df3 <- readRDS('../../00_data/lung/E35_lung_epithelium_251222.rds')
df3$celltype <- 'E35 Lung Epithelium'

df4 <- readRDS('../../00_data/lung/E35_lung_mesechymal_251220.rds')
df4$celltype <- 'E35 Lung Mesenchyme'

df5 <- readRDS('../../00_data/lung/E45_lung_epithelium_251220.rds')
df5$celltype <- 'E45 Lung Epithelium'

df6 <- readRDS('../../00_data/lung/E45_lung_mesenchymal_251220.rds')
df6$celltype <- 'E45 Lung Mesenchyme'

lung_combined <- merge(x = df1, 
                       y = c(df2,df3, df4, df5, df6), 
                       add.cell.ids = c("E25_Epi","E25_Mes", "E35_Epi", "E35_Mes", "E45_Epi", "E45_Mes"), 
                       project = "Chicken_Lung_Dev")
DefaultAssay(lung_combined) <- 'Spatial'
lung_combined <- NormalizeData(lung_combined,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)

DimPlot(lung_combined, group.by = 'celltype') 

library(scRNAtoolVis)
library(scales)
library(ggsci)
library(ggforce)
library(ggrepel)
obj <- lung_combined
head(obj)
obj@meta.data <- obj@meta.data %>% 
  mutate(celltype_lung = case_when(
    grepl("E25 Lung Mesenchyme", celltype) ~ "Lung Mesenchyme",
    grepl("E35 Lung Mesenchyme", celltype) ~ "Lung Mesenchyme",
    grepl("E45 Lung Mesenchyme", celltype) ~ "Lung Mesenchyme",
    TRUE ~ "Lung Epithelium"
  ))

obj@meta.data <- obj@meta.data %>% 
  mutate(Stage = case_when(
    grepl("E25 Lung Mesenchyme", celltype) ~ "E25",
    grepl("E25 Lung Epithelium", celltype) ~ "E25",
    grepl("E35 Lung Mesenchyme", celltype) ~ "E35",
    grepl("E35 Lung Epithelium", celltype) ~ "E35",
    TRUE ~ "E45"
  ))

meta_info <- cbind(obj@meta.data, Embeddings(obj, "umap"))
bad_cells1 <- rownames(meta_info)[
  meta_info$celltype_lung == "Lung Epithelium" &  # 锁定特定类型
  meta_info$umap_1 < 5 
]
bad_cells2 <- rownames(meta_info)[
  meta_info$celltype_lung == "Lung Mesenchyme" &  # 锁定特定类型
  meta_info$umap_1 > 4.5 
]
bad_cells <- c(bad_cells1, bad_cells2)
lung_clean <- subset(obj, cells = bad_cells, invert = TRUE)

meta <- cbind(lung_clean@meta.data, lung_clean@reductions$umap@cell.embeddings)
#plot1
meta_type_med1 <- meta %>% group_by(celltype) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p1<-ggplot(meta, aes(x = umap_1, y = umap_2)) +
  #stat_unchull(aes(fill=celltype,color=celltype),delta=0.15,alpha=0.01,lty=2,size=1,show.legend=F)+
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype), size = 0.2, show.legend = FALSE) +
  geom_text_repel(aes(label=celltype, color=celltype),fontface="bold", data=meta_type_med1, show.legend=F, size=4) +
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(),
        axis.title.x = element_blank(),
        axis.title = element_blank())
p1
#plot2
head(obj)
lung_epi <- subset(lung_clean, subset = celltype_lung == "Lung Epithelium")
meta2 <- cbind(lung_epi@meta.data, lung_epi@reductions$umap@cell.embeddings)
meta_type_med2 <- meta %>%
  filter(celltype_lung != "Lung Mesenchyme") %>%
  group_by(celltype_lung) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p2<-ggplot(meta2, aes(x = umap_1, y = umap_2)) +
  geom_mark_hull(data=meta,
    aes( color = celltype_lung), 
    concavity = 5,
    radius = unit(3, "mm"),
    expand = unit(3, "mm"),
    alpha = 0.1, 
    size = 1, 
    linetype = 2, 
    show.legend = F
  ) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype_lung), size = 0.2, show.legend = FALSE) +
  geom_text_repel(aes(label=celltype_lung, color=celltype_lung),fontface="bold", data=meta_type_med2, show.legend=F, size=4) +
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

p2

#plot3
head(obj)
lung_mes <- subset(lung_clean, subset = celltype_lung == "Lung Mesenchyme")
meta3 <- cbind(lung_mes@meta.data, lung_mes@reductions$umap@cell.embeddings)
meta_type_med3 <- meta %>%
  filter(celltype_lung != "Lung Epithelium") %>%
  group_by(celltype_lung) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p3<-ggplot(meta3, aes(x = umap_1, y = umap_2)) +
  geom_mark_hull(data=meta,
    aes( color = celltype_lung), 
    concavity = 5,
    radius = unit(3, "mm"),
    expand = unit(3, "mm"),
    alpha = 0.1, 
    size = 1, 
    linetype = 2, 
    show.legend = F
  ) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype_lung), size = 0.2, show.legend = FALSE) +
  geom_text_repel(aes(label=celltype_lung, color=celltype_lung),fontface="bold", data=meta_type_med3, show.legend=F, size=4) +
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())

p3

#plot4
head(obj)
lung_25 <- subset(lung_clean, subset = Stage == "E25")
meta4 <- cbind(lung_25@meta.data, lung_25@reductions$umap@cell.embeddings)
meta_type_med4 <- meta %>%
  filter(!Stage %in% c("E35","E45")) %>%
  group_by(celltype_lung) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p4<-ggplot(meta4, aes(x = umap_1, y = umap_2)) +
  geom_mark_hull(data=meta,
    aes( color = celltype_lung), 
    concavity = 5,
    radius = unit(3, "mm"),
    expand = unit(3, "mm"),
    alpha = 0.1, 
    size = 1, 
    linetype = 2, 
    show.legend = F
  ) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype_lung), size = 0.2, show.legend = FALSE) +
  #geom_text_repel(aes(label=celltype_lung, color=celltype_lung),fontface="bold", data=meta_type_med4, show.legend=F, size=5) +
  annotate(geom='text',label='E25',x=-5,y=7,size=6,fontface = "bold")+
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(),
        axis.title = element_blank())

p4

#plot5
head(obj)
lung_35 <- subset(lung_clean, subset = Stage == "E35")
meta5 <- cbind(lung_35@meta.data, lung_35@reductions$umap@cell.embeddings)
meta_type_med5 <- meta %>%
  filter(!Stage %in% c("E25","E45")) %>%
  group_by(celltype_lung) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p5<-ggplot(meta5, aes(x = umap_1, y = umap_2)) +
  geom_mark_hull(data=meta,
    aes( color = celltype_lung), 
    concavity = 5,
    radius = unit(3, "mm"),
    expand = unit(3, "mm"),
    alpha = 0.1, 
    size = 1, 
    linetype = 2, 
    show.legend = F
  ) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype_lung), size = 0.2, show.legend = FALSE) +
  #geom_text_repel(aes(label=celltype_lung, color=celltype_lung),fontface="bold", data=meta_type_med5, show.legend=F, size=5) +
  annotate(geom='text',label='E35',x=-5,y=7,size=6,fontface = "bold")+
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.title.y = element_blank(),
        axis.title = element_blank())

p5

#plot6
head(obj)
lung_45 <- subset(lung_clean, subset = Stage == "E45")
meta6 <- cbind(lung_45@meta.data, lung_45@reductions$umap@cell.embeddings)
meta_type_med6 <- meta %>%
  filter(!Stage %in% c("E25","E35")) %>%
  group_by(celltype_lung) %>%
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)
p6<-ggplot(meta6, aes(x = umap_1, y = umap_2)) +
  geom_mark_hull(data=meta,
    aes( color = celltype_lung), 
    concavity = 5,
    radius = unit(3, "mm"),
    expand = unit(3, "mm"),
    alpha = 0.1, 
    size = 1, 
    linetype = 2, 
    show.legend = F
  ) +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  geom_point(aes(color = celltype_lung), size = 0.2, show.legend = FALSE) +
  #geom_text_repel(aes(label=celltype_lung, color=celltype_lung),fontface="bold", data=meta_type_med6, show.legend=F, size=5) +
  annotate(geom='text',label='E45',x=-5,y=7,size=6,fontface = "bold")+
  scale_color_npg() +
  scale_fill_npg() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.title.y = element_blank(),
        axis.title = element_blank())

p6

combined_plot <- (p1 + p2 +p3)/(p4+p5+p6)
combined_plot
pdf('lung_umap.pdf',width = 12,height = 8)
combined_plot
dev.off()



text_colors <- pal_npg()(6)
p7 <- DotPlot(obj, features = c('NKX2-1','RIPPLY3','BMP4','ZFPM2'), group.by = 'celltype') +
  NoLegend() +
  scale_y_discrete(position = "right") + # Y轴放右边
  
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "italic", hjust = 0.5),
    axis.text.y.right = element_text(
      color = text_colors, 
      face = "bold"        
    )
  )

p7
pdf('lung_marker_gene.pdf',width = 6,height = 8)
p7
dev.off()
