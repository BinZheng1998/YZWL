setwd('~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/plot/')
library(ggplot2)
library(Seurat)
library(patchwork)
library(aplot)
data1 <- readRDS('../E2.5_foregut.rds')
#data1@meta.data$celltype_gut <- recode(data1@meta.data$celltype_gut,"Foregut Mesothelial cells" = "Foregut Mesothelial Cells")
head(data1@meta.data)
table(data1$celltype_gut)
#data1@meta.data[31] <- NULL
all_markers <- FindAllMarkers(object = data1,group.by = 'celltype_gut',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

marker_genes <-c('FILIP1','RAB3IP','CREB5','ALDH1A2','AQP1','AQP3',
                 'RIPPLY3','SHH','MME','WNT2','WNT2B','SMOC2','BARX1','NKX3-2',
                 'UTS2B','SLC51B','FGG','FGA','ABCC9','VWDE',
                 'SOX2','CLDN1',
                 'PDX1','RFX6','GUCA2B','SLC7A9','INS','GCG','SOX17','HHEX')
vln.dat=FetchData(data1,c(marker_genes,"celltype_gut"))
vln.dat$Cell <- rownames(vln.dat)
res <- reshape2::melt(vln.dat, id.vars = c("Cell","celltype_gut"), 
                      measure.vars = marker_genes,
                      variable.name = "gene", 
                      value.name = "Expr") %>%
  group_by(celltype_gut,gene) %>% #分组
  mutate(fillcolor=mean(Expr)) #计算均值

table(data1$celltype_gut)
celltype <- c('Foregut Cells1','Foregut Cells2','Foregut Cells3','Lung Bud1','Lung Bud2',
              'Proventriculus Progenitor Cells','Gizzard Progenitor Cells','Gizzard-Duodenum Progenitor Cells',
              'Hepatocyte Progenitor Cells','Hepatocyte Mesenchymal Progenitor Cells',
              'Foregut Mucosa Epithelial Progenitor Cells',
              'Duodenum Progenitor Cells','Duodenum-Midgut Progenitor Cells',
              'Pancreatic Bud','Hepatopancreatic Ductal Progenitor Cells'
)
res$celltype_gut <- factor(res$celltype_gut,levels = celltype)

table(data1$celltype_gut)

my_colors <- c(
  "Duodenum Enterocyte-Neurons Cells" = "#1F77B4", # 蓝色
  "Duodenum Mesenchymal Progenitor Cells" = "#FF7F0E",
  "Duodenum Mucosa Epithelial Cells" = "#2CA02C", # 橙色
  "Duodenum Progenitor Cells" = "#D62728", # 绿色
  "Duodenum-Midgut Progenitor Cells" = "#9467BD", # 红色
  "Foregut Cells1" = "#8C564B",  # 紫色
  "Foregut Cells2" = "#E377C2", # 最深蓝
  "Foregut Cells3" = "#7F7F5F", 
  "Foregut Mesothelial Cells" = "#BCBD22",
  "Foregut Mucosa Epithelial Progenitor Cells" = "#17BECF",
  "Gizzard Mesenchymal Progenitor Cells" = "#AEC7E8",  # 最浅蓝
  "Gizzard Mucosa Epithelial Cells" = "#FFBB78", # 最深绿
  "Gizzard Progenitor Cells" = "#98DF8A", 
  "Gizzard-Duodenum Mesenchymal Progenitor Cells" = "#FF9896",
  "Gizzard-Duodenum Progenitor Cells" = "#C5B0D5",
  "Pancreatic Bud" = "#C49C94",  # 最浅绿
  "Spleno-Pancreatic Mesenchymal Progenitor Cells" = "#F7B6D2", # 最深紫
  "Proventriculus Mesothelial Cells" = "#C7C7C7", 
  "Proventriculus Mucosa Epithelial Cells" = "#DBDB8D",
  "Proventriculus Progenitor Cells" = "#9EDAE5",
  "Hepatocyte Mesenchymal Progenitor Cells" = "#393B79",  # 最浅紫
  "Hepatocyte Progenitor Cells" = "#637939", # 最深橙（近乎棕）
  "Hepatopancreatic Ductal Progenitor Cells" = "#8C6D31", 
  "Lung Bud1" = "#843C39",
  "Lung Bud2" = "#7B4173",
  "Lung Mesenchymal Progenitor Cells" = "#5254A3", 
  "Lung Primary Bronchi" = "#6BAED6",
  "Lung Progenitor Cells" = "#9E9AC8",
  "Lung Secondary Bronchi" = "#969696",
  "Trachea Mesenchymal Progenitor Cells" = "#6B6ECF",
  "Colacal Epithelial Cells" = "#FD8D3C",
  "Gut Mesothelial Cells" = "#E6550D",
  "Hindgut Mucosa Epithelial Cells" = "#31A399",
  "Hindgut Progenitor Cells" = "#756BB1",
  "Midgut-Hindgut Progenitor Cells" = "#636383"
  
)

metadata <- read.table('celltype_tissue.txt',sep = '\t',header = T,fill = T)
metadata1 <- metadata[,c(1,3)]
res1 <- merge(metadata1,res,by = 'celltype_gut',all=T)

#res$cluster_id <- as.numeric(factor(res$celltype_gut)) - 1
point_data <- data.frame(table(res$celltype_gut))
point_data$Var1 <- factor(point_data$Var1,levels = celltype)
colnames(point_data)[1] <- 'celltype_gut'
point_data <- merge(point_data,metadata1,by='celltype_gut',all=T)
#point_data$cluster_id <- as.numeric(factor(point_data$Var1)) - 1
point_data$x <- 0

unique(res$gene)
p_main <- ggplot(res, aes(x = Expr, y = celltype_gut, fill = celltype_gut)) +
  geom_violin(scale = 'width', color = 'white', size = 0.45, alpha = 0.8) +
  facet_wrap(~ gene, scales = 'free_x', nrow = 1) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 2)) +
  theme_bw() +
  labs(x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", size = 1),
    legend.position = 'none',
    strip.background = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    #axis.text.x = element_text(size = 10, color = 'black'),
    strip.text.x = element_text(angle = 90, size = 10, hjust = 0, vjust = 0.5,color = "black",face = 'italic'),
    strip.placement = "outside",
    panel.spacing.x = unit(0.1, "lines"),
    panel.border = element_rect(size = 1, color = "black")
  )
p_main
# 点图（独立 y 轴图例）
p_dots <- ggplot(point_data, aes(x = x, y = celltype_gut)) +
  geom_point(aes(fill = celltype_gut), size = 7, shape = 21, color = "white") +
  geom_text(aes(label = cluster_id), size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(expand = c(0, 0), breaks = NULL) +  # 隐藏 x 轴
  scale_y_discrete(position = "left") +
  theme_void() +
  theme(
    axis.text = element_blank(),
    legend.position = "none"
  )

p_combined <- p_main %>% insert_left(p_dots,width = 0.02)
#p_combined <-  p_dots %>% insert_left(p_main,width = 0.99)
p_combined
ggsave('E2.5_cellmarker_violin.pdf',p_combined,width = 20,height = 7,dpi = 300)

###########################2
cell_number <- data.frame(table(data1$celltype_gut))
colnames(cell_number) <- c('celltype_gut','number')
library(dplyr)
cell_number <- merge(cell_number,point_data,by='celltype_gut')
cell_number$cluster_id <- as.character(cell_number$cluster_id)

cell_number <- cell_number %>%
  mutate(Tissue = case_when(
    grepl("\\b(15|16)\\b", cluster_id) ~ "Pancreas",
    grepl("\\b(20|21|22)\\b", cluster_id) ~ "Liver",
    grepl("\\b(30|31|32|33|34)\\b", cluster_id) ~ "Mid-Hindgut",
    grepl("\\b(0|1|2|3|4)\\b", cluster_id) ~ "Duodenum", 
    grepl("\\b(10|11|12|13|14)\\b", cluster_id) ~ "Gizzard",
    grepl("\\b(17|18|19)\\b", cluster_id) ~ "Proventriculus",
    grepl("\\b(23|24|25|26|27|28|29)\\b", cluster_id) ~ "Lung",
    grepl("\\b(5|6|7|8|9)\\b", cluster_id) ~ "Foregut",
    TRUE ~ "Other"
  ))
cell_number$prop <- cell_number$number/sum(cell_number$number)
color <- data.frame(my_colors)
color$celltype_gut <- rownames(color)

cell_number <- merge(cell_number,color,by='celltype_gut')

library(WeightedTreemaps)
head(cell_number)
library(dplyr)
cell_number1 <- cell_number %>%
  group_by(Tissue) %>%
  mutate(Prop_Sum = sum(prop),  # 计算同一组织的prop总和
         Tissue1 = paste0(Tissue, "\n", sprintf("%.1f%%", Prop_Sum * 100))) %>%  # 创建格式化标签
  ungroup()

tm <- voronoiTreemap(
  cell_number1,
  levels = c("Tissue1", "cluster_id"),
  fun = sum,
  cell_size = "prop",
  shape = "rounded_rect",#rounded_rect
  seed = 123
)


pdf('E2.5_treemap.pdf',width = 8,height = 6)
drawTreemap(
  tm,
  color_type = "categorical",
  color_level = 2,
  color_palette = cell_number1$my_colors,
  legend = F,
  label_autoscale=F,
  label_size = 3,
  label_level = 1:2,
  label_color = c('black','white'),
  border_size = c(6,1),
  border_level = 1:2
)
dev.off()
