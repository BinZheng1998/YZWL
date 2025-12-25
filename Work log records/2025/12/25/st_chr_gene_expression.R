setwd('~/project/03_3D_spatial/02_result/251127_chr_gene_expression/')
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(Matrix) # 处理稀疏矩阵
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

df1 <- readRDS('../../00_data/E2.5_sp_new_seurat_v5_cluster.rds')
DefaultAssay(df1) <- "Spatial"
df1 <- DietSeurat(df1, assays = "Spatial", dimreducs = NULL, scale.data = FALSE)
df2 <- readRDS('../../00_data/E3.5_sp_new_seurat_v5_cluster.rds')
DefaultAssay(df2) <- "Spatial"
df2 <- DietSeurat(df2, assays = "Spatial", dimreducs = NULL, scale.data = FALSE)
df3 <- readRDS('../../00_data/E4.5_sp_new_seurat_v5_cluster.rds')
DefaultAssay(df3) <- "Spatial"
df3 <- DietSeurat(df3, assays = "Spatial", dimreducs = NULL, scale.data = FALSE)

pseudo_bulk <- AggregateExpression(df3, 
                                   group.by = "orig.ident", 
                                   assays = "Spatial", 
                                   slot = "counts", 
                                   return.seurat = FALSE)
count_mat <- pseudo_bulk$Spatial
head(count_mat)
cpm_mat <- t(t(count_mat) / colSums(count_mat) * 1e6)
cpm_df <- as.data.frame(cpm_mat) %>% 
  rownames_to_column("gene_id")

library(data.table)
dt <- fread('~/project/00_ref/20250529_chicken.gtf', header = FALSE, sep = "\t")
gene_map_dt <- dt[V3 == "gene", .(
  chromosome = V1,
  gene_id = sub('.*gene_id "([^"]+)".*', '\\1', V9),
  gene_name = sub('.*gene_name "([^"]+)".*', '\\1', V9)
)]
gene_map <- as.data.frame(gene_map_dt)
merged_data <- inner_join(cpm_df, gene_map, by = "gene_id")

# 1. 准备绘图数据 (修正版)
plot_data_box <- merged_data %>%
  # 【关键修改】: 不猜列名，而是排除掉非数值的元数据列
  # 意思就是：除了 gene_id, chromosome, gene_name 之外，剩下的列都当作样本表达量合并
  pivot_longer(cols = -c(gene_id, chromosome, gene_name), 
               names_to = "sample", 
               values_to = "cpm_val") %>%
  
  # A. 先取 Log (基因层面的 Log)
  mutate(log_tpm = log2(cpm_val + 1)) %>%
  
  # B. 定义染色体类型
  mutate(type = case_when(
    chromosome %in% c(paste0("chr", 1:9), "chrZ", "chrW") ~ "Macro",
    chromosome %in% c(paste0("chr", 29:38), "chr16") ~ "Dot",
    TRUE ~ "Micro"
  )) %>%
  mutate(chromosome = gsub("chr", "", chromosome)) %>%
  filter(chromosome %in% gsub("chr", "", target_chroms))

# 2. 排序 (按中位数排序)
plot_data_box$chromosome <- reorder(plot_data_box$chromosome, plot_data_box$log_tpm, FUN = median)
plot_data_box$type <- factor(plot_data_box$type, levels = c("Macro", "Micro", "Dot"))

# 3. 绘图 (Boxplot 展示的是该染色体上成百上千个基因的分布)
p_box3 <- ggplot(plot_data_box, aes(x = chromosome, y = log_tpm)) +
  geom_boxplot(outlier.shape = NA, aes(color = type), fill = "white", width = 0.7) +
  
  # 可选：如果你想看每个箱子的平均值在哪里，可以加上这一行
  # stat_summary(fun = mean, geom = "point", shape=18, size=2, color="blue") +
  
  scale_color_manual(values = c(
    "Macro" = "#4E84C4", 
    "Micro" = "#F28E79", 
    "Dot"   = "#D14949"
  )) +
  theme_classic() +
  labs(color = "", y = "log2 CPM+1", x = "Chromosome") +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color = "black", size = 10, angle = 0),
    axis.text.y = element_text(color = "black", size = 12)
  )
print(p_box3)

pdf('E25_chr_gene_expression.pdf',width = 12,height = 6)
print(p_box1)
dev.off()

pdf('E35_chr_gene_expression.pdf',width = 12,height = 6)
print(p_box2)
dev.off()

pdf('E45_chr_gene_expression.pdf',width = 12,height = 6)
print(p_box3)
dev.off()

library(patchwork)
pdf('E25_E35_E45_chr_gene_expression.pdf',width = 12,height = 12)
p_box1/p_box2/p_box3
dev.off()
