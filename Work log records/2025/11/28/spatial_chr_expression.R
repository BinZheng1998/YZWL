setwd('~/project/03_3D_spatial/02_result/251127_chr_gene_expression/')
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(Matrix) # 处理稀疏矩阵
library(rtracklayer)

df1 <- readRDS('../../00_data/E2.5_sp_new_seurat_v5_cluster.rds')
df2 <- readRDS('../../00_data/E3.5_sp_new_seurat_v5_cluster.rds')
df3 <- readRDS('../../00_data/E4.5_sp_new_seurat_v5_cluster.rds')
head(df)
pseudo_bulk <- AggregateExpression(df, 
                                   group.by = "orig.ident", 
                                   assays = "Spatial", 
                                   slot = "counts", 
                                   return.seurat = FALSE)
count_mat <- pseudo_bulk$Spatial
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
chrom_result <- merged_data %>%
  group_by(chromosome) %>%
  summarise(across(where(is.numeric), mean)) %>%
  column_to_rownames("chromosome")
target_chroms <- c(paste0("chr", 1:38))
final_result <- chrom_result[rownames(chrom_result) %in% target_chroms, ]



library(ggplot2)
plot_data_box <- final_result %>%
  rownames_to_column(var = "chromosome") %>%
  pivot_longer(cols = -chromosome, names_to = "sample", values_to = "val") %>%
  mutate(log_tpm = log2(val + 1)) %>%
  mutate(type = case_when(
    chromosome %in% c(paste0("chr", 1:9), "chrZ", "chrW") ~ "Macro",
    chromosome %in% c(paste0("chr", 29:38),"chr16") ~ "Dot",
    TRUE ~ "Micro"
  )) %>%
  mutate(chromosome = gsub("chr", "", chromosome))

plot_data_box$chromosome <- reorder(plot_data_box$chromosome, plot_data_box$log_tpm, FUN = median)
plot_data_box$type <- factor(plot_data_box$type, levels = c("Macro", "Micro", "Dot"))


p_E25 <- ggplot(plot_data_box, aes(x = chromosome, y = log_tpm)) +
  geom_boxplot(outlier.shape = NA, aes(color = type), fill = "white", width = 0.7) +
   #geom_jitter(width = 0.2, aes(color = type), size = 1.5, alpha = 0.7) + # 如果样本少建议加上
  scale_color_manual(values = c(
    "Macro" = "#4E84C4", 
    "Micro" = "#F28E79", 
    "Dot"   = "#D14949"
  )) +
  theme_classic() +
  labs(color = "", y = "log2(CPM + 1)", x = "") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12)
  )

pdf('E25_chr_gene_expression.pdf',width = 12,height = 6)
print(p_E25)
dev.off()



#E35
p_E35 <- ggplot(plot_data_box, aes(x = chromosome, y = log_tpm)) +
  geom_boxplot(outlier.shape = NA, aes(color = type), fill = "white", width = 0.7) +
   #geom_jitter(width = 0.2, aes(color = type), size = 1.5, alpha = 0.7) + # 如果样本少建议加上
  scale_color_manual(values = c(
    "Macro" = "#4E84C4", 
    "Micro" = "#F28E79", 
    "Dot"   = "#D14949"
  )) +
  theme_classic() +
  labs(color = "", y = "log2(CPM + 1)", x = "") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12)
  )
pdf('E35_chr_gene_expression.pdf',width = 12,height = 6)
print(p_E35)
dev.off()


#E45
p_E45 <- ggplot(plot_data_box, aes(x = chromosome, y = log_tpm)) +
  geom_boxplot(outlier.shape = NA, aes(color = type), fill = "white", width = 0.7) +
   #geom_jitter(width = 0.2, aes(color = type), size = 1.5, alpha = 0.7) + # 如果样本少建议加上
  scale_color_manual(values = c(
    "Macro" = "#4E84C4", 
    "Micro" = "#F28E79", 
    "Dot"   = "#D14949"
  )) +
  theme_classic() +
  labs(color = "", y = "log2(CPM + 1)", x = "Chromosome") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12)
  )
pdf('E45_chr_gene_expression.pdf',width = 12,height = 6)
print(p_E45)
dev.off()

library(patchwork)

p_E25/p_E35/p_E45
