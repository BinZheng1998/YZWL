setwd('~/project/04_rnaseq/chicken/05_result/260113_volcano/')
library(DESeq2)
library(ggplot2)
library(sva)
library(limma)
library(patchwork)
library(reshape2)
data <- read.csv('../../03_stringtie_2.2.1/re_quant/gene_count_matrix.csv')
data$gene_id <- gsub(".*\\|", "", data$gene_id)
data <- as.data.frame(data)
data1 <- data[!duplicated(data$gene_id), ]
rownames(data1) <- data1$gene_id
data2 <- data1[,-1]
data3 <- data2[, colSums(data2, na.rm = TRUE) >= 1000000]
detected_genes_gt1 <- colSums(data3 > 1)
summary(detected_genes_gt1)
hist(detected_genes_gt1, 
     main = "Histogram of Detected Genes per Sample",
     xlab = "Number of Detected Genes (TPM > 1)",
     breaks = 30) 
#去除极端值
Q1 <- quantile(detected_genes_gt1, 0.25)
Q3 <- quantile(detected_genes_gt1, 0.75)
IQR_val <- IQR(detected_genes_gt1)
whisker_lower <- Q1 - 1.5 * IQR_val
whisker_upper <- Q3 + 1.5 * IQR_val
samples_to_keep <- names(detected_genes_gt1[ detected_genes_gt1 <= whisker_upper])
data4 <- data3[, samples_to_keep]

df <- read.table("../../metadata/chicken_rnaseq_metadata.txt",sep = '\t',header = T,fill=T)
df <- df[df$Filter == "PASS",]
df2 <- df[, c(2, 4,8,11, 14,16)]
head(df2)
colnames(df2) <- c('Sample','Project','Breed', 'BW', 'Tissue','Instrument')
df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose","Pituitary"),]
df2 <- df2[df2$BW %in% c("High","Low"),]
df2$Group <- paste0(df2$BW,sep = '_',df2$Tissue)
df3 <- df2[df2$Tissue == "Pituitary",]

valid_cols <- intersect(df3$Sample, colnames(data4))
data5 <- data4[, valid_cols]
keep <- rowSums(data5 > 5) >= 10 
data5 <- data5[keep, ]
data6 <- data5 +1 

head(data6)
coldata <- df3[df3$Sample %in% colnames(data6),]
rownames(coldata) <- coldata[,1]
coldata$Group <- as.factor(coldata$Group)
dds<-DESeqDataSetFromMatrix(countData=data6,colData=coldata,design=~Group)
dds1 <- DESeq(dds)
res <- results(dds1, contrast = c("Group", "High_Pituitary", "Low_Pituitary"))
res2 <- as.data.frame(res)
res2$gene <- rownames(res2)
write.table(res2,'Pituitary_BW_DESeq2_raw.txt',sep = '\t',row.names = F)

#res2 <- res2[!grepl("STRG.", res2$gene), ]


dds <- DESeqDataSetFromMatrix(countData = data6, 
                              colData = coldata, 
                              design = ~ Group)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
dat_for_sva <- assay(vsd)
mod <- model.matrix(~ Group, data = coldata)
mod0 <- model.matrix(~ 1, data = coldata)
n.sv <- num.sv(dat_for_sva, mod, method="be")
n.sv

svobj <- sva(dat_for_sva, mod, mod0, n.sv = n.sv)
#svobj$sv

sv_data <- as.data.frame(svobj$sv)
colnames(sv_data) <- paste0("SV", 1:ncol(sv_data))
sv_data$Group <- coldata$Group 
r_sq_values <- sapply(1:n.sv, function(i) {
  fit <- lm(sv_data[,i] ~ as.numeric(as.factor(sv_data$Group)))
  summary(fit)$r.squared
})
bad_svs <- which(r_sq_values > 0.5) # 阈值 0.5，如果相关性超过这个值，说明该SV可能是生物学信号
cat("警告：以下 SV 与你的分组高度相关，建议剔除：\n")
if(length(bad_svs) > 0) {
  print(names(r_sq_values)[bad_svs])
  print(paste("R-squared:", round(r_sq_values[bad_svs], 2)))
} else {
  cat("好消息！没有发现与分组严重冲突的 SV。\n")
}
barplot(r_sq_values, names.arg = paste0("SV", 1:n.sv), 
        las=2, main="SV vs Group Correlation", ylab="R-squared")
abline(h=0.5, col="red", lty=2)



sv_data <- as.data.frame(svobj$sv)
colnames(sv_data) <- paste0("SV", 1:ncol(sv_data))
coldata_with_sv <- cbind(coldata, sv_data)
sv_names <- paste0("SV", 1:svobj$n.sv)
design_formula <- as.formula(paste0("~ ", paste(c(sv_names, "Group"), collapse = " + ")))
dds_sva <- DESeqDataSetFromMatrix(countData = data6, # 这里用原始 counts
                                  colData = coldata_with_sv, 
                                  design = design_formula)
dds_sva <- DESeq(dds_sva)
res <- results(dds_sva, contrast = c("Group", "High_Pituitary", "Low_Pituitary"))
res1 <- as.data.frame(res)
res1$gene <- rownames(res1)
res_pituitary <- res1
write.table(res1,'Pituitary_BW_DESeq2_rmbatcheffect.txt',sep = '\t',row.names = F)



sv_matrix <- svobj$sv
dat_cleaned <- removeBatchEffect(dat_for_sva, 
                                 covariates = sv_matrix, 
                                 design = mod)


plot_pca <- function(mat, info, title) {
  pca <- prcomp(t(mat))
  df <- as.data.frame(pca$x)
  df$Group <- info$Group
  df$Project <- info$Project
  perc <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  ggplot(df, aes(PC1, PC2, shape = Group, color = Project)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = title, 
         x = paste0("PC1 (", perc[1], "%)"), 
         y = paste0("PC2 (", perc[2], "%)")) +
    theme_bw()
}
pca <- prcomp(t(dat_for_sva))
pcadf <- as.data.frame(pca$x)

p1 <- plot_pca(dat_for_sva, coldata, "Raw")+theme(legend.position = "none")
p2 <- plot_pca(dat_cleaned, coldata, "After SVA (22 SVs removed)")
pdf('20260113_SVA_Pituitary.pdf',width = 12,height = 5)
p1 + p2
dev.off()






#####volcnao
library(tidyverse)
library(ggrepel)
library(ggprism)
data_list <- list(
  "Hypothalamus"  = res_hypothalamus,
  "Liver"  = res_liver,
  "Muscle"   = res_muscle,
  "Adipose" = res_adipose,
  "Pituitary" = res_pituitary
)

df <- map_dfr(data_list, function(x) {
  x <- as.data.frame(x)
  if (!"gene" %in% names(x)) {
    x <- rownames_to_column(x, "gene")
  }
  
  x %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    rename(
      avg_log2FC = log2FoldChange,
      p_val_adj = padj
    ) %>%
    select(gene, avg_log2FC, p_val_adj)
    
}, .id = "cluster") %>%
  mutate(
    chr  = as.character(cluster),
    type = if_else(avg_log2FC > 0, "High BW", "Low BW")
  )

df <- df[!grepl("STRG.", df$gene), ]
df <- df[!grepl("MTgene", df$gene), ]

df$cluster <- factor(df$cluster, levels = c("Hypothalamus", "Liver", "Muscle", "Adipose","Pituitary"))
cluster_colors <- c("#f39b7fff","#8491b4ff","#91d1c2ff","#7e6148ff","#4dbbd5ff") 

# 水平背景 (X轴标签条)
bg_df <- tibble(
  cluster = levels(df$cluster),
  xmin = seq_along(levels(df$cluster)) - 0.48,
  xmax = seq_along(levels(df$cluster)) + 0.48,
  ymin = -0.5,
  ymax = 0.5)

# 垂直背景 (灰色背景柱)
bg_vertical <- df %>%
  group_by(cluster) %>%
  summarise(
    ymin = min(avg_log2FC, na.rm = TRUE) - 0.1,
    ymax = max(avg_log2FC, na.rm = TRUE) + 0.1,
    .groups = "drop") %>%
  mutate(
    xmin = seq_along(cluster) - 0.48,
    xmax = seq_along(cluster) + 0.48)

# 标签筛选 (筛选 Top 5 上调和下调)
label_df <- df %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1) %>%
  group_by(chr) %>%
  group_modify(~ bind_rows(
    slice_max(.x, avg_log2FC, n = 5),
    slice_min(.x, avg_log2FC, n = 5))) %>%
  ungroup()


df1 <- df[abs(df$avg_log2FC) > 0.6,]
limit_val <- ceiling(max(abs(df$avg_log2FC)+0.1, na.rm = TRUE))

my_jitter <- position_jitter(width = 0.4, seed = 123)
pdf('chicken_multiple_tissue_volcano.pdf',width = 12,height = 7)
ggplot(df1, aes(cluster, avg_log2FC, color = type)) +
  geom_rect(data = bg_vertical,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey95") +
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = cluster),
            inherit.aes = FALSE) +
  geom_jitter(stroke = 0, width = 0.4) + 
  geom_text_repel(
    data = label_df,
    aes(cluster, avg_log2FC, label = gene),
    size = 3, color = "black", box.padding = 0.2,
    fontface = "italic",
    max.overlaps = 20) + 
  geom_text(aes(cluster, 0, label = chr),position = my_jitter,
            size = 5, color = "white", show.legend = FALSE) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_color_manual(values = c("#3c5488ff", "#dc0000ff")) +
  scale_y_continuous(
    limits = c(-limit_val, limit_val),
    guide = "prism_offset_minor") +
  labs(x = NULL, y = "log2 Fold Change") +
  guides(color = guide_legend(
    override.aes = list(size = 5, shape = 19))) +
  theme_prism(base_line_size = 0.3) +
  theme(
    panel.background = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.y  = element_text(size = 10),
    axis.title   = element_text(size = 11),
    legend.position = c(0.08, 0.9),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(margin = margin(l = 0)))
dev.off()
library(ggsci)
library(scales) # 需要加载这个包
show_col(pal_npg("nrc")(10))
