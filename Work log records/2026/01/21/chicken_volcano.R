setwd('~/project/04_rnaseq/chicken/05_result/260121_volcano/')
#####volcnao
library(tidyverse)
library(ggrepel)
library(ggprism)
res_liver <- read.table('../260120_rmbatcheffect_PCA/Liver_BW_DESeq2_rmbatcheffect.txt',sep = '\t',header = T)
res_hypothalamus <- read.table('../260120_rmbatcheffect_PCA/Hypothalamus_BW_DESeq2_rmbatcheffect.txt',sep = '\t',header = T)
res_muscle <- read.table('../260120_rmbatcheffect_PCA/Muscle_BW_DESeq2_rmbatcheffect.txt',sep = '\t',header = T)
res_adipose <- read.table('../260120_rmbatcheffect_PCA/Adipose_BW_DESeq2_rmbatcheffect.txt',sep = '\t',header = T)
res_pituitary <- read.table('../260120_rmbatcheffect_PCA/Pituitary_BW_DESeq2_rmbatcheffect.txt',sep = '\t',header = T)
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
  ymin = -1,
  ymax = 1)

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


df1 <- df[abs(df$avg_log2FC) > 1,]
limit_val <- ceiling(max(abs(df$avg_log2FC)+0.1, na.rm = TRUE))
#my_jitter <- position_jitter(width = 0.4, seed = 123)

pdf('chicken_multiple_tissue_volcano.pdf',width = 12,height = 7)
ggplot(df1, aes(cluster, avg_log2FC, color = type)) +

  geom_rect(data = bg_vertical,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey95") +
  geom_jitter(stroke = 0, width = 0.4) + 
  geom_rect(data = bg_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                fill = cluster),
            inherit.aes = FALSE) +

  geom_text_repel(
    data = label_df,
    aes(cluster, avg_log2FC, label = gene),
    size = 3, color = "black", box.padding = 0.2,
    fontface = "italic",
    max.overlaps = 20) + 
  geom_text(aes(cluster, 0, label = chr),size = 4, color = "white", show.legend = FALSE) +
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
