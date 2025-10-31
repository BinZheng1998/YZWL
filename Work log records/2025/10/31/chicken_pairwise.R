setwd('~/project/01_evolution/convergent_evolution/07_result/20251029_BWhigh_vs_low_pairwise_species/01result/chicken/test/')
# 加载必要的包
library(dplyr)
library(readr)

#pos
files <- list.files(pattern = "_pos_DCMS_regions.txt$", full.names = TRUE)
all_data <- list()
for (i in seq_along(files)) {
  file_path <- files[i]
  data <- read.table(file_path, header=T,sep='\t')
  data$file_id <- basename(file_path)
  data <- data %>% dplyr::select(CHROM, BIN_START, BIN_END, p_values, threshold_005, file_id)
  
  all_data[[i]] <- data
  print(paste("已读取文件:", basename(file_path)))
}

full_data <- bind_rows(all_data)
full_data$region_key <- paste(full_data$CHROM, full_data$BIN_START, full_data$BIN_END, sep = "_")
full_data$is_significant <- full_data$p_values < full_data$threshold_005

pos <- full_data %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    total_files = n_distinct(file_id),
    sig_files = sum(is_significant),
    .groups = 'drop'
  ) %>%
  arrange(CHROM, BIN_START) 
pos$prop <- pos$sig_files/pos$total_files
print(head(pos, 10))
write_tsv(pos, "chicken_pairwise_pos_summary.tsv", na = "")

#neg

files <- list.files(pattern = "_neg_DCMS_regions.txt$", full.names = TRUE)
all_data <- list()
for (i in seq_along(files)) {
  file_path <- files[i]
  data <- read.table(file_path, header=T,sep='\t')
  data$file_id <- basename(file_path)
  data <- data %>% dplyr::select(CHROM, BIN_START, BIN_END, p_values, threshold_005, file_id)
  
  all_data[[i]] <- data
  print(paste("已读取文件:", basename(file_path)))
}

full_data <- bind_rows(all_data)
full_data$region_key <- paste(full_data$CHROM, full_data$BIN_START, full_data$BIN_END, sep = "_")
full_data$is_significant <- full_data$p_values < full_data$threshold_005

neg <- full_data %>%
  group_by(CHROM, BIN_START, BIN_END) %>%
  summarise(
    total_files = n_distinct(file_id),
    sig_files = sum(is_significant),
    .groups = 'drop'
  ) %>%
  arrange(CHROM, BIN_START) 
neg$prop <- neg$sig_files/neg$total_files
print(head(neg, 10))
write_tsv(neg, "chicken_pairwise_neg_summary.tsv", na = "")


#画图
n_pos <- nrow(pos)  
n_top_pos <- ceiling(n_pos * 0.05) 
pos$is_top5 <- "normal" 
top_pos <- head(pos[order(-pos$prop), ], n_top_pos) 
pos$is_top5[pos$prop %in% top_pos$prop] <- "top_pos" 
top5pos <- pos[pos$is_top5 == "top_pos",]
write_tsv(top5pos, "chicken_pairwise_pos_top5_regions.tsv", na = "")

n_neg <- nrow(neg)
n_top_neg <- ceiling(n_neg * 0.05)
neg$is_top5 <- "normal"
top_neg <- head(neg[order(-neg$prop), ], n_top_neg)
neg$is_top5[neg$prop %in% top_neg$prop] <- "top_neg"
neg$prop <- -neg$prop
top5neg <- neg[neg$is_top5 == "top_neg",]
write_tsv(top5neg, "chicken_pairwise_neg_top5_regions.tsv", na = "")

summary_table <- rbind(pos,neg)
summary_table$CHROM <- gsub("chr", "", summary_table$CHROM)
summary_table$POS <- (summary_table$BIN_START+summary_table$BIN_END)/2
summary_table1 <- summary_table[, c(1, 8, 6,7)]
SNP <- paste(summary_table1$CHROM, summary_table1$POS, sep = ":")
head(SNP)
summary_table2 <- cbind(SNP, summary_table1)
#rm(list = c("SNP","summary_table1"))
colnames(summary_table2) <- c("SNP", "CHR", "POS", "Prop","top5")

# 确保 POS 和 DCMS 为数值型
summary_table2$POS <- as.numeric(summary_table2$POS)
summary_table2$DCMS <- as.numeric(summary_table2$Prop)

chr_len <- summary_table2 %>% 
  group_by(CHR) %>% 
  summarise(chr_len = max(POS))
chr_len <- chr_len %>% arrange(as.numeric(CHR))

# 计算染色体累积位置
chr_pos <- chr_len %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len)

# 合并染色体位置信息
Fst_pos <- chr_pos %>%
  left_join(summary_table2, ., by = "CHR") %>%
  arrange(CHR, POS) %>%
  mutate(BPcum = POS + total)

# 设置科学计数法选项
options(scipen = 999)

# 计算横轴标签位置
X_axis_agg <- aggregate(BPcum ~ CHR, data = Fst_pos, FUN = function(x) (max(x) + min(x)) / 2)
X_axis <- X_axis_agg %>% arrange(as.numeric(CHR))

# 设置染色体数量和顺序
nCHR <- length(unique(Fst_pos$CHR))
Fst_pos <- Fst_pos %>% arrange(as.numeric(CHR))
Fst_pos$CHR <- factor(Fst_pos$CHR, levels = unique(Fst_pos$CHR))


# 创建颜色分组
Fst_pos$color_group <- paste(Fst_pos$CHR, Fst_pos$top5, sep = "_")

# 创建颜色映射
color_values <- rep(c("grey50", "grey80"), ceiling(nCHR / 2))[1:nCHR]  # 普通点交替灰色
names(color_values) <- paste(unique(Fst_pos$CHR), "normal", sep = "_")  # 命名普通点

# 为 top 5% 数据点分配颜色，按染色体交替深浅
color_values_top <- rep(c("dark", "light"), ceiling(nCHR / 2))[1:nCHR]  # 交替标记
names(color_values_top) <- unique(Fst_pos$CHR)
top_pos_colors <- ifelse(color_values_top == "dark", "#3c5488ff", "#8491b4ff")
names(top_pos_colors) <- paste(names(color_values_top), "top_pos", sep = "_")
top_neg_colors <- ifelse(color_values_top == "dark", "#dc0000ff", "#f39b7fff")
names(top_neg_colors) <- paste(names(color_values_top), "top_neg", sep = "_")
color_values <- c(color_values, top_pos_colors, top_neg_colors)

X_axis$CHR[X_axis$CHR %in% c("11","13","15","16","17","19","20","22","23","24","25",
                             "27","28","29","30","32","33","34","35","36","37")] <- ""

p <- ggplot(Fst_pos, aes(x = BPcum, y = Prop)) +
  geom_point(aes(color = color_group), size = 1) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(label = X_axis$CHR, breaks = X_axis$BPcum, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-1, 1), 
                     breaks = c(-1,-0.75,-0.5,-0.25,0, 0.25, 0.5, 0.75, 1),
                     labels = c("1","0.75","0.5","0.25","0", "0.25", "0.5", "0.75", "1")) +
  geom_hline(yintercept = 0, color = "white", linewidth = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 15),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("Chromosome") +
  ylab("Proportion(%)")

p
ggsave(p,filename = 'chicken_BW_pairwise_DCMS.pdf',dpi = 500,height = 6,width = 16)
