setwd('~/project/01_evolution/convergent_evolution/07_result/20251110_sheep_BW_new_analysis/')
library(dplyr)
library(ggplot2)

# 读取正向和负向 DCMS 数据
pos_dcms <- read.table('DCMS/sheep_BW_DCMS_251110_pos_DCMS_regions.txt', sep = '\t', header = TRUE)
pos_dcms$log10p <- -log10(pos_dcms$p_values)
n_pos <- nrow(pos_dcms)  # 正向点数量
n_top_pos <- ceiling(n_pos * 0.05)  # top 5% 点数，向上取整
pos_dcms$is_top5 <- "normal"  # 默认值
top_pos <- head(pos_dcms[order(-pos_dcms$dcms), ], n_top_pos)  # 按 dcms 降序取 top 5%
pos_dcms$is_top5[pos_dcms$dcms %in% top_pos$dcms] <- "top_pos"  # 标记为 top_pos
#pos_dcms <- pos_dcms[pos_dcms$dcms > 0, ]

# 读取并处理 neg_dcms 数据
neg_dcms <- read.table('DCMS/sheep_BW_DCMS_251110_neg_DCMS_regions.txt', sep = '\t', header = TRUE)
neg_dcms$log10p <- -log10(neg_dcms$p_values)
n_neg <- nrow(neg_dcms)  # 负向点数量
n_top_neg <- ceiling(n_neg * 0.05)  # top 5% 点数，向上取整
neg_dcms$is_top5 <- "normal"  # 默认值
top_neg <- head(neg_dcms[order(-neg_dcms$dcms), ], n_top_neg)  # 按 dcms 升序取 top 5%（绝对值最大）
neg_dcms$is_top5[neg_dcms$dcms %in% top_neg$dcms] <- "top_neg"  # 标记为 top_neg
neg_dcms$log10p <- -neg_dcms$log10p

# 合并正向和负向数据
dcms <- rbind(pos_dcms, neg_dcms)
dcms2 <- dcms[, c(1, 2, 3, 12,13)]
dcms2$POS <- (dcms2$BIN_START + dcms2$BIN_END) / 2
dcms3 <- dcms2[, c(1, 6, 4,5)]
dcms3$CHROM <- gsub("chr", "", dcms3$CHROM)
SNP <- paste(dcms3[, 1], dcms3[, 2], sep = ":")
dcms4 <- cbind(SNP, dcms3)
colnames(dcms4) <- c("SNP", "CHR", "POS", "Log10P","top5")

# 确保 POS 和 DCMS 为数值型
dcms4$POS <- as.numeric(dcms4$POS)
dcms4$DCMS <- as.numeric(dcms4$Log10P)

# 将 is_top5 列合并到 dcms4
dcms4$is_top5 <- dcms$is_top5

# 计算染色体长度
chr_len <- dcms4 %>% 
  group_by(CHR) %>% 
  summarise(chr_len = max(POS))
chr_len <- chr_len %>% arrange(as.numeric(CHR))

# 计算染色体累积位置
chr_pos <- chr_len %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len)

# 合并染色体位置信息
Fst_pos <- chr_pos %>%
  left_join(dcms4, ., by = "CHR") %>%
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
Fst_pos$color_group <- paste(Fst_pos$CHR, Fst_pos$is_top5, sep = "_")

# 创建颜色映射
color_values <- rep(c("grey50", "grey80"), ceiling(nCHR / 2))[1:nCHR]  # 普通点交替灰色
names(color_values) <- paste(unique(Fst_pos$CHR), "normal", sep = "_")  # 命名普通点

# 为 top 5% 数据点分配颜色，按染色体交替深浅
color_values_top <- rep(c("dark", "light"), ceiling(nCHR / 2))[1:nCHR]  # 交替标记
names(color_values_top) <- unique(Fst_pos$CHR)

# 创建 top_pos 和 top_neg 的颜色映射
top_pos_colors <- ifelse(color_values_top == "dark", "#3c5488ff", "#8491b4ff")
names(top_pos_colors) <- paste(names(color_values_top), "top_pos", sep = "_")
top_neg_colors <- ifelse(color_values_top == "dark", "#dc0000ff", "#f39b7fff")
names(top_neg_colors) <- paste(names(color_values_top), "top_neg", sep = "_")

# 合并所有颜色映射
color_values <- c(color_values, top_pos_colors, top_neg_colors)

# 绘制散点图
p <- ggplot(Fst_pos, aes(x = BPcum, y = Log10P)) +
  geom_point(aes(color = color_group), size = 0.5) +
  scale_color_manual(values = color_values) +
  scale_x_continuous(label = X_axis$CHR, breaks = X_axis$BPcum, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-21, 21), 
                     breaks = c(-20,-15, -10, -5, 0,5,10,15,20),
                     labels = c("20","15","10", "5","0", "5","10", "15","20")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title = element_text(size = 20),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlab("") +
  ylab("DCMS")

p
ggsave(p,filename = 'sheep_BW_DCMS_251110.pdf',dpi = 500,height = 6,width = 16)
