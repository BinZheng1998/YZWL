setwd('~/project/01_evolution/convergent_evolution/07_result/20251010_DCMS_plot_redraw/manhattan_plot/')
library(dplyr)
library(ggplot2)

# 读取正向和负向 DCMS 数据
pos_dcms <- read.table('../../20250801_BW_DCMS/cattle/cattle_BW_pos_DCMS_regions.txt', sep = '\t', header = TRUE)
pos_dcms$log10p <- -log10(pos_dcms$p_values)
n_pos <- nrow(pos_dcms)  # 正向点数量
n_top_pos <- ceiling(n_pos * 0.05)  # top 5% 点数，向上取整
pos_dcms$is_top5 <- "normal"  # 默认值
top_pos <- head(pos_dcms[order(-pos_dcms$dcms), ], n_top_pos)  # 按 dcms 降序取 top 5%
pos_dcms$is_top5[pos_dcms$dcms %in% top_pos$dcms] <- "top_pos"  # 标记为 top_pos
#pos_dcms <- pos_dcms[pos_dcms$dcms > 0, ]

# 读取并处理 neg_dcms 数据
neg_dcms <- read.table('../../20250801_BW_DCMS/cattle/cattle_BW_neg_DCMS_regions.txt', sep = '\t', header = TRUE)
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

#ggsave(p,filename = 'sheep_BW_DCMS_251110.pdf',dpi = 500,height = 6,width = 16)

library(GenomicRanges)
library(IRanges)
library(ggsci)
######chicken
chr_length <- read.table('../../20250723_chr_region_plot/cattle_chr_length.txt',sep = '\t',header = T)
chr_length$Chr <- paste0('chr',chr_length$Chr)

#pos
data <- read.table('../../20250801_BW_DCMS/cattle/cattle_BW_pos_DCMS_regions.txt',sep = '\t',header = T)
data1 <- data[data$p_values < data$threshold_005,]
data1$length <- data1$BIN_END-data1$BIN_START
gr_pos <- GRanges(
  seqnames = data1$CHROM,
  ranges = IRanges(start = data1$BIN_START, end = data1$BIN_END)
)
gr_reduced_pos <- reduce(gr_pos)
result_pos <- data.frame(
  chromosome = as.character(seqnames(gr_reduced_pos)),
  length = width(gr_reduced_pos)
)
result1_pos <- aggregate(length ~ chromosome, data = result_pos, sum)
result2_pos <- merge(chr_length,result1_pos,by.x='Chr',by.y='chromosome')
result2_pos$Selection <- 'pos'
result2_pos$prop <- result2_pos$length/result2_pos$End

#neg
df <- read.table('../../20250801_BW_DCMS/cattle/cattle_BW_neg_DCMS_regions.txt',sep = '\t',header = T)
df1 <- df[df$p_values < df$threshold_005,]
df1$length <- df1$BIN_END-df1$BIN_START
gr_neg <- GRanges(
  seqnames = df1$CHROM,
  ranges = IRanges(start = df1$BIN_START, end = df1$BIN_END)
)
gr_reduced_neg <- reduce(gr_neg)
result_neg <- data.frame(
  chromosome = as.character(seqnames(gr_reduced_neg)),
  length = width(gr_reduced_neg)
)
result1_neg <- aggregate(length ~ chromosome, data = result_neg, sum)
result2_neg <- merge(chr_length,result1_neg,by.x='Chr',by.y='chromosome')
result2_neg$Selection <- 'neg'
result2_neg$prop <- result2_neg$length/result2_neg$End

result3 <- rbind(result2_neg,result2_pos)
result3$prop <- result3$prop*100
result3$Chr <- factor(result3$Chr,chr_length$Chr)




# --- 1. 加载库 ---
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

# --- 2. 数据准备 (Data Prep) ---
df <- result3 %>%
  select(CHR = Chr, End, Method = Selection, value = prop) %>%
  mutate(Start = 0,
         CHR = as.character(CHR))

master_chr_order <- paste0("chr", 1:38)
chrs_in_your_data <- unique(df$CHR)
final_chr_levels <- master_chr_order[master_chr_order %in% chrs_in_your_data]

chr_info <- df %>%
  select(CHR, End) %>%
  distinct() %>%
  mutate(CHR = factor(CHR, levels = final_chr_levels)) %>%
  arrange(CHR) %>%
  mutate(
    start = lag(cumsum(as.numeric(End)), default = 0),
    position = start + (End / 2),
    end_pos = start + End
  )

df_with_pos <- df %>%
  left_join(chr_info, by = "CHR")

df_final <- df_with_pos %>%
  mutate(CHR = factor(CHR, levels = final_chr_levels)) %>%
  arrange(Method, CHR) %>%
  group_by(Method) %>%
  mutate(
    Value_next_chr = lead(value),
    Value_sliding_avg = (value + Value_next_chr) / 2
  ) %>%
  ungroup()

df_original <- df_final %>%
  select(
    CHR, 
    Method, 
    X_Position = start,
    Value_stacked = value
  ) %>% 
  mutate(DataType = "Original")

df_average <- df_final %>%
  select(
    CHR, 
    Method, 
    X_Position = position, 
    Value_stacked = Value_sliding_avg 
  ) %>%
  mutate(DataType = "Average")

df_endpoint <- df_final %>%
  filter(is.na(Value_sliding_avg)) %>% 
  select(
    CHR, Method,
    X_Position = end_pos,
    Value_stacked = value
  ) %>%
  mutate(DataType = "Average")

df_stacked <- bind_rows(df_original, df_average, df_endpoint) %>%
  filter(!is.na(Value_stacked)) %>%
  arrange(Method, X_Position, DataType)

# --- 3. 绘图参数 (Plot Parameters) ---
max(Fst_pos$Log10P)
min(Fst_pos$Log10P)
max(df_stacked$Value_stacked)
params_top <- list(
  dcms_breaks = c(0, 3, 6, 9, 12, 15),
  dcms_labels = c("0", "3", "6", "9", "12", "15"),
  dcms_limit = 15, 
  prop_breaks = c(0, 3,6,9, 12, 15),
  prop_labels = c("0", "3","6", "9", "12", "15"),
  prop_limit = 15, 
  scaling_factor = 15 / 15
)
params_top$sec_breaks_primary <- params_top$prop_breaks * params_top$scaling_factor

params_bottom <- list(
  dcms_breaks = c(0, 10, 20, 30, 40),
  dcms_labels = c("0", "10", "20", "30", "40"),
  dcms_limit = 40,
  prop_breaks = c(0, 3,6,9, 12, 15),
  prop_labels = c("0", "3","6", "9", "12", "15"),
  prop_limit = 15,
  scaling_factor = 40 / 15
)
params_bottom$sec_breaks_primary <- params_bottom$prop_breaks * params_bottom$scaling_factor

df_stacked_signed <- df_stacked %>%
  mutate(
    Y_scaled_pos = ifelse(Method == "pos", Value_stacked * params_top$scaling_factor, NA),
    Y_scaled_neg = ifelse(Method == "neg", Value_stacked * params_bottom$scaling_factor, NA)
  )

# --- 4. 创建 p_top ---
p_top <- ggplot() + 
  geom_point(
    data = Fst_pos %>% filter(Log10P >= 0), 
    aes(x = BPcum, y = Log10P, color = color_group), size = 1
  ) +
  scale_color_manual(values = color_values, drop = FALSE) +
  stat_smooth(
    data = df_stacked_signed %>% filter(Method == "pos"),
    aes(x = X_Position, y = Y_scaled_pos, ymin = 0, ymax = after_stat(y)),
    geom = "ribbon", fill = "grey50", alpha = 0.15, se = FALSE, 
    inherit.aes = FALSE, method = "loess", span = 0.1, 
    color = "white", linewidth = 0.8
  ) +
  scale_x_continuous(
    label = NULL, 
    breaks = chr_info$start,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    name = "",
    breaks = params_top$dcms_breaks,
    labels = params_top$dcms_labels,
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ . / params_top$scaling_factor,
      #name = "Proportion (%)",
      breaks = params_top$prop_breaks, 
      labels = params_top$prop_labels
    )
  ) +
  coord_cartesian(ylim = c(0, params_top$dcms_limit),expand = FALSE )+
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    legend.position = "none", axis.line = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5, unit = "pt")
  )
p_top

p_bottom <- ggplot() + 
  geom_point(
    data = Fst_pos %>% filter(Log10P < 0), 
    aes(x = BPcum, y = -Log10P, color = color_group),
    size = 1
  ) +
  scale_color_manual(values = color_values, drop = FALSE) +
  stat_smooth(
    data = df_stacked_signed %>% filter(Method == "neg"),
    aes(x = X_Position, y = Y_scaled_neg, ymin = 0, ymax = after_stat(y)),
    geom = "ribbon", fill = "grey50", alpha = 0.15, se = FALSE, 
    inherit.aes = FALSE, method = "loess", 
    span = 0.2, 
    color = "white", linewidth = 0.8
  ) +
  scale_x_continuous(
    label = NULL, 
    breaks = chr_info$start,
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_reverse(
    name = "",
    breaks = params_bottom$dcms_breaks,
    labels = params_bottom$dcms_labels,
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ . / params_bottom$scaling_factor,
     # name = "Proportion (%)",
      breaks = params_bottom$prop_breaks,
      labels = params_bottom$prop_labels
    )
  ) +
  coord_cartesian(ylim = c(params_bottom$dcms_limit, 0)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right = element_text(color = "black"),
    legend.position = "none", axis.line = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
  )
p_bottom

max_x_limit <- max(chr_info$end_pos)
p_labels <- ggplot(chr_info) +
  geom_text(aes(x = position, 
                y = 0, 
                label = stringr::str_remove(CHR, "chr")),
            size = 7, color = "black") +
  scale_x_continuous(
    label = NULL, 
    breaks = chr_info$start, 
    expand = expansion(mult = c(0, 0)),
    limits = c(0, max_x_limit) 
  ) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 5.5, b = 0, l = 5.5, unit = "pt"))
p_labels

cattle_plot <- p_top / p_labels / p_bottom +
  plot_layout(heights = c(10,0.5, 10))

cattle_plot
ggsave('cattle_BW_manhattan_251112.pdf',cattle_plot,width = 24,height = 10)
