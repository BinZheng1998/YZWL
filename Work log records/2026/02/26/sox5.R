setwd('~/project/01_evolution/convergent_evolution/07_result/20250808_selection_region_plot/BW/SOX5/')
library(ggplot2)
library(dplyr)
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr1",]
min_site<-67069623-250000
max_site <- 67717269+250000
df1 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df1$Selection <- "Neg"
df1$Top5[df1$p_values < df1$threshold_005] <- "Neg"
head(df1)
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr1",]
min_site<-67069623-250000
max_site <- 67717269+250000
df2 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df2$Selection <- "Pos"
df2$Top5[df2$p_values < df2$threshold_005] <- "Pos"
res1 <- rbind(df1,df2)

res1_plot <- res1 %>%
  mutate(
    Top5_label = ifelse(Top5 %in% c("Neg", "Pos"), Selection, NA)
  )

p1 <- ggplot(res1_plot, aes(x = BIN_START, y = dcms)) +
  geom_point(data = filter(res1_plot, is.na(Top5_label)), 
             aes(shape = Selection), color = "grey80", size = 1, stroke = 0.5) +
  geom_point(data = filter(res1_plot, !is.na(Top5_label)), 
             aes(color = Top5_label, shape = Selection), size = 1.5) +
  theme_classic() +
  scale_shape_manual(
    name = "Selection Direction",
    values = c("Neg" = 17, "Pos" = 16), 
    labels = c("Neg" = "Negative Selection", "Pos" = "Positive Selection")
  ) +
  scale_color_manual(
    name = "Top 5% Significance",
    values = c("Neg" = "#E41A1C", "Pos" = "#377EB8"),
    labels = c("Neg" = "Top 5% (Neg)", "Pos" = "Top 5% (Pos)"),
    na.translate = FALSE # 核心：不显示 NA 的图例
  ) +
  theme(axis.title.y = element_text(face = "italic")) +
  labs(y = 'DCMS', x = 'Chromosome 1 position') +
  annotate("segment", x = 67069623, y = 6.6, xend = 67717269, yend = 6.6,
           arrow = arrow(type = "closed", length = unit(0.05, "inches")),
           color = "grey40", linewidth = 0.35) +
  annotate('text', x = (67069623 + 67717269) / 2, y = 6.9, 
           label = 'SOX5', size = 3, fontface = 'italic')

p1
ggsave('chicken_SOX5_dcms.pdf',p1,dpi = 300,width = 6,height = 4)


#cattle
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/cattle/DCMS/cattle_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr5",]
min_site<-86798896-250000
max_site <- 87265093+250000
df1 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df1$Selection <- "Neg"
df1$Top5[df1$p_values < df1$threshold_005] <- "Neg"
head(df1)
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/cattle/DCMS/cattle_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr5",]
min_site<-86798896-250000
max_site <- 87265093+250000
df2 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df2$Selection <- "Pos"
df2$Top5[df2$p_values < df2$threshold_005] <- "Pos"
res1 <- rbind(df1,df2)

res1_plot <- res1 %>%
  mutate(
    Top5_label = ifelse(Top5 %in% c("Neg", "Pos"), Selection, NA)
  )

p2 <- ggplot(res1_plot, aes(x = BIN_START, y = dcms)) +
  geom_point(data = filter(res1_plot, is.na(Top5_label)), 
             aes(shape = Selection), color = "grey80", size = 1, stroke = 0.5) +
  geom_point(data = filter(res1_plot, !is.na(Top5_label)), 
             aes(color = Top5_label, shape = Selection), size = 1.5) +
  theme_classic() +
  scale_shape_manual(
    name = "Selection Direction",
    values = c("Neg" = 17, "Pos" = 16), 
    labels = c("Neg" = "Negative Selection", "Pos" = "Positive Selection")
  ) +
  scale_color_manual(
    name = "Top 5% Significance",
    values = c("Neg" = "#E41A1C", "Pos" = "#377EB8"),
    labels = c("Neg" = "Top 5% (Neg)", "Pos" = "Top 5% (Pos)"),
    na.translate = FALSE # 核心：不显示 NA 的图例
  ) +
  theme(axis.title.y = element_text(face = "italic")) +
  labs(y = 'DCMS', x = 'Chromosome 5 position') +
  annotate("segment", x = 86798896, y = 4.6, xend = 87265093, yend = 4.6,
           arrow = arrow(type = "closed", length = unit(0.05, "inches")),
           color = "grey40", linewidth = 0.35) +
  annotate('text', x = (86798896 + 87265093) / 2, y = 4.9, 
           label = 'SOX5', size = 3, fontface = 'italic')

p2
ggsave('cattle_SOX5_dcms.pdf',p2,dpi = 300,width = 6,height = 4)


#pig
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/pig/DCMS/pig_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr5",]
min_site<-49159950-250000
max_site <- 50166463+250000
df1 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df1$Selection <- "Neg"
df1$Top5[df1$p_values < df1$threshold_005] <- "Neg"
head(df1)
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/pig/DCMS/pig_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr5",]
min_site<-49159950-250000
max_site <- 50166463+250000
df2 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df2$Selection <- "Pos"
df2$Top5[df2$p_values < df2$threshold_005] <- "Pos"
res1 <- rbind(df1,df2)

res1_plot <- res1 %>%
  mutate(
    Top5_label = ifelse(Top5 %in% c("Neg", "Pos"), Selection, NA)
  )

p2 <- ggplot(res1_plot, aes(x = BIN_START, y = dcms)) +
  geom_point(data = filter(res1_plot, is.na(Top5_label)), 
             aes(shape = Selection), color = "grey80", size = 1, stroke = 0.5) +
  geom_point(data = filter(res1_plot, !is.na(Top5_label)), 
             aes(color = Top5_label, shape = Selection), size = 1.5) +
  theme_classic() +
  scale_shape_manual(
    name = "Selection Direction",
    values = c("Neg" = 17, "Pos" = 16), 
    labels = c("Neg" = "Negative Selection", "Pos" = "Positive Selection")
  ) +
  scale_color_manual(
    name = "Top 5% Significance",
    values = c("Neg" = "#E41A1C", "Pos" = "#377EB8"),
    labels = c("Neg" = "Top 5% (Neg)", "Pos" = "Top 5% (Pos)"),
    na.translate = FALSE # 核心：不显示 NA 的图例
  ) +
  theme(axis.title.y = element_text(face = "italic")) +
  labs(y = 'DCMS', x = 'Chromosome 5 position') +
  annotate("segment", x = 49159950, y = 4.6, xend = 50166463, yend = 4.6,
           arrow = arrow(type = "closed", length = unit(0.05, "inches")),
           color = "grey40", linewidth = 0.35) +
  annotate('text', x = (49159950 + 50166463) / 2, y = 4.9, 
           label = 'SOX5', size = 3, fontface = 'italic')

p2
ggsave('pig_SOX5_dcms.pdf',p2,dpi = 300,width = 6,height = 4)

#sheep
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/sheep/DCMS/sheep_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr3",]
min_site<-205132332-250000
max_site <- 205599843+250000
df1 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df1$Selection <- "Neg"
df1$Top5[df1$p_values < df1$threshold_005] <- "Neg"
head(df1)
df <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/sheep/DCMS/sheep_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
df <- df[df$CHROM == "chr3",]
min_site<-205132332-250000
max_site <- 205599843+250000
df2 <- df[df$BIN_START > min_site & df$BIN_END < max_site,]
df2$Selection <- "Pos"
df2$Top5[df2$p_values < df2$threshold_005] <- "Pos"
res1 <- rbind(df1,df2)

res1_plot <- res1 %>%
  mutate(
    Top5_label = ifelse(Top5 %in% c("Neg", "Pos"), Selection, NA)
  )

p4 <- ggplot(res1_plot, aes(x = BIN_START, y = dcms)) +
  geom_point(data = filter(res1_plot, is.na(Top5_label)), 
             aes(shape = Selection), color = "grey80", size = 1, stroke = 0.5) +
  geom_point(data = filter(res1_plot, !is.na(Top5_label)), 
             aes(color = Top5_label, shape = Selection), size = 1.5) +
  theme_classic() +
  scale_shape_manual(
    name = "Selection Direction",
    values = c("Neg" = 17, "Pos" = 16), 
    labels = c("Neg" = "Negative Selection", "Pos" = "Positive Selection")
  ) +
  scale_color_manual(
    name = "Top 5% Significance",
    values = c("Neg" = "#E41A1C", "Pos" = "#377EB8"),
    labels = c("Neg" = "Top 5% (Neg)", "Pos" = "Top 5% (Pos)"),
    na.translate = FALSE # 核心：不显示 NA 的图例
  ) +
  theme(axis.title.y = element_text(face = "italic")) +
  labs(y = 'DCMS', x = 'Chromosome 3 position') +
  annotate("segment", x = 205132332, y = 9.6, xend = 205599843, yend = 9.6,
           arrow = arrow(type = "closed", length = unit(0.05, "inches")),
           color = "grey40", linewidth = 0.35) +
  annotate('text', x = (205132332 + 205599843) / 2, y = 9.9, 
           label = 'SOX5', size = 3, fontface = 'italic')

p4
ggsave('sheep_SOX5_dcms.pdf',p4,dpi = 300,width = 6,height = 4)

pdf('FTO_4species_neg_dcms.pdf',width = 14,height = 3)
p1|p2|p3|p4
dev.off()
