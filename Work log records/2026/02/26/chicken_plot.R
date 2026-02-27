setwd('~/project/01_evolution/convergent_evolution/07_result/20250808_selection_region_plot/BW/multic_genes/')
library(ggplot2)
library(dplyr)
library(patchwork)

fst <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/fst/chicken_BW_10kbwindow_5kbstep.windowed.weir.fst',sep = '\t',header = T)

pi1 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWhigh_10kbwindow_5kbstep.windowed.pi',sep = '\t',header = T)
pi1$Group <- 'High BW'

pi2 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWlow_10kbwindow_5kbstep_noRJF.windowed.pi',sep = '\t',header = T)
pi2$Group <- 'Low BW'
pi <- rbind(pi1,pi2)

tajimaD1 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/tajimaD/chicken_BWlow_10kbwindow_noRJF.Tajima.D',sep = '\t',header = T)
tajimaD1$Group <- 'Low BW'

tajimaD2 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/tajimaD/chicken_BWhigh_10kbwindow.Tajima.D',sep = '\t',header = T)
tajimaD2$Group <- 'High BW'
tajimaD <- rbind(tajimaD1,tajimaD2)

df1 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
df1$Group <- "Low BW"
df1$Top5[df1$p_values < df1$threshold_005] <- "Neg"
head(df1)
df2 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
df2$Group <- "High BW"
df2$Top5[df2$p_values < df2$threshold_005] <- "Pos"
res1 <- rbind(df1,df2)


plot_selection_signals <- function(gene_name, chrom, start, end, flank = 250000) {
  
  # 计算绘图区间
  min_s <- start - flank
  max_s <- end + flank
  
  # --- 1. 数据筛选辅助函数 ---
  filter_data <- function(data) {
    data %>% filter(CHROM == chrom, BIN_START > min_s, BIN_START < max_s)
  }

  # --- 2. 准备各层数据 ---
  plot_fst <- filter_data(fst)
  plot_pi  <- filter_data(pi)      # 假设 pi 已经是 rbind 后的全局变量
  plot_td  <- filter_data(tajimaD) # 假设 tajimaD 已经是 rbind 后的全局变量
  plot_dcms <- filter_data(res1)   # 假设 res1 包含 DCMS 的 High/Low 结果

  # --- 3. 绘图模板 (保持样式一致) ---
  base_theme <- theme_classic() + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "right")

  # P1: Fst
  p1 <- ggplot(plot_fst, aes(BIN_START, WEIGHTED_FST)) +
    geom_line(color = "black") +
    labs(y = "Fst", title = paste("Region:", gene_name)) + base_theme

  # P2: Pi
  p2 <- ggplot(plot_pi, aes(BIN_START, PI, color = Group)) +
    geom_line() +
    labs(y = "Pi") + base_theme

  # P3: Tajima's D
  p3 <- ggplot(plot_td, aes(BIN_START, TajimaD, color = Group)) +
    geom_line() +
    labs(y = "Tajima's D") + base_theme

  # P4: DCMS
  p4 <- ggplot(plot_dcms, aes(BIN_START, dcms, color = Group)) +
    geom_line() +
    labs(y = "DCMS") + base_theme

  # P5: Gene Structure
  p_gene <- ggplot() +
    annotate("segment", x = start, xend = end, y = 1, yend = 1,
             arrow = arrow(type = "closed", length = unit(0.1, "inches")),
             color = "black", linewidth = 1) +
    annotate("text", x = (start + end)/2, y = 1.2, label = gene_name, 
             fontface = "italic", size = 4) +
    theme_void() +
    theme(axis.text.x = element_text(color = "black"),
          axis.line.x = element_line(color = "black")) +
    scale_x_continuous(labels = function(x) paste0(round(x/1e6, 2), " Mb"),
                       limits = c(min_s, max_s)) +
    ylim(0.75, 1.5)

  # --- 4. 拼图并返回 ---
  # 使用 patchwork 控制高度比例，让基因图层扁一点
  final_plot <- (p1 / p2 / p3 / p4 / p_gene) + 
    plot_layout(heights = c(1, 1, 1, 1, 0.5), guides = 'collect')
  
  return(final_plot)
}


library(rtracklayer)
library(dplyr)
get_gene_coords <- function(anno_path, target_ids) {
  message("Loading annotation file...")
  
  # rtracklayer 会根据后缀自动识别 GFF 或 GTF
  anno_data <- import(anno_path)
  anno_df <- as.data.frame(anno_data)
  
  # 兼容性处理列名
  # GFF3 常用 'ID' 和 'Name'，GTF 常用 'gene_id' 和 'gene_name'
  colnames(anno_df) <- gsub("gene_id", "ID", colnames(anno_df))
  colnames(anno_df) <- gsub("gene_name", "Name", colnames(anno_df))
  
  res_df <- anno_df %>%
    filter(type == "gene") %>%
    # 查找 ID 列或 Name 列（不区分大小写更稳妥）
    filter(ID %in% target_ids | Name %in% target_ids) %>%
    mutate(
      # 确保有名字显示
      final_name = coalesce(as.character(Name), as.character(ID)),
      # 染色体前缀统一化
      seqnames = as.character(seqnames),
      chr = if_else(grepl("^chr", seqnames, ignore.case = TRUE), 
                    seqnames, 
                    paste0("chr", seqnames))
    ) %>%
    select(name = final_name, chr, s = start, e = end) %>%
    distinct(name, .keep_all = TRUE)
  
  if (nrow(res_df) == 0) {
    warning("No genes found. Check if target_ids match the 'ID' or 'Name' field in your file.")
  }
  
  return(split(res_df, seq(nrow(res_df))) %>% lapply(as.list))
}

gff_gtf_file <- "~/ref/huxu_v23_ref/20250514_final_ref/20250529_chicken.gtf" # 替换为你的 GTF/GFF 路径

my_target_genes <- c("SOX5", "IGF1", "LCORL", "ADAMTS17", "BRCA1","CACNA1D","CACNA2D3","CCSER1","ELP4","FHIT","FMNL2",
                      "FNDC3B","LPP","LRRC7","MAGI2","MECOM","METTL15P1","NPAS3","OSBPL9","PAX6",
                    "R3HDM2","RBFOX1","RORA","RUNDC3B","SEMA3C","SEMA3D","SLC4A10","SLCO1A2","SORD","LOC107050760","SUPT3H",
                  "THSD7A","TOX","UNC13C","ZBTB20","ZMIZ1")

genes_to_plot <- get_gene_coords(gff_gtf_file, my_target_genes)

for (g in genes_to_plot) {
  p <- plot_selection_signals(g$name, g$chr, g$s, g$e)
  ggsave(paste0("chicken_Plot_", g$name, ".pdf"), p, width = 12, height = 10)
}
