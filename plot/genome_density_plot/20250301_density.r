setwd('E:/02-群体进化/07-结果/20250301-选择区间可视化/')
library(ggplot2)

#chicken
df <- read.table('E:/02-群体进化/02-鸡/phased_breed/converge/10kb_window_5kb_step/fst/chicken_fst_95%_regions.txt',header = T)
df <- df[,c(2,3,4)]
colnames(df)[1] <- "chr"
df$value <- 1
length <- read_excel('chicken_chr_length.xlsx')

library(dplyr)
library(tidyr)

fill_missing_intervals <- function(df1, df2) {
  df <- df1 %>% 
    left_join(df2, by = "chr") %>%
    arrange(chr, start)
  
  result <- list()
  
  for (chrom in unique(df$chr)) {
    chrom_df <- filter(df, chr == chrom)
    chrom_length <- chrom_df$length[1]
    
    all_positions <- data.frame(pos = 1:chrom_length, chr = chrom)
    
    all_positions$value <- 0
    
    for (i in 1:nrow(chrom_df)) {
      start_pos <- chrom_df$start[i]
      end_pos <- chrom_df$end[i]
      all_positions$value[start_pos:end_pos] <- 1
    }
    
    intervals <- all_positions %>%
      mutate(change = value != lag(value, default = -1)) %>%
      mutate(group = cumsum(change)) %>%
      group_by(chr, group, value) %>%
      summarise(start = min(pos), end = max(pos), .groups = "drop") %>%
      select(chr, start, end, value)
    
    result[[as.character(chrom)]] <- intervals
  }
  
  final_df <- bind_rows(result)
  return(final_df)
}

filled_df <- fill_missing_intervals(df, length)


intervals <- filled_df[,c(1,2,3,4)]
chr_lengths <- length[,c(1,2)]

intervals$chr <- factor(intervals$chr, levels = paste0("chr", 40:1))
chr_lengths$chr <- factor(chr_lengths$chr, levels = paste0("chr", 40:1))


chr_lengths <- chr_lengths %>%
  arrange(chr) %>%
  mutate(y_pos = row_number()) 


intervals <- intervals %>%
  left_join(chr_lengths, by = "chr")


p<-ggplot() +
  geom_rect(data = chr_lengths,
            aes(xmin = 0, xmax = length,
                ymin = y_pos - 0.4, ymax = y_pos + 0.4),
            fill = "white",color = "black") +
  geom_rect(data = intervals %>% filter(value == 1),
            aes(xmin = start, xmax = end,
                ymin = y_pos - 0.4, ymax = y_pos + 0.4),
            fill='white',color = "red",linewidth = 0.1) +
  geom_text(data = chr_lengths,
            aes(x = -0.001 * max(length), y = y_pos, label = chr),
            hjust = 1, size = 3) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(breaks = chr_lengths$y_pos, labels = chr_lengths$chr) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  xlab('')
p
ggsave("chicken_selection_fst.png",p,dpi=400,width = 16,height = 8)
