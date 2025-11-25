setwd('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/tables/')
library(dplyr)

#pos
df1 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_pos_DCMS_top_005_raw_regions_50kb_genes.txt',sep = '\t',header = T)
ref <- read.table('../chicken/DCMS/chicken_BW_10kwindow_5kstep_pos_DCMS_regions.txt',sep = '\t',header = T)
colnames(ref)[1:3] <- c("chromosome", "start", "end")
result <- inner_join(df1, ref, by = c("chromosome", "start", "end"))
res1 <- result[,c("chromosome", "start", "end","FST","PI","dcms","p_values","genes")]
res1$Species <- 'Chicken'
res1$Selection <- 'Positive Selection'

df2 <- read.table('../chicken/tajimaD/BW_high_chicken_10kwindow.Tajima.D',sep = '\t',header = T)
df2$BIN_START <- df2$BIN_START +1
df2$end <- df2$BIN_START+10000
colnames(df2)[c(1, 2, 4,5)] <- c("chromosome", "start",'High_BW_TajimaD', "end")
df3 <- read.table('../chicken/tajimaD/BW_low_chicken_10kwindow.Tajima.D',sep = '\t',header = T)
df3$BIN_START <- df3$BIN_START +1
df3$end <- df3$BIN_START+10000
colnames(df3)[c(1, 2, 4,5)] <- c("chromosome", "start",'Low_BW_TajimaD', "end")
df4 <- cbind(df2,df3[,c(4)])
colnames(df4)[6] <- "Low_BW_TajimaD"

res2 <- left_join(res1,df4,by = c("chromosome", "start", "end"))
res3 <- res2[c("chromosome", "start", "end","FST","PI","dcms","p_values","genes","Species","Selection","High_BW_TajimaD","Low_BW_TajimaD")]
colnames(res3)[5] <- "PI Ratio"
res4 <- res3[res3$`PI Ratio` >1 ,]

pi1 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWhigh_10kbwindow_5kbstep.windowed.pi',sep = '\t',header = T)
pi1$BIN_END <- pi1$BIN_END +1
colnames(pi1)[5] <- 'High_BW_PI'
pi2 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWlow_10kbwindow_5kbstep_noRJF.windowed.pi',sep = '\t',header = T)
pi2$BIN_END <- pi2$BIN_END +1
colnames(pi2)[5] <- 'Low_BW_PI'
pi <- inner_join(pi1,pi2,by=c("CHROM","BIN_START","BIN_END"))
colnames(pi)[1:3]<- c("chromosome", "start", "end")

res5 <- left_join(res4,pi,by = c("chromosome", "start", "end"))
res6 <- res5[c("chromosome", "start", "end","FST","PI Ratio","dcms","p_values","genes","Species","Selection","High_BW_TajimaD","Low_BW_TajimaD",'High_BW_PI','Low_BW_PI')]
write.table(res6,'chicken_BW_pos_top5_res.txt',row.names = F,sep = '\t')

#neg
df1 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/chicken/DCMS/chicken_BW_10kwindow_5kstep_neg_DCMS_top_005_raw_regions_50kb_genes.txt',sep = '\t',header = T)
ref <- read.table('../chicken/DCMS/chicken_BW_10kwindow_5kstep_neg_DCMS_regions.txt',sep = '\t',header = T)
colnames(ref)[1:3] <- c("chromosome", "start", "end")
result <- inner_join(df1, ref, by = c("chromosome", "start", "end"))
res1 <- result[,c("chromosome", "start", "end","FST","PI","dcms","p_values","genes")]
res1$Species <- 'Chicken'
res1$Selection <- 'Natural Selection'

df2 <- read.table('../chicken/tajimaD/BW_high_chicken_10kwindow.Tajima.D',sep = '\t',header = T)
df2$BIN_START <- df2$BIN_START +1
df2$end <- df2$BIN_START+10000
colnames(df2)[c(1, 2, 4,5)] <- c("chromosome", "start",'High_BW_TajimaD', "end")
df3 <- read.table('../chicken/tajimaD/BW_low_chicken_10kwindow.Tajima.D',sep = '\t',header = T)
df3$BIN_START <- df3$BIN_START +1
df3$end <- df3$BIN_START+10000
colnames(df3)[c(1, 2, 4,5)] <- c("chromosome", "start",'Low_BW_TajimaD', "end")
df4 <- cbind(df2,df3[,c(4)])
colnames(df4)[6] <- "Low_BW_TajimaD"

res2 <- left_join(res1,df4,by = c("chromosome", "start", "end"))
res3 <- res2[c("chromosome", "start", "end","FST","PI","dcms","p_values","genes","Species","Selection","High_BW_TajimaD","Low_BW_TajimaD")]
colnames(res3)[5] <- "PI Ratio"
res4 <- res3[res3$`PI Ratio` <1 ,]

pi1 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWhigh_10kbwindow_5kbstep.windowed.pi',sep = '\t',header = T)
pi1$BIN_END <- pi1$BIN_END +1
colnames(pi1)[5] <- 'High_BW_PI'
pi2 <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/1.10kbwindow/pi/chicken_BWlow_10kbwindow_5kbstep_noRJF.windowed.pi',sep = '\t',header = T)
pi2$BIN_END <- pi2$BIN_END +1
colnames(pi2)[5] <- 'Low_BW_PI'
pi <- inner_join(pi1,pi2,by=c("CHROM","BIN_START","BIN_END"))
colnames(pi)[1:3]<- c("chromosome", "start", "end")

res5 <- left_join(res4,pi,by = c("chromosome", "start", "end"))
res6 <- res5[c("chromosome", "start", "end","FST","PI Ratio","dcms","p_values","genes","Species","Selection","High_BW_TajimaD","Low_BW_TajimaD",'High_BW_PI','Low_BW_PI')]
write.table(res6,'chicken_BW_neg_top5_res.txt',row.names = F,sep = '\t')

data1 <- read.table('chicken_BW_neg_top5_res.txt',sep = '\t',header = T)
data2 <- read.table('chicken_BW_pos_top5_res.txt',sep = '\t',header = T)
data <- rbind(data1,data2)
write.table(data,'chicken_BW_all_top5_res.txt',row.names = F,sep = '\t')
