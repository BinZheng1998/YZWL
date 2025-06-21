setwd('~/project/01_evolution/convergent_evolution/07_result/20250616_manhattan/BW/pig/')
library(readxl)
library(dplyr)
library(tidyverse)

pos_dcms <- read.table('../../../20250605_DCMS_res/BW/pig/pig_BW_pos_DCMS_regions.txt',sep = '\t',header = T)
pos_dcms <- pos_dcms[pos_dcms$dcms > 0,]
neg_dcms <- read.table('../../../20250605_DCMS_res/BW/pig/pig_BW_neg_DCMS_regions.txt',sep = '\t',header = T)
neg_dcms <- neg_dcms[neg_dcms$dcms > 0,]
neg_dcms$dcms <- -neg_dcms$dcms

dcms <- rbind(pos_dcms,neg_dcms)
dcms2 <- dcms[,c(1,2,3,6)]
dcms2$POS <- (dcms2$BIN_START+dcms2$BIN_END)/2
dcms3 <- dcms2[,c(1,5,4)]
dcms3$CHROM <- gsub("chr","",dcms3$CHROM)
SNP<-paste(dcms3[,1],dcms3[,2],sep = ":")
dcms4=cbind(SNP,dcms3)
colnames(dcms4)<-c("SNP","CHR","POS","DCMS")

dcms4$POS <- as.numeric(dcms4$POS)
dcms4$DCMS <- as.numeric(dcms4$DCMS)

chr_len <- dcms4 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(POS))
chr_len
chr_len <- chr_len %>% arrange(as.numeric(CHR))
chr_len

chr_pos <- chr_len  %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  dplyr::select(-chr_len)           #去掉数据集中的chr_len这一列
chr_pos

Fst_pos <- chr_pos %>%
  left_join(dcms4, ., by="CHR") %>%
  arrange(CHR, POS) %>%
  mutate( BPcum = POS + total)
Fst_pos

options(scipen = 999)
X_axis_agg <- aggregate(BPcum ~ CHR, data = Fst_pos, FUN = function(x) (max(x) + min(x)) / 2)
X_axis <- X_axis_agg %>% arrange(as.numeric(CHR))
X_axis

nCHR<-length(unique(Fst_pos$CHR))
nCHR
Fst_pos <- Fst_pos %>% arrange(as.numeric(CHR))
unique(Fst_pos$CHR)
Fst_pos$CHR <- factor(Fst_pos$CHR,levels = unique(Fst_pos$CHR))

#X_axis$CHR[X_axis$CHR %in% c("11","13","15","16","17","19","20","22","23","24","25",
#                             "27","28","29","30","32","33","34","35","36","37")] <- ""

p<-ggplot(Fst_pos, aes(x = BPcum, y = DCMS)) +
  geom_point(aes(color=as.factor(CHR)), size = 0.5) +
  #scale_color_manual(values = rep(c("#4292c6","#08306b"),nCHR)) +
  scale_color_manual(values = rep(c("darkgrey","black"),nCHR)) +
  scale_x_continuous(label = X_axis$CHR, breaks=X_axis$BPcum,expand = expansion(mult=c(0.01,0.01))) + 
  scale_y_continuous(limits = c(-13,13)) + 
  geom_hline(yintercept = 0,color="grey")+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        axis.text = element_text(size=15),
        axis.title = element_text(size=15),
        #axis.ticks.x = element_blank(),
        axis.line = element_line(color = "grey"),
        axis.ticks = element_line(color = "grey"),
        legend.position="none",
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  xlab("Chromosome")+
  ylab("DCMS")
p
ggsave(p,filename = 'Pig_BW_DCMS.png',dpi = 500,height = 6,width = 16)
