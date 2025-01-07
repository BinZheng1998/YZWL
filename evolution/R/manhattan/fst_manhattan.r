setwd('E:/02-群体进化/02-鸡/phased_breed/BW_selection/plot/')
library(qqman)
library(CMplot)
library(readxl)

Fstfile <- read.table('../fst/chicken_BW_fst_10kbwindow_5kbstep.windowed.weir.fst',header = T,stringsAsFactors = F)
Fstfile2 <- Fstfile %>% filter(CHROM != 'Z')
Fstfile3 <- Fstfile2 %>% filter(CHROM != 'W')
Fstfile4 <- Fstfile3 %>% filter(CHROM != 'MT')

Fstfile4 <- Fstfile4[,c(1,2,5)]
Fstfile4$WEIGHTED_FST[Fstfile4$WEIGHTED_FST < 0] <- 0

SNP<-paste(Fstfile4[,1],Fstfile4[,2],sep = ":")
Fstfile4=cbind(SNP,Fstfile4)

colnames(Fstfile4)<-c("SNP","CHR","POS","Fst")

###删除ZW和MT染色体
library(dplyr)
library(tidyverse)
head(gwasResults)
head(Fstfile4)
Fstfile4$POS <- as.numeric(Fstfile4$POS)
Fstfile4$Fst <- as.numeric(Fstfile4$Fst)

# 1)  计算chr长度
chr_len <- Fstfile4 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(POS))
chr_len
chr_len <- chr_len %>% arrange(as.numeric(CHR))
chr_len

# 2） 计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)           #去掉数据集中的chr_len这一列
chr_pos


#3)   计算累计SNP的位置
Fst_pos <- chr_pos %>%
  left_join(Fstfile4, ., by="CHR") %>%
  arrange(CHR, POS) %>%
  mutate( BPcum = POS + total)
Fst_pos

ggplot(Fst_pos, aes(x = BPcum, y = Fst)) +
  geom_point(aes(color=as.factor(CHR)))

#计算普通fst的top5%阈值
options(scipen = 999)
sort_data <- as.numeric(sort(Fstfile4$Fst,decreasing = T))
percentile_5 <- quantile(sort_data, probs = 0.95)
percentile_5


FST_res <- Fstfile4[Fstfile4$Fst >= percentile_5,]
FST_res <- na.omit(FST_res)
write.table(FST_res,'FST_BW_95%_res.txt',sep = '\t',row.names = F)

head(Fst_pos)
library(dplyr)
#X_axis <-  Fst_pos %>% group_by(CHR) %>% summarize(center=max(BPcum)/2) 之前的代码报错了
X_axis_agg <- aggregate(BPcum ~ CHR, data = Fst_pos, FUN = function(x) (max(x) + min(x)) / 2)

X_axis <- X_axis_agg %>% arrange(as.numeric(CHR))
X_axis

nCHR<-length(unique(Fst_pos$CHR))
nCHR

Fst_pos <- Fst_pos %>% arrange(as.numeric(CHR))
unique(Fst_pos$CHR)
Fst_pos$CHR <- factor(Fst_pos$CHR,levels = unique(Fst_pos$CHR))
Fst_pos$zFST <- (as.numeric(Fst_pos$Fst)-mean(as.numeric(Fst_pos$Fst)))/sd(as.numeric(Fst_pos$Fst))

sort_data2 <- as.numeric(sort(Fst_pos$zFST,decreasing = T))
percentile_5 <- quantile(sort_data2, probs = 0.95)
percentile_5

p<-ggplot(Fst_pos, aes(x = BPcum, y = zFST)) +
  geom_point(aes(color=as.factor(CHR)), size = 0.5) +
  #scale_color_manual(values = rep(c("#4292c6","#08306b"),nCHR)) +
  scale_color_manual(values = rep(c("#000078","#8784C6"),nCHR)) +
  #设置横坐标刻度
  scale_x_continuous(label = X_axis$CHR, breaks=X_axis$BPcum, expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,30), expand = c(0,0)) + 
  geom_hline(yintercept = percentile_5,lty='dashed',color="red")+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        #axis.text.x = element_blank(),
        #axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        legend.position="none",
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  labs(x = "Chromosome")
p
ggsave(p,filename = 'chicken_BW_FST.png',dpi = 500,height = 4,width = 16)
