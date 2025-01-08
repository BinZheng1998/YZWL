library(readxl)
library(dplyr)
library(tidyverse)

piH <- read.table("pi_high_10kbwindow_5kb_step.windowed.pi",header = T)
piL <- read.table('pi_low_10kbwindow_5kb_step.windowed.pi',header = T)

pi <- piL%>% 
  left_join(piH,by=c("CHROM","BIN_START","BIN_END")) %>% 
  dplyr::rename(piL=PI.x, piH=PI.y)%>%
  mutate(piHVSL=piH/piL,
         piLVSH=piL/piH)
###删除ZW和MT染色体
pi <- pi[pi$CHROM != "Z",]
pi <- pi[pi$CHROM != "W",]
pi <- pi[pi$CHROM != "MT",]

pi95 <-quantile(pi$piLVSH,0.95,na.rm=TRUE)

Fstfile <- pi[,c(1,2,9)]

SNP<-paste(Fstfile[,1],Fstfile[,2],sep = ":")
Fstfile=cbind(SNP,Fstfile)
colnames(Fstfile)<-c("SNP","CHR","POS","PI")


Fstfile$POS <- as.numeric(Fstfile$POS)
Fstfile$PI <- as.numeric(Fstfile$PI)

# 1)  计算chr长度
chr_len <- Fstfile %>% 
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
  left_join(Fstfile, ., by="CHR") %>%
  arrange(CHR, POS) %>%
  mutate( BPcum = POS + total)
Fst_pos

ggplot(Fst_pos, aes(x = BPcum, y = PI)) +
  geom_point(aes(color=as.factor(CHR)))

#计算top5%阈值
options(scipen = 999)
sort_data <- as.numeric(sort(Fstfile$PI,decreasing = T))
percentile_5 <- quantile(sort_data, probs = 0.95)
percentile_5


PI_res <- Fstfile[Fstfile$PI>= percentile_5,]
PI_res <- na.omit(PI_res)
write.table(PI_res,'PI_BW_95%_res.txt',sep = '\t',row.names = F)

#X_axis <-  Fst_pos %>% group_by(CHR) %>% summarize(center=max(BPcum)/2) 之前的代码报错了
X_axis_agg <- aggregate(BPcum ~ CHR, data = Fst_pos, FUN = function(x) (max(x) + min(x)) / 2)

X_axis <- X_axis_agg %>% arrange(as.numeric(CHR))
X_axis

nCHR<-length(unique(Fst_pos$CHR))
nCHR

Fst_pos <- Fst_pos %>% arrange(as.numeric(CHR))
unique(Fst_pos$CHR)
Fst_pos$CHR <- factor(Fst_pos$CHR,levels = unique(Fst_pos$CHR))
Fst_pos$log2pi <- log2(as.numeric(Fst_pos$PI))

sort_data2 <- as.numeric(sort(Fst_pos$log2pi,decreasing = T))
percentile_5 <- quantile(sort_data2, probs = 0.95)
percentile_5

p<-ggplot(Fst_pos, aes(x = BPcum, y = log2pi)) +
  geom_point(aes(color=as.factor(CHR)), size = 0.5) +
  #scale_color_manual(values = rep(c("#4292c6","#08306b"),nCHR)) +
  scale_color_manual(values = rep(c("#000078","#8784C6"),nCHR)) +
  #设置横坐标刻度
  scale_x_continuous(label = X_axis$CHR, breaks=X_axis$BPcum, expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,10), expand = c(0,0)) + 
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
ggsave(p,filename = 'BW_PI.png',dpi = 500,height = 4,width = 16)
