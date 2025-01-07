setwd('E:/02-群体进化/02-鸡/selection_BW/10kb_window_5kb_step/plot/')
library(qqman)
library(CMplot)
library(readxl)
library(dplyr)
library(tidyverse)
Fstfile <- read.table('../xpnsl/chicken_BW_xpnsl_5kb.txt',header = F,stringsAsFactors = F)
Fstfile <- Fstfile[Fstfile$V1 != "chrZ",]
Fstfile <- Fstfile[Fstfile$V1 != "chrW",]
Fstfile <- Fstfile[Fstfile$V1 != "chrMT",]

Fstfile$V9[is.na(Fstfile$V9)] <- -1
Fstfile$V10[is.na(Fstfile$V10)] <- -1
Fstfile <- Fstfile %>% mutate(XPNSL = rowMeans(select(., V9, V10), na.rm = TRUE))
Fstfile <- Fstfile[,c(1,2,11)]

SNP<-paste(Fstfile[,1],Fstfile[,2],sep = ":")
Fstfile=cbind(SNP,Fstfile)
colnames(Fstfile)<-c("SNP","CHR","POS","XPNSL")

Fstfile$POS <- as.numeric(Fstfile$POS)
Fstfile$XPNSL <- as.numeric(Fstfile$XPNSL)
#去掉chr字符
Fstfile$CHR <- as.numeric(gsub("chr", "", Fstfile$CHR))

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


ggplot(Fst_pos, aes(x = BPcum, y = XPNSL)) +
  geom_point(aes(color=as.factor(CHR)))

#计算普通fst的top5%阈值
options(scipen = 999)
sort_data <- as.numeric(sort(Fstfile$XPNSL,decreasing = T))
percentile_5 <- quantile(sort_data, probs = 0.95)
percentile_5

#只保留正向选择的结果
Fst_pos <- Fst_pos[Fst_pos$XPNSL>0,]

xpnsl_res <- Fstfile[Fstfile$XPNSL >= percentile_5,]
xpnsl_res <- na.omit(xpnsl_res)
write.table(xpnsl_res,'chicken_BW_XPNSL_95%_res.txt',sep = '\t',row.names = F)

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

p<-ggplot(Fst_pos, aes(x = BPcum, y = XPNSL)) +
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
ggsave(p,filename = 'chicken_BW_XPNSL.png',dpi = 500,height = 4,width = 16)
