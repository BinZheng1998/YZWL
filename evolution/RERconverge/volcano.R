setwd("~/project/10_RER/result/RER/")
library(RERconverge)
library(ape)
data <- read.table('RERconverge_res_20250708.txt',sep = '\t',header = T)
head(data)
data$sig <- ifelse(data$permpval < 0.05,ifelse(data$Rho > 0, "Up", "Down"),"no")
data$logP <- -log10(data$permpval)
cut_off <- -log10(0.05)
ggplot(data,aes(x=Rho,y=logP))+
  geom_point(aes(fill=sig),color='white',size=2,shape=21)+
  scale_fill_manual(values = c("Down" = "blue", "no" = "NA", "Up" = "red"))+
  scale_x_continuous(limits = c(-0.6,0.6))+
  geom_hline(yintercept = cut_off,lty=2,color='black')+
  theme_bw()+
  theme(legend.position = "none")+
  ylab('-log10(Pvalue)')

data$label <- ifelse(data$permpval < 0.01 & abs(data$Rho) >= 0.4, as.character(data$human), "")
library(ggrepel)
library(ggsci)
library("scales")
mypal =pal_npg("nrc", alpha =0.7)(9)
mypal
show_col(mypal)
p<-ggplot(data, aes(x = Rho, y = logP)) +
  geom_point(data = subset(data, sig == "no"),
    aes(x = Rho, y = logP),shape = 21,color = "grey",size = 3) +
  geom_point(data = subset(data, sig != "no"),
    aes(x = Rho, y = logP, color = sig),size = 3,) +
  scale_color_manual(values = c("Down" = "#3C5488B2", "Up" = "#E64B35B2")) +
  scale_x_continuous(limits = c(-0.6,0.6))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(x = "RER Rho",y = "-log10(P-value)") +
  theme_classic() +
  theme(legend.position = "none",panel.grid = element_blank())+
  geom_text_repel(data = data,aes(x = Rho, y = logP, label = label),size = 3, show.legend = F)
p
ggsave('RER_volcano.pdf',p,dpi=500,width = 4,height = 4)
