setwd(dir = '~/project/01_evolution/convergent_evolution/07_result/20250908_PCA/plink/')
library(ggplot2)
library(dplyr)
library(randomcoloR)
pc <- read.table('chicken_2036samples.eigenvec')
#pc <- read.table('chicken_smallchr.eigenvec')
#pc <- read.table('chicken_midchr.eigenvec')
#pc <- read.table('chicken_bigchr.eigenvec')
pc1 <- pc[,c(2,3,4)]
colnames(pc1) <- c('sample','pc1','pc2')
#group <- read.table('../vcf2cluster/sample_group.txt',sep = '\t')
group <- read.table('sample_group2.txt',sep = '\t')
colnames(group) <- c('sample','group')
pc1$group <- group$group[match(pc1$sample,group$sample)]
pc2 <- pc1[pc1$pc1 < 0 & pc1$pc2 < 0.02,]
ggplot(data = pc2,aes(x = pc1,y = pc2,color = group))+
  #stat_ellipse(aes(fill = group),type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 1)+
  #labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","orange","pink"))+
  scale_color_manual(values = distinctColorPalette(length(unique(pc2$group))))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))

ggsave(p.pca1,filename = "PCA.pdf")

RJF <- pc1[pc1$group == "Gallus_gallus_gallus",]
RJF_group <- read.table('../metadata/RJF_group.txt',sep = '\t',header = T)
RJF$group1 <- RJF_group$group[match(RJF$sample,RJF_group$sample)]
RJF <- na.omit(RJF)
ggplot(data = RJF,aes(x = pc1,y = pc2,color = group1))+
  #stat_ellipse(aes(fill = group),type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 1)+
  #labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  #scale_fill_manual(values = c("purple","orange","pink"))+
  #scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
