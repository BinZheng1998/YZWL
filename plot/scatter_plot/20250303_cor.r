library(ggplot2)
library(readxl)
library(ggpubr)
#pig
df <- read_excel('E:/02-群体进化/03-猪/pig/body_weight_pig.xlsx',sheet = 2)

head(df)
p<-ggplot(df,aes(x=Male,y=Female))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm',color='blue',formula = y ~ x)+
  stat_cor(method = 'spearman')+
  theme_bw()+
  xlab('Male Pig Body Weight(kg)')+
  ylab('Female Pig Body Weight(kg)')
p
ggsave("pig_BW.png",p,dpi = 500,height = 4,width = 4)
