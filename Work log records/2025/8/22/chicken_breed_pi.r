setwd('~/project/01_evolution/convergent_evolution/07_result/20250813_chicken_breed_piANDallele_fequency/')
library(dplyr)
data <- read_tsv('res/pi_means_per_chromosome.tsv')
#data1 <- t(data)
#colnames(data1) <- data1[1,]
#data1 <- data1[-1,]
long_data <- data %>%
  pivot_longer(
    cols = -CHROM,          # 选择除chr列外的所有列
    names_to = "species", # 新列名：存储原列名（物种名）
    values_to = "value"   # 新列名：存储原单元格值
  )
long_data <- long_data %>%
  mutate(type = case_when(
    species == "Gallus_gallus_gallus" ~ "RJF",
    species %in% c("Gallus_gallus_bankiva","Gallus_gallus_murghi","Gallus_gallus_spadiceus") ~ "subRJF",
    species %in% c("Gallus_varius","Gallus_sonneratii", "Gallus_lafayettii") ~ "otherJF",
    TRUE ~ "Domestic"  # 剩余所有情况
  ))
long_data$CHROM <- gsub('chr','Chr',long_data$CHROM)
chromosome <- paste0('Chr',1:38)
long_data$CHROM <- factor(long_data$CHROM,levels = chromosome)
long_data_sorted <- long_data %>%
  arrange(factor(type, levels = c( "Domestic","otherJF", "RJF")))

p<-ggplot()+
  geom_point(long_data_sorted,mapping=aes(y=CHROM,x=value,color=type))+
  scale_color_manual(values = c('RJF'='red','otherJF'='blue','subRJF'='green','Domestic'='grey90'))+
  theme(panel.background = element_blank())+
  labs(y='',x='PI')
p
ggsave('chicken_breeds_pi.pdf',height = 8,width = 14,dpi = 300)


data1 <- read_tsv('RJF_res/pi_means_per_chromosome.tsv')
long_data1 <- data1 %>%
  pivot_longer(
    cols = -CHROM,          # 选择除chr列外的所有列
    names_to = "species", # 新列名：存储原列名（物种名）
    values_to = "value"   # 新列名：存储原单元格值
  )
long_data1$CHROM <- gsub('chr','Chr',long_data1$CHROM)
long_data1$CHROM <- factor(long_data1$CHROM,levels = chromosome)
p1<-ggplot()+
  geom_point(long_data1,mapping=aes(y=CHROM,x=value,color=species))+
  scale_color_manual(values = c('RJF_Indonesia'='red','RJF_China'='blue','RJF_Thailand'='green','RJF_Singapore'='grey90'))+
  theme(panel.background = element_blank())+
  labs(y='',x='PI')
p1
ggsave('RJF_breeds_pi.pdf',height = 8,width = 14,dpi = 300)
