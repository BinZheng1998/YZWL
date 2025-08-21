setwd(dir = '~/project/01_evolution/convergent_evolution/07_result/20250805_annovar/domestication/')
library(dplyr)
library(tidyr)
library(stringr)
library(aplot)
library(funkyheatmap)
library(ggsci)
data <- read.table('chicken/chicken.avinput.variant_function',sep = '\t')
longdata1 <- data %>%
  group_by(V3) %>%
  summarise(
    total_snp =n(),
    nonCDS_SNP = sum(!V1 %in% c('exonic', 'exonic;splicing', 'intronic')),
    Intron_SNP = sum(V1 == 'intronic'),
    Exon_SNP = sum(V1 %in% c('exonic', 'exonic;splicing')),
    nonCDS_ratio = round(nonCDS_SNP/total_snp,4),
    intron_ratio = round(Intron_SNP/total_snp,4),
    exon_ratio = round(Exon_SNP/total_snp,4)
  ) %>%
  rename(Chromosome = V3)
longdata1$Chromosome <- gsub('chr','Chr',longdata1$Chromosome)
longdata1 <- longdata1 %>%pivot_longer(cols = c('nonCDS_SNP','Exon_SNP','Intron_SNP'),names_to = 'Region',values_to = 'SNP_count')
longdata1 <- longdata1[!longdata1$Chromosome %in% c('ChrZ','ChrW'),]
chromosome <- paste0("Chr",1:38)
longdata1$Chromosome <- factor(longdata1$Chromosome,levels = chromosome)
longdata1 <- longdata1[order(longdata1$Chromosome), ]
exon1 <- longdata1[longdata1$Region=="Exon_SNP",]
exon1 <- exon1[,c(1,5,6)]
intron1<-longdata1[longdata1$Region=="Intron_SNP",]
intron1 <- intron1[,c(1,4,6)]
colnames(intron1)[2] <- 'ratio'
nocds1 <- longdata1[longdata1$Region=="nonCDS_SNP",]
nocds1 <- nocds1[,c(1,3,6)]
colnames(nocds1)[2] <- 'ratio'
data1 <- rbind(intron1,nocds1)

#
data2 <- read.table('chicken/chicken.avinput.exonic_variant_function',sep = '\t')
longdata2 <- data2 %>%
  group_by(V4) %>%
  summarise(
    total_snp =n(),
    Nonsynonymous = sum(V2 %in% c('nonsynonymous SNV','stopgain','stoploss')),
    Synonymous = sum(V2 %in% c('synonymous SNV')),
    #Stopgain = sum(V2 == 'stopgain'),
    #Stoploss = sum(V2 == 'stoploss'),
    Unknown = sum(V2 == 'unknown'),
    Nonsynonymous_ratio = round(Nonsynonymous/total_snp,4),
    Synonymous_ratio = round(Synonymous/total_snp,4),
    Unknown_ratio = round(Unknown/total_snp,4),
  ) %>%
  rename(Chromosome = V4)
longdata2$Chromosome <- gsub('chr','Chr',longdata2$Chromosome)

longdata2 <-longdata2 %>%pivot_longer(cols = c('Nonsynonymous','Synonymous','Unknown'),names_to = 'Type', values_to = 'SNP_count')
longdata2 <- longdata2[!longdata2$Chromosome %in% c('ChrZ','ChrW'),]
longdata2$Chromosome <- factor(longdata2$Chromosome,levels = chromosome)
longdata2 <- longdata2[order(longdata2$Chromosome), ]

Nonsynonymous1 <- longdata2[longdata2$Type=="Nonsynonymous",]
Nonsynonymous1 <- Nonsynonymous1[,c(1,3,6)]
colnames(Nonsynonymous1)[2] <- 'ratio'
Nonsynonymous1$exon_ratio <-  exon1$exon_ratio
Nonsynonymous1$ratio2 <- Nonsynonymous1$exon_ratio * Nonsynonymous1$ratio
  
synonymous1 <- longdata2[longdata2$Type=="Synonymous",]
synonymous1 <- synonymous1[,c(1,4,6)]
colnames(synonymous1)[2] <- 'ratio'
synonymous1$exon_ratio <-  exon1$exon_ratio
synonymous1$ratio2 <- synonymous1$exon_ratio * synonymous1$ratio

unknown1 <- longdata2[longdata2$Type=="Unknown",]
unknown1 <- unknown1[,c(1,5,6)]
colnames(unknown1)[2] <- 'ratio'
unknown1$exon_ratio <-  exon1$exon_ratio
unknown1$ratio2 <- unknown1$exon_ratio * unknown1$ratio
data2 <- rbind(unknown1,Nonsynonymous1,synonymous1)
data2 <- data2[,c(1,5,3)]

head(data1)
colnames(data1)<-c('Chromosome','ratio','Type')
data1$ratio <- data1$ratio*-1
head(data1)
data2$ratio <- data2$ratio2*5
data2 <- data2[c(1,4,3)]
head(data2)
data <- rbind(data1,data2)

data$Type <- factor(data$Type,levels = c('Unknown','Synonymous','Nonsynonymous','nonCDS_SNP','Intron_SNP'))
chromosome <- paste0("Chr",1:38)
data$Chromosome <- factor(data$Chromosome,levels = rev(chromosome))

data_lable1 <- data[data$Type %in% c('Unknown','Synonymous','Nonsynonymous'),]
data_lable1$ratio <- 100* data_lable1$ratio

data_lable2 <- data[data$Type %in% c('nonCDS_SNP','Intron_SNP'),]
data_lable2$ratio <- -100* data_lable2$ratio
data$ratio<- 100* data$ratio

p <- ggplot(data, aes(y = Chromosome, x = ratio, fill = Type)) +
  geom_col(position = "stack", width = 0.9) +
  labs(x = '', y = '') +
  scale_x_continuous(
    breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100),
    labels = c('100%', '75%', '50%', '25%', 0, '5%', '10%', '15%', "20%"),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c(
    'Unknown' = alpha('#3c5488ff', 0.2),
    'Synonymous' = alpha('#3c5488ff', 0.6),
    'Nonsynonymous' = alpha('#3c5488ff', 1),
    'Intron_SNP' = alpha('#e64b35ff', 1),
    'nonCDS_SNP' = alpha('#e64b35ff', 0.6)),name='') +
  geom_text(data = data_lable1,
    aes(y = Chromosome, x = ratio,
        label = sprintf("%0.0f", round(ratio, digits = 1))),
    position = position_stack(vjust = 0.05),  # 关键修改
    color = "black",  # 建议使用对比色
    hjust=0,
    size = 3
  ) +
  geom_text(data = data_lable2,
            aes(y = Chromosome, x = -ratio,
                label = sprintf("%0.0f", round(ratio, digits = 1))),
            position = position_stack(vjust = 0.95),  # 关键修改
            color = "black",  # 建议使用对比色
            hjust=1,
            size = 3
  ) +
  theme(legend.position = 'top',
    axis.ticks = element_blank(),
    axis.text = element_text(color = 'black'),
    axis.text.y = element_text(hjust = 0),
    panel.background = element_blank()
  )
p
ggsave('chicken_chr_SNPs_20250820.pdf',p,height = 12,width = 7,dpi = 300)
