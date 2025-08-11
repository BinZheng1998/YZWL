setwd(dir = '~/project/01_evolution/convergent_evolution/07_result/20250805_annovar/domestication/')
chicken <- read.table('chicken/chicken.avinput.variant_function',sep = '\t')
cattle <- read.table('cattle/cattle.avinput.variant_function',sep = '\t')
sheep <- read.table('sheep/sheep.avinput.variant_function',sep = '\t')
pig <- read.table('pig/pig.avinput.variant_function',sep = '\t',fill = T)
dog <- read.table('dog/dog.avinput.variant_function',sep = '\t',fill = T)

chicken_CDS <- chicken[chicken$V1 %in% c('exonic','exonic;splicing','intronic'),]
#非编码区域包含 非编码RNA，5’和3’区，基因间区
cattle_CDS  <- cattle[cattle$V1 %in% c('exonic','exonic;splicing','intronic'),]
pig_CDS  <- pig[pig$V1 %in% c('exonic','exonic;splicing','intronic'),]
sheep_CDS  <- sheep[sheep$V1 %in% c('exonic','exonic;splicing','intronic'),]
dog_CDS  <- dog[dog$V1 %in% c('exonic','exonic;splicing','intronic'),]

chicken_exon <- chicken[chicken$V1 %in% c('exonic','exonic;splicing'),]
#非编码区域包含 非编码RNA，5’和3’区，基因间区
cattle_exon <- cattle[cattle$V1 %in% c('exonic','exonic;splicing'),]
pig_exon <- pig[pig$V1 %in% c('exonic','exonic;splicing'),]
sheep_exon <- sheep[sheep$V1 %in% c('exonic','exonic;splicing'),]
dog_exon <- dog[dog$V1 %in% c('exonic','exonic;splicing'),]

chicken_intron <- chicken[chicken$V1 %in% c('intronic'),]
#非编码区域包含 非编码RNA，5’和3’区，基因间区
cattle_intron <- cattle[cattle$V1 %in% c('intronic'),]
pig_intron <- pig[pig$V1 %in% c('intronic'),]
sheep_intron <- sheep[sheep$V1 %in% c('intronic'),]
dog_intron <- dog[dog$V1 %in% c('intronic'),]

data <- data.frame(species = c('Chicken','Cattle','Pig','Sheep','Dog',
                               'Chicken','Cattle','Pig','Sheep','Dog',
                               'Chicken','Cattle','Pig','Sheep','Dog'),
                   SNP = c(nrow(chicken_exon),nrow(cattle_exon),nrow(pig_exon),nrow(sheep_exon),nrow(dog_exon),
                           nrow(chicken_intron),nrow(cattle_intron),nrow(pig_intron),nrow(sheep_intron),nrow(dog_intron),
                           nrow(chicken)-nrow(chicken_CDS),nrow(cattle)-nrow(cattle_CDS),nrow(pig)-nrow(pig_CDS),nrow(sheep)-nrow(sheep_CDS),nrow(dog)-nrow(dog_CDS)),
                   Regions = c('Exonic','Exonic','Exonic','Exonic','Exonic','Intronic','Intronic','Intronic','Intronic','Intronic','Noncoding','Noncoding','Noncoding','Noncoding','Noncoding'))
data
data$SNP <- as.numeric(data$SNP)
library(ggplot2)
library(ggsci)
1-nrow(chicken_CDS)/nrow(chicken)
1-nrow(cattle_CDS)/nrow(cattle)
1-nrow(pig_CDS)/nrow(pig)
1-nrow(sheep_CDS)/nrow(sheep)
p <- ggplot(data, aes( x = species, weight = SNP, fill = Regions))+
  geom_bar(position = "fill")+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x='',y='SNPs Percent(%)')+
  annotate(geom = 'text',x=1,y=0.25,label=round(1-nrow(cattle_CDS)/nrow(cattle),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.25,label=round(1-nrow(chicken_CDS)/nrow(chicken),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.25,label=round(1-nrow(pig_CDS)/nrow(pig),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.25,label=round(1-nrow(sheep_CDS)/nrow(sheep),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.25,label=round(1-nrow(dog_CDS)/nrow(dog),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.75,label=round(nrow(cattle_intron)/nrow(cattle),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.75,label=round(nrow(chicken_intron)/nrow(chicken),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.75,label=round(nrow(pig_intron)/nrow(pig),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.75,label=round(nrow(sheep_intron)/nrow(sheep),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.75,label=round(nrow(dog_intron)/nrow(dog),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.975,label=round(nrow(cattle_exon)/nrow(cattle),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.975,label=round(nrow(chicken_exon)/nrow(chicken),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.975,label=round(nrow(pig_exon)/nrow(pig),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.975,label=round(nrow(sheep_exon)/nrow(sheep),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.975,label=round(nrow(dog_exon)/nrow(dog),digits = 4)*100)+
  theme(axis.text = element_text(color='black'))
p
ggsave('domestic_SNPs_regions.pdf',p,dpi = 300,width = 4.5,height = 3.5)




chicken1 <- read.table('chicken/chicken.avinput.exonic_variant_function',sep = '\t')
cattle1 <- read.table('cattle/cattle.avinput.exonic_variant_function',sep = '\t')
sheep1 <- read.table('sheep/sheep.avinput.exonic_variant_function',sep = '\t')
pig1 <- read.table('pig/pig.avinput.exonic_variant_function',sep = '\t',fill = T)
dog1 <- read.table('dog/dog.avinput.exonic_variant_function',sep = '\t',fill = T)
data1 <- data.frame(species = c('Chicken','Chicken','Chicken','Chicken','Chicken',
                                'Cattle','Cattle','Cattle','Cattle','Cattle',
                                'Pig','Pig','Pig','Pig','Pig',
                                'Sheep','Sheep','Sheep','Sheep','Sheep',
                                'Dog','Dog','Dog','Dog','Dog'),
                    Type = c('Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow'),
                    SNP = c(nrow(chicken1[chicken1$V2=='nonsynonymous SNV',]),
                            nrow(chicken1[chicken1$V2=='synonymous SNV',]),
                            nrow(chicken1[chicken1$V2=='stopgain',]),
                            nrow(chicken1[chicken1$V2=='stoploss',]),
                            nrow(chicken1[chicken1$V2=='unknown',]),
                            nrow(cattle1[cattle1$V2=='nonsynonymous SNV',]),
                            nrow(cattle1[cattle1$V2=='synonymous SNV',]),
                            nrow(cattle1[cattle1$V2=='stopgain',]),
                            nrow(cattle1[cattle1$V2=='stoploss',]),
                            nrow(cattle1[cattle1$V2=='unknown',]),
                            nrow(pig1[pig1$V2=='nonsynonymous SNV',]),
                            nrow(pig1[pig1$V2=='synonymous SNV',]),
                            nrow(pig1[pig1$V2=='stopgain',]),
                            nrow(pig1[pig1$V2=='stoploss',]),
                            nrow(pig1[pig1$V2=='unknown',]),
                            nrow(sheep1[sheep1$V2=='nonsynonymous SNV',]),
                            nrow(sheep1[sheep1$V2=='synonymous SNV',]),
                            nrow(sheep1[sheep1$V2=='stopgain',]),
                            nrow(sheep1[sheep1$V2=='stoploss',]),
                            nrow(sheep1[sheep1$V2=='unknown',]),
                            nrow(dog1[dog1$V2=='nonsynonymous SNV',]),
                            nrow(dog1[dog1$V2=='synonymous SNV',]),
                            nrow(dog1[dog1$V2=='stopgain',]),
                            nrow(dog1[dog1$V2=='stoploss',]),
                            nrow(dog1[dog1$V2=='unknown',])
                    ))
head(data1)
p1 <- ggplot(data1, aes( x = species, weight = SNP, fill = Type))+
  geom_bar(position = "fill")+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  labs(x='',y='Mutation Percent(%)')+
  annotate(geom = 'text',x=1,y=0.4,label=round(nrow(cattle1[cattle1$V2=='synonymous SNV',])/nrow(cattle1),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.4,label=round(nrow(chicken1[chicken1$V2=='synonymous SNV',])/nrow(chicken1),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.4,label=round(nrow(pig1[pig1$V2=='synonymous SNV',])/nrow(pig1),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.4,label=round(nrow(sheep1[sheep1$V2=='synonymous SNV',])/nrow(sheep1),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.4,label=round(nrow(dog1[dog1$V2=='synonymous SNV',])/nrow(sheep1),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.85,label=round(nrow(cattle1[cattle1$V2=='nonsynonymous SNV',])/nrow(cattle1),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.85,label=round(nrow(chicken1[chicken1$V2=='nonsynonymous SNV',])/nrow(chicken1),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.85,label=round(nrow(pig1[pig1$V2=='nonsynonymous SNV',])/nrow(pig1),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.85,label=round(nrow(sheep1[sheep1$V2=='nonsynonymous SNV',])/nrow(sheep1),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.85,label=round(nrow(dog1[dog1$V2=='nonsynonymous SNV',])/nrow(sheep1),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.05,label=round(nrow(cattle1[cattle1$V2=='unknown',])/nrow(cattle1),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.05,label=round(nrow(chicken1[chicken1$V2=='unknown',])/nrow(chicken1),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.05,label=round(nrow(pig1[pig1$V2=='unknown',])/nrow(pig1),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.05,label=round(nrow(sheep1[sheep1$V2=='unknown',])/nrow(sheep1),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.05,label=round(nrow(dog1[dog1$V2=='unknown',])/nrow(sheep1),digits = 4)*100)+
  theme(axis.text = element_text(color='black'))
p1
ggsave('domestic_Mutation_types.pdf',p1,dpi = 300,width = 4.5,height = 3.5)
