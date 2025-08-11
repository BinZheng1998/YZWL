setwd(dir = '~/project/01_evolution/convergent_evolution/07_result/20250805_annovar/domestication/')
library(GenomicRanges)
library(dplyr)
chicken <- read.table('chicken/chicken.avinput.exonic_variant_function',sep = '\t')
colnames(chicken)[4:6] <- c("chr", "start", "end")
chicken_annovar_gr <- makeGRangesFromDataFrame(chicken,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = TRUE)

cattle <- read.table('cattle/cattle.avinput.exonic_variant_function',sep = '\t')
colnames(cattle)[4:6] <- c("chr", "start", "end")
cattle_annovar_gr <- makeGRangesFromDataFrame(cattle,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = TRUE)

sheep <- read.table('sheep/sheep.avinput.exonic_variant_function',sep = '\t')
colnames(sheep)[4:6] <- c("chr", "start", "end")
sheep_annovar_gr <- makeGRangesFromDataFrame(sheep,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = TRUE)

pig <- read.table('pig/pig.avinput.exonic_variant_function',sep = '\t',fill = T)
colnames(pig)[4:6] <- c("chr", "start", "end")
pig_annovar_gr <- makeGRangesFromDataFrame(pig,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = TRUE)

dog <- read.table('dog/dog.avinput.exonic_variant_function',sep = '\t',fill = T)
colnames(dog)[4:6] <- c("chr", "start", "end")
dog_annovar_gr <- makeGRangesFromDataFrame(dog,seqnames.field = "chr",start.field = "start",end.field = "end",keep.extra.columns = TRUE)



#chicken pos
chicken_pos <- read.table('../../20250805_domestication_DCMS/chicken/chicken_domestication_pos_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(chicken_pos)[1:3] <- c("chr", "start", "end")
chicken_pos_bed_gr <- makeGRangesFromDataFrame(chicken_pos,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(chicken_annovar_gr, chicken_pos_bed_gr)
chicken_pos_res <- chicken[queryHits(overlaps), ]

#chicken neg
chicken_neg <- read.table('../../20250805_domestication_DCMS/chicken/chicken_domestication_neg_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(chicken_neg)[1:3] <- c("chr", "start", "end")
chicken_neg_bed_gr <- makeGRangesFromDataFrame(chicken_neg,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(chicken_annovar_gr, chicken_neg_bed_gr)
chicken_neg_res <- chicken[queryHits(overlaps), ]

#cattle pos
cattle_pos <- read.table('../../20250805_domestication_DCMS/cattle/cattle_domestication_pos_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(cattle_pos)[1:3] <- c("chr", "start", "end")
cattle_pos$chr <- gsub("chr","",cattle_pos$chr)
cattle_pos_bed_gr <- makeGRangesFromDataFrame(cattle_pos,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(cattle_annovar_gr, cattle_pos_bed_gr)
cattle_pos_res <- cattle[queryHits(overlaps), ]

#cattle neg
cattle_neg <- read.table('../../20250805_domestication_DCMS/cattle/cattle_domestication_neg_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(cattle_neg)[1:3] <- c("chr", "start", "end")
cattle_neg$chr <- gsub("chr","",cattle_neg$chr)
cattle_neg_bed_gr <- makeGRangesFromDataFrame(cattle_neg,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(cattle_annovar_gr, cattle_neg_bed_gr)
cattle_neg_res <- cattle[queryHits(overlaps), ]

#sheep pos
sheep_pos <- read.table('../../20250805_domestication_DCMS/sheep/sheep_domestication_pos_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(sheep_pos)[1:3] <- c("chr", "start", "end")
sheep_pos$chr <- gsub("chr","",sheep_pos$chr)
sheep_pos_bed_gr <- makeGRangesFromDataFrame(sheep_pos,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(sheep_annovar_gr, sheep_pos_bed_gr)
sheep_pos_res <- sheep[queryHits(overlaps), ]

#sheep neg
sheep_neg <- read.table('../../20250805_domestication_DCMS/sheep/sheep_domestication_neg_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(sheep_neg)[1:3] <- c("chr", "start", "end")
sheep_neg$chr <- gsub("chr","",sheep_neg$chr)
sheep_neg_bed_gr <- makeGRangesFromDataFrame(sheep_neg,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(sheep_annovar_gr, sheep_neg_bed_gr)
sheep_neg_res <- sheep[queryHits(overlaps), ]

#pig pos
pig_pos <- read.table('../../20250805_domestication_DCMS/pig/pig_domestication_pos_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(pig_pos)[1:3] <- c("chr", "start", "end")
pig_pos$chr <- gsub("chr","",pig_pos$chr)
pig_pos_bed_gr <- makeGRangesFromDataFrame(pig_pos,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(pig_annovar_gr, pig_pos_bed_gr)
pig_pos_res <- pig[queryHits(overlaps), ]

#pig neg
pig_neg <- read.table('../../20250805_domestication_DCMS/pig/pig_domestication_neg_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(pig_neg)[1:3] <- c("chr", "start", "end")
pig_neg$chr <- gsub("chr","",pig_neg$chr)
pig_neg_bed_gr <- makeGRangesFromDataFrame(pig_neg,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(pig_annovar_gr, pig_neg_bed_gr)
pig_neg_res <- pig[queryHits(overlaps), ]

#dog pos
dog_pos <- read.table('../../20250805_domestication_DCMS/dog/dog_domestication_pos_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(dog_pos)[1:3] <- c("chr", "start", "end")
#dog_pos$chr <- gsub("chr","",dog_pos$chr)
dog_pos_bed_gr <- makeGRangesFromDataFrame(dog_pos,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(dog_annovar_gr, dog_pos_bed_gr)
dog_pos_res <- dog[queryHits(overlaps), ]

#dog neg
dog_neg <- read.table('../../20250805_domestication_DCMS/dog/dog_domestication_neg_DCMS_top_001_merged_regions.txt',sep = '\t')
colnames(dog_neg)[1:3] <- c("chr", "start", "end")
#dog_neg$chr <- gsub("chr","",dog_neg$chr)
dog_neg_bed_gr <- makeGRangesFromDataFrame(dog_neg,seqnames.field = "chr",start.field = "start",end.field = "end")
overlaps <- findOverlaps(dog_annovar_gr, dog_neg_bed_gr)
dog_neg_res <- dog[queryHits(overlaps), ]


data1 <- data.frame(species = c('Chicken','Chicken','Chicken','Chicken','Chicken',
                                'Cattle','Cattle','Cattle','Cattle','Cattle',
                                'Pig','Pig','Pig','Pig','Pig',
                                'Sheep','Sheep','Sheep','Sheep','Sheep',
                                'Dog','Dog','Dog','Dog','Dog',
                                'Chicken','Chicken','Chicken','Chicken','Chicken',
                                'Cattle','Cattle','Cattle','Cattle','Cattle',
                                'Pig','Pig','Pig','Pig','Pig',
                                'Sheep','Sheep','Sheep','Sheep','Sheep',
                                'Dog','Dog','Dog','Dog','Dog'),
                    selection = c("pos","pos",'pos','pos','pos',
                                  "pos","pos",'pos','pos','pos',
                                  "pos","pos",'pos','pos','pos',
                                  "pos","pos",'pos','pos','pos',
                                  "pos","pos",'pos','pos','pos',
                                  'neg','neg','neg','neg','neg',
                                  'neg','neg','neg','neg','neg',
                                  'neg','neg','neg','neg','neg',
                                  'neg','neg','neg','neg','neg',
                                  'neg','neg','neg','neg','neg'
                    ),
                    Type = c('Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow',
                             'Nonsynonymous','Synonymous','Stopgain','Stoploss','Unknow'),
                    SNP = c(nrow(chicken_pos_res[chicken_pos_res$V2=='nonsynonymous SNV',]),
                            nrow(chicken_pos_res[chicken_pos_res$V2=='synonymous SNV',]),
                            nrow(chicken_pos_res[chicken_pos_res$V2=='stopgain',]),
                            nrow(chicken_pos_res[chicken_pos_res$V2=='stoploss',]),
                            nrow(chicken_pos_res[chicken_pos_res$V2=="unknown",]),
                            nrow(cattle_pos_res[cattle_pos_res$V2=='nonsynonymous SNV',]),
                            nrow(cattle_pos_res[cattle_pos_res$V2=='synonymous SNV',]),
                            nrow(cattle_pos_res[cattle_pos_res$V2=='stopgain',]),
                            nrow(cattle_pos_res[cattle_pos_res$V2=='stoploss',]),
                            nrow(cattle_pos_res[cattle_pos_res$V2=='unknown',]),
                            nrow(pig_pos_res[pig_pos_res$V2=='nonsynonymous SNV',]),
                            nrow(pig_pos_res[pig_pos_res$V2=='synonymous SNV',]),
                            nrow(pig_pos_res[pig_pos_res$V2=='stopgain',]),
                            nrow(pig_pos_res[pig_pos_res$V2=='stoploss',]),
                            nrow(pig_pos_res[pig_pos_res$V2=='unknown',]),
                            nrow(sheep_pos_res[sheep_pos_res$V2=='nonsynonymous SNV',]),
                            nrow(sheep_pos_res[sheep_pos_res$V2=='synonymous SNV',]),
                            nrow(sheep_pos_res[sheep_pos_res$V2=='stopgain',]),
                            nrow(sheep_pos_res[sheep_pos_res$V2=='stoploss',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='unknown',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='nonsynonymous SNV',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='synonymous SNV',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='stopgain',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='stoploss',]),
                            nrow(dog_pos_res[dog_pos_res$V2=='unknown',]),
                            ############
                            nrow(chicken_neg_res[chicken_neg_res$V2=='nonsynonymous SNV',]),
                            nrow(chicken_neg_res[chicken_neg_res$V2=='synonymous SNV',]),
                            nrow(chicken_neg_res[chicken_neg_res$V2=='stopgain',]),
                            nrow(chicken_neg_res[chicken_neg_res$V2=='stoploss',]),
                            nrow(chicken_neg_res[chicken_neg_res$V2=='unknown',]),
                            nrow(cattle_neg_res[cattle_neg_res$V2=='nonsynonymous SNV',]),
                            nrow(cattle_neg_res[cattle_neg_res$V2=='synonymous SNV',]),
                            nrow(cattle_neg_res[cattle_neg_res$V2=='stopgain',]),
                            nrow(cattle_neg_res[cattle_neg_res$V2=='stoploss',]),
                            nrow(cattle_neg_res[cattle_neg_res$V2=='unknown',]),
                            nrow(pig_neg_res[pig_neg_res$V2=='nonsynonymous SNV',]),
                            nrow(pig_neg_res[pig_neg_res$V2=='synonymous SNV',]),
                            nrow(pig_neg_res[pig_neg_res$V2=='stopgain',]),
                            nrow(pig_neg_res[pig_neg_res$V2=='stoploss',]),
                            nrow(pig_neg_res[pig_neg_res$V2=='unknown',]),
                            nrow(sheep_neg_res[sheep_neg_res$V2=='nonsynonymous SNV',]),
                            nrow(sheep_neg_res[sheep_neg_res$V2=='synonymous SNV',]),
                            nrow(sheep_neg_res[sheep_neg_res$V2=='stopgain',]),
                            nrow(sheep_neg_res[sheep_neg_res$V2=='stoploss',]),
                            nrow(sheep_neg_res[sheep_neg_res$V2=='unknown',]),
                            nrow(dog_neg_res[dog_neg_res$V2=='nonsynonymous SNV',]),
                            nrow(dog_neg_res[dog_neg_res$V2=='synonymous SNV',]),
                            nrow(dog_neg_res[dog_neg_res$V2=='stopgain',]),
                            nrow(dog_neg_res[dog_neg_res$V2=='stoploss',]),
                            nrow(dog_neg_res[dog_neg_res$V2=='unknown',])
                    ))

data1
data1$SNP <- as.numeric(data1$SNP)
library(ggplot2)
library(ggsci)
pos_data <- data1[data1$selection == "pos",]
p1 <- ggplot(pos_data, aes( x = species, weight = SNP, fill = Type))+
  geom_bar(position = "fill")+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0),breaks = c(0,0.25,0.50,0.75,1),labels = c("0","25","50","75","100"))+
  theme_classic()+
  labs(x='',y='POS SNPs Percent(%)')+
  annotate(geom = 'text',x=1,y=0.25,label=round(nrow(cattle_pos_res[cattle_pos_res$V2 %in% c('synonymous SNV'),])/nrow(cattle_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.25,label=round(nrow(chicken_pos_res[chicken_pos_res$V2 %in% c('synonymous SNV'),])/nrow(chicken_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.25,label=round(nrow(pig_pos_res[pig_pos_res$V2 %in% c('synonymous SNV'),])/nrow(pig_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.25,label=round(nrow(sheep_pos_res[sheep_pos_res$V2 %in% c('synonymous SNV'),])/nrow(sheep_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.25,label=round(nrow(dog_pos_res[dog_pos_res$V2 %in% c('synonymous SNV'),])/nrow(dog_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.85,label=round(nrow(cattle_pos_res[cattle_pos_res$V2 %in% c('nonsynonymous SNV'),])/nrow(cattle_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.85,label=round(nrow(chicken_pos_res[chicken_pos_res$V2 %in% c('nonsynonymous SNV'),])/nrow(chicken_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.85,label=round(nrow(pig_pos_res[pig_pos_res$V2 %in% c('nonsynonymous SNV'),])/nrow(pig_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.85,label=round(nrow(sheep_pos_res[sheep_pos_res$V2 %in% c('nonsynonymous SNV'),])/nrow(sheep_pos_res),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.85,label=round(nrow(dog_pos_res[dog_pos_res$V2 %in% c('nonsynonymous SNV'),])/nrow(dog_pos_res),digits = 4)*100)+
  coord_flip()+
  theme(axis.text = element_text(color='black',size=12),legend.position = 'none')
p1
neg_data <- data1[data1$selection == "neg",]
p2 <- ggplot(neg_data, aes( x = species, weight = SNP, fill = Type))+
  geom_bar(position = "fill")+
  scale_fill_npg()+
  scale_y_continuous(expand = c(0,0),breaks = c(0,0.25,0.50,0.75,1),labels = c("0","25","50","75","100"))+
  theme_classic()+
  labs(x='',y='NEG SNPs Percent(%)')+
  annotate(geom = 'text',x=1,y=0.25,label=round(nrow(cattle_neg_res[cattle_neg_res$V2 %in% c('synonymous SNV'),])/nrow(cattle_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.25,label=round(nrow(chicken_neg_res[chicken_neg_res$V2 %in% c('synonymous SNV'),])/nrow(chicken_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.25,label=round(nrow(pig_neg_res[pig_neg_res$V2 %in% c('synonymous SNV'),])/nrow(pig_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.25,label=round(nrow(sheep_neg_res[sheep_neg_res$V2 %in% c('synonymous SNV'),])/nrow(sheep_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.25,label=round(nrow(dog_neg_res[dog_neg_res$V2 %in% c('synonymous SNV'),])/nrow(dog_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=1,y=0.85,label=round(nrow(cattle_neg_res[cattle_neg_res$V2 %in% c('nonsynonymous SNV'),])/nrow(cattle_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=2,y=0.85,label=round(nrow(chicken_neg_res[chicken_neg_res$V2 %in% c('nonsynonymous SNV'),])/nrow(chicken_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=4,y=0.85,label=round(nrow(pig_neg_res[pig_neg_res$V2 %in% c('nonsynonymous SNV'),])/nrow(pig_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=5,y=0.85,label=round(nrow(sheep_neg_res[sheep_neg_res$V2 %in% c('nonsynonymous SNV'),])/nrow(sheep_neg_res),digits = 4)*100)+
  annotate(geom = 'text',x=3,y=0.85,label=round(nrow(dog_neg_res[dog_neg_res$V2 %in% c('nonsynonymous SNV'),])/nrow(dog_neg_res),digits = 4)*100)+
  
  coord_flip()+
  theme(axis.text = element_text(color='black',size=12),
        #legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y=element_blank())
p2
p<-p1+p2
p
ggsave('domestication_selection_regions_Mutation_types.pdf',p,dpi = 300,width = 6.5,height = 3)
