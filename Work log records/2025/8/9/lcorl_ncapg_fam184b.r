setwd('~/project/01_evolution/convergent_evolution/07_result/20250808_selection_region_plot/BW/lcorl-ncapg-fam184b/')
chicken <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/2.snp.site/fst/chicken_BW.weir.fst',sep = '\t',header = T)
chicken1 <- chicken[chicken$CHROM == "chr4",]
# FAM184B 75567781 75615607
min_site<-75567781-250000
max_site <- 75615607+250000
chicken2 <- chicken1[chicken1$POS > min_site & chicken1$POS < max_site,]
head(chicken2)
p1<-ggplot(chicken2,aes(POS,WEIR_AND_COCKERHAM_FST))+
  geom_point(color='grey',size=1)+
  theme_classic()+
  theme(axis.title.y = element_text(face = "italic"))+
  scale_x_continuous(
    breaks = c(75400000, 75500000, 75600000, 75700000,75800000),
    labels = c( "75.4", "75.5", "75.6", "75.7","75.8"))+
  labs(y='Fst',x='Chromosome 4 position(Mb)')+
  annotate("segment", x = 75567781, y = 0.24,xend = 75615607, yend = 0.24,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=75567781-50000,y=0.24,label='FAM184B',size=3,fontface='italic')+
  annotate("segment", x = 75468951, y = 0.26,xend = 75541506, yend = 0.26,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=75468951-50000,y=0.26,label='LCORL',size=3,fontface='italic')+
  annotate("segment", x = 75567728, y = 0.26,xend = 75544703, yend = 0.26,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=75567728+50000,y=0.26,label='NCAPG',size=3,fontface='italic')
p1
ggsave('chicken.pdf',p1,dpi = 300,width = 4,height = 3)
#cattle
cattle <- read.table('~/project/01_evolution/convergent_evolution/05_cattle/selection_BW/DP4_2allele_50miss_maf005/2.snp.site/fst/cattle_BW.weir.fst',sep = '\t',header = T)
cattle1 <- cattle[cattle$CHROM == "6",]
# FAM184B 38824539 38945766
min_site<-38824539-250000
max_site <- 38945766+250000
cattle2 <- cattle1[cattle1$POS > min_site & cattle1$POS < max_site,]
head(cattle2)
p2<-ggplot(cattle2,aes(POS,WEIR_AND_COCKERHAM_FST))+
  geom_point(color='grey',size=1)+
  theme_classic()+
  theme(axis.title.y = element_text(face = "italic"))+
  scale_x_continuous(
    breaks = c(38500000, 38600000, 38700000, 38800000,38900000,39000000,39100000),
    labels = c( "38.5", "38.6", "38.7", "38.8","38.9", "39.0","39.1"))+
  labs(y='Fst',x='Chromosome 6 position(Mb)')+
  annotate("segment", x = 38945766, y = 0.8,xend = 38824539, yend = 0.8,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=38945766+70000,y=0.8,label='FAM184B',size=3,fontface='italic')+
  annotate("segment", x = 39201942, y = 0.75,xend = 39051450, yend = 0.75,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=39201942-75000,y=0.72,label='LCORL',size=3,fontface='italic')+
  annotate("segment", x = 38976555, y = 0.75,xend = 39022642, yend = 0.75,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=38976555+25000,y=0.72,label='NCAPG',size=3,fontface='italic')
p2
ggsave('cattle.pdf',p2,dpi = 300,width = 4,height = 3)

#sheep
sheep <- read.table('~/project/01_evolution/convergent_evolution/02_sheep/selection_analysis_BW/DP4_2allele_50miss_maf005/2.snp.site/sheep_BW.weir.fst',sep = '\t',header = T)
sheep1 <- sheep[sheep$CHROM == "6",]
# FAM184B 41940280 42059944
min_site<-41940280-250000
max_site <- 42059944+250000
sheep2 <- sheep1[sheep1$POS > min_site & sheep1$POS < max_site,]
head(sheep2)
p3<-ggplot(sheep2,aes(POS,WEIR_AND_COCKERHAM_FST))+
  geom_point(color='grey',size=1)+
  theme_classic()+
  theme(axis.title.y = element_text(face = "italic"))+
  scale_x_continuous(
    breaks = c(41800000, 41900000, 42000000, 42100000,42200000,42300000),
    labels = c( "41.8", "41.9", "42.0", "42.1","42.2","42.3"))+
  labs(y='Fst',x='Chromosome 6 position(Mb)')+
  annotate("segment", x = 42059944, y = 0.35,xend = 41940280, yend = 0.35,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=42059944+70000,y=0.35,label='FAM184B',size=3,fontface='italic')+
  annotate("segment", x = 42306535, y = 0.33,xend = 42136396, yend = 0.33,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=42306535-75000,y=0.31,label='LCORL',size=3,fontface='italic')+
  annotate("segment", x = 42091768, y = 0.33,xend = 42136696, yend = 0.33,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=42091768+25000,y=0.31,label='NCAPG',size=3,fontface='italic')
p3
ggsave('sheep.pdf',p3,dpi = 300,width = 4,height = 3)
#pig
pig <- read.table('~/project/01_evolution/convergent_evolution/03_pig/selection_analysis_BW/DP4_2allele_50miss_maf005/2.snp.site/fst/pig_BW.weir.fst',sep = '\t',header = T)
pig1 <- pig[pig$CHROM == "8",]
# FAM184B 12635359 12741482
min_site<-12635359-1000000
max_site <- 12741482+1000000
pig2 <- pig1[pig1$POS > min_site & pig1$POS < max_site,]
head(pig2)
p4<-ggplot(pig2,aes(POS,WEIR_AND_COCKERHAM_FST))+
  geom_point(color='grey',size=1)+
  theme_classic()+
  theme(axis.title.y = element_text(face = "italic"))+
  scale_x_continuous(
    breaks = c(12400000, 12500000, 12600000, 12700000,12800000,12900000),
    labels = c( "12.4", "12.5", "12.6", "12.7","12.8","12.9"))+
  labs(y='Fst',x='Chromosome 8 position(Mb)')+
  annotate("segment", x = 12741482, y = 0.92,xend = 12635359, yend = 0.92,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=12741482+70000,y=0.92,label='FAM184B',size=3,fontface='italic')+
  annotate("segment", x = 12969370, y = 0.87,xend = 12806878, yend = 0.87,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=12969370-75000,y=0.82,label='LCORL',size=3,fontface='italic')+
  annotate("segment", x = 12759225, y = 0.87,xend = 12807164, yend = 0.87,arrow = arrow( type = "closed",length = unit(0.05, "inches")),color = "grey40",linewidth = 0.35)+
  annotate('text',x=12759225+25000,y=0.82,label='NCAPG',size=3,fontface='italic')
p4
#ggsave('pig2.pdf',p4,dpi = 300,width = 4,height = 3)
ggsave('pig2.pdf',p4,dpi = 300,width = 8,height = 4.5)
