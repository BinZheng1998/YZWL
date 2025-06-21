setwd(dir = "~/project/01_evolution/convergent_evolution/07_result/20250616_chr_regions_plot/")
library(RIdeogram)

#chicken
chicken <- read.table('chicken_chr_length.txt',sep = '\t',header = T)
ideogram(karyotype = chicken)
df1 <- read.table("../20250605_DCMS_res/convergent/chicken/chicken_convergent_pos_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df1) <- c("Chr","Start","End")
df1$Chr <- gsub("chr","",df1$Chr)
df1$Value <- 1
df2 <- read.table("../20250605_DCMS_res/convergent/chicken/chicken_convergent_neg_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df2) <- c("Chr","Start","End")
df2$Chr <- gsub("chr","",df2$Chr)
df2$Value <- 1

ideogram(karyotype = chicken,overlaid = df1,label = df2,
         label_type = "heatmap",colorset1 = c("white","#e34a33"),colorset2 = c("white","#2c7fb8"),
         output = 'chicken_convergent_top001.svg')

#cattle
cattle <- read.table('cattle_chr_length.txt',sep = '\t',header = T)
#ideogram(karyotype = cattle)
df1 <- read.table("../20250605_DCMS_res/convergent/cattle/cattle_convergent_pos_DCMS_top_005_merged_regions.txt",sep = '\t')
colnames(df1) <- c("Chr","Start","End")
df1$Chr <- gsub("chr","",df1$Chr)
df1$Value <- 1
df2 <- read.table("../20250605_DCMS_res/convergent/cattle/cattle_convergent_neg_DCMS_top_005_merged_regions.txt",sep = '\t')
colnames(df2) <- c("Chr","Start","End")
df2$Chr <- gsub("chr","",df2$Chr)
df2$Value <- 1

ideogram(karyotype = cattle,overlaid = df1,label = df2,
         label_type = "heatmap",colorset1 = c("white","#e34a33"),colorset2 = c("white","#2c7fb8"),
         output = 'cattle_convergent_top005.svg')

#pig
pig <- read.table('pig_chr_length.txt',sep = '\t',header = T)
#ideogram(karyotype = cattle)
df1 <- read.table("../20250605_DCMS_res/convergent/pig/pig_convergent_pos_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df1) <- c("Chr","Start","End")
df1$Chr <- gsub("chr","",df1$Chr)
df1$Value <- 1
df2 <- read.table("../20250605_DCMS_res/convergent/pig/pig_convergent_neg_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df2) <- c("Chr","Start","End")
df2$Chr <- gsub("chr","",df2$Chr)
df2$Value <- 1

ideogram(karyotype = pig,overlaid = df1,label = df2,
         label_type = "heatmap",colorset1 = c("white","#e34a33"),colorset2 = c("white","#2c7fb8"),
         output = 'pig_convergent_top001.svg')


#sheep
sheep <- read.table('sheep_chr_length.txt',sep = '\t',header = T)
#ideogram(karyotype = cattle)
df1 <- read.table("../20250605_DCMS_res/convergent/sheep/sheep_convergent_pos_DCMS_top_005_merged_regions.txt",sep = '\t')
colnames(df1) <- c("Chr","Start","End")
df1$Chr <- gsub("chr","",df1$Chr)
df1$Value <- 1
df2 <- read.table("../20250605_DCMS_res/convergent/sheep/sheep_convergent_neg_DCMS_top_005_merged_regions.txt",sep = '\t')
colnames(df2) <- c("Chr","Start","End")
df2$Chr <- gsub("chr","",df2$Chr)
df2$Value <- 1

ideogram(karyotype = sheep,overlaid = df1,label = df2,
         label_type = "heatmap",colorset1 = c("white","#e34a33"),colorset2 = c("white","#2c7fb8"),
         output = 'sheep_convergent_top005.svg')


#dog
dog <- read.table('dog_chr_length.txt',sep = '\t',header = T)
#ideogram(karyotype = cattle)
df1 <- read.table("../20250605_DCMS_res/convergent/dog/dog_convergent_pos_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df1) <- c("Chr","Start","End")
df1$Chr <- gsub("chr","",df1$Chr)
df1$Value <- 1
df2 <- read.table("../20250605_DCMS_res/convergent/dog/dog_convergent_neg_DCMS_top_001_merged_regions.txt",sep = '\t')
colnames(df2) <- c("Chr","Start","End")
df2$Chr <- gsub("chr","",df2$Chr)
df2$Value <- 1

ideogram(karyotype = dog,overlaid = df1,label = df2,
         label_type = "heatmap",colorset1 = c("white","#e34a33"),colorset2 = c("white","#2c7fb8"),
         output = 'dog_convergent_top001.svg')
