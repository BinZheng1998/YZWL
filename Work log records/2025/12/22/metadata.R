metadata <- read.table('../pig/1.data/pigGTEx_metadata.txt',sep = '\t',header = T)
df2 <- metadata[, c(2, 9, 19)]
colnames(df2) <- c('Sample', 'BW', 'Tissue')

df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose","Pituitary"),]
df3 <- df2[df2$BW %in% c("high","low"),]
df4 <- df3[df3$Tissue == "Pituitary",]
table(df4$BW)

metadata <- read.table('../cattle/metadata/cattle_metadata.txt',sep = '\t',header = T)
df2 <- metadata[, c(1, 5, 6)]
colnames(df2) <- c('Sample', 'BW', 'Tissue')

df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose","Pituitary"),]
df3 <- df2[df2$BW %in% c("high","low"),]
df4 <- df3[df3$Tissue == "Muscle",]
table(df4$BW)

metadata <- read.table('../chicken/metadata/chicken_rnaseq_metadata.txt',sep = '\t',header = T,fill=T)
df2 <- metadata[, c(2, 11, 14)]
colnames(df2) <- c('Sample', 'BW', 'Tissue')

df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose","Pituitary"),]
df3 <- df2[df2$BW %in% c("High","Low"),]
df4 <- df3[df3$Tissue == "Adipose",]
table(df4$BW)


metadata <- read.table('../sheep/1.data/metadata/sheepGTEx_metadata.txt',sep = '\t',header = T,fill=T)
df2 <- metadata[, c(1, 14, 36)]
colnames(df2) <- c('Sample', 'BW', 'Tissue')

df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose","Pituitary"),]
df3 <- df2[df2$BW %in% c("High","Low"),]
df4 <- df3[df3$Tissue == "Pituitary",]
table(df4$BW)
