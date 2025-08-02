setwd('~/project/01_evolution/convergent_evolution/07_result/20250721_freq/')
wild <- read.table('wild_chr29_freq.frq',sep = '\t')
wild$type <- 'wild'
domestic <- read.table('domestic_chr29_freq.frq',sep = '\t')
domestic$type <- 'domestic'
data <- rbind(wild,domestic)
data$V5 <- gsub("^[A-Z]:","",data$V5)
data$V6 <- gsub("^[A-Z]:","",data$V6)
data$V5 <- as.numeric(data$V5)
data$V6 <- as.numeric(data$V6)
t.test(V6 ~ type,data=data)



results <- data.frame(
  chromosome = character(),
  p_value = numeric(),
  t_statistic = numeric(),
  df = numeric(),
  wild_mean = numeric(),
  domestic_mean = numeric(),
  stringsAsFactors = FALSE
)

for (chr in 1:38) {
  wild_file <- paste0("wild_chr", chr, "_freq.frq")
  domestic_file <- paste0("domestic_chr", chr, "_freq.frq")

  if (!file.exists(wild_file) | !file.exists(domestic_file)) {
    warning(paste("Files for chr", chr, "not found. Skipping..."))
    next
  }
  wild <- read.table(wild_file, sep = "\t", header = FALSE)
  domestic <- read.table(domestic_file, sep = "\t", header = FALSE)
  wild$type <- "wild"
  domestic$type <- "domestic"
  data <- rbind(wild, domestic)
  data$V5 <- as.numeric(gsub("^[A-Z]:", "", data$V5))  # 假设V5是等位基因频率列
  data$V6 <- as.numeric(gsub("^[A-Z]:", "", data$V6))  # 假设V6是次要等位基因频率
  t_test <- t.test(V6 ~ type, data = data)
  p_value <- t_test$p.value
  t_statistic <- t_test$statistic
  df <- t_test$parameter
  wild_mean <- mean(data$V6[data$type == "wild"], na.rm = TRUE)
  domestic_mean <- mean(data$V6[data$type == "domestic"], na.rm = TRUE)
  
  results <- rbind(results, data.frame(
    chromosome = paste0("chr", chr),
    p_value = p_value,
    t_statistic = t_statistic,
    df = df,
    wild_mean = wild_mean,
    domestic_mean = domestic_mean
  ))
}
write.csv(results, "wild_vs_domestic_t_test_results.csv", row.names = FALSE)

results$FC <- results$domestic_mean/results$wild_mean
results$chromosome <- factor(results$chromosome,levels = unique(results$chromosome))
p<-ggplot(results,aes(x= chromosome,y= FC))+
  geom_segment(aes(x= chromosome,xend = chromosome,y= 0,yend= FC))+
  geom_point(size=4,color='grey')+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.5))+
  labs(x='',y='ALT Allele Freqency Fold Change')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,color = 'black',size = 10),
        axis.text.y = element_text(color='black',size=10),
        panel.border=element_rect(size=1),
        panel.grid = element_blank())
p
ggsave('ALT_AF_foldChange.pdf',p,dpi = 300,width = 8,height = 3.5)

res1 <- results
res1$type <- 'wild'
res1$af <- res1$wild_mean
res2 <- results
res2$type <- 'domestic'
res2$af <- res2$domestic_mean
res <- rbind(res1,res2)
res3 <- res
p<-ggplot(res,aes(chromosome,af))+
  geom_col(aes(fill=type),position = 'dodge')+
  scale_fill_manual(values = c('#E64b35b2','#4dbbd5b2'))+
  scale_y_continuous(limits = c(0,0.23),expand = c(0,0),sec.axis = sec_axis(~.*7.5,name = 'ALT Allele Frequency Ratio'))+    
  geom_point(res3,mapping=aes(chromosome,FC/7.5),size=3,color='grey')+
  geom_line(res3,mapping=aes(chromosome,FC/7.5,group=type),color='grey')+
  theme_bw()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 60,hjust = 1),
        panel.grid = element_blank())+
  labs(x='',y='ALT Allele Frequency')+
  geom_rect(aes(xmin = 30,xmax = 31,ymin = 0.21,ymax = 0.22),fill='#4dbbd5b2')+
  annotate(geom = 'text',x=33,y=0.215,label='Wild',size=4)+  
  geom_rect(aes(xmin = 30,xmax = 31,ymin = 0.195,ymax = 0.205),fill='#E64b35b2')+
  annotate(geom = 'text',x=33,y=0.20,label='Domestic',size=4)+
  geom_point(aes(x = 30.5,y=0.185),color='grey',size=4)+
  annotate(geom = 'text',x=34,y=0.185,label='ALT AF Ratio',size=4)
p  
ggsave('ALT_AF_foldChange2.pdf',p,dpi = 300,width = 9,height = 4)
