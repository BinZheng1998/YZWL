setwd('~/project/01_evolution/convergent_evolution/07_result/20250718_small_chr_result/')
repeatGC <- read.table('repeat_GC.txt',sep = '\t',header = T)
chr29_repeat <- repeatGC[repeatGC$Chromosome == "chr29",]
colnames(chr29_repeat) <- c('chr','start','end','lowercase','gc')

chr29_snp <- read.table('chr29_snps_5kb_stats.txt',sep = '\t',header = T)
colnames(chr29_snp) <- c('chr','start','end','snp')
chr29_snp$end <- chr29_snp$end+1

chr29_indel <- read.table('chr29_indels_5kb_stats.txt',sep = '\t',header = T)
chr29_indel <- chr29_indel[,-8]
chr29_indel$indel <- chr29_indel$插入总长度+chr29_indel$删除总长度
chr29_indel <- chr29_indel[,c(1,2,3,8)]
colnames(chr29_indel) <- c('chr','start','end','indel')
chr29_indel$end <- chr29_indel$end+1

#OR gene
ORgene <- read.table('OR_gene_loc.txt',sep = '\t',header = T)
colnames(ORgene) <- c('chr','start','end','gene')
ORgene$chr <- paste0('chr',ORgene$chr)
ORgene <- ORgene[ORgene$chr != "chr16",]
#ORgene$OR <- 1
chr29_length <- 6840000
all_positions <- data.frame(
  chr = "chr29",
  start = 0,
  end = chr29_length
)

interval_subtract <- function(df1, df2) {
  result <- data.frame()
  for (i in 1:nrow(df1)) {
    current_start <- df1$start[i]
    current_end <- df1$end[i]
    for (j in 1:nrow(df2)) {
      subtract_start <- df2$start[j]
      subtract_end <- df2$end[j]
      if (current_start < subtract_start) {
        new_row <- data.frame(
          chr = df1$chr[i],
          start = current_start,
          end = pmin(subtract_start - 1, current_end),
          gene = 0
        )
        result <- rbind(result, new_row)
        current_start <- pmax(subtract_end + 1, current_start)
      }
      if (current_start > subtract_end) {
        next
      }
      current_start = pmax(current_start, subtract_end + 1)
    }
    if (current_start <= current_end) {
      new_row <- data.frame(
        chr = df1$chr[i],
        start = current_start,
        end = current_end,
        gene = 0
      )
      result <- rbind(result, new_row)
    }
  }
  return(result)
}
non_gene_intervals <- interval_subtract(all_positions, ORgene)
chr29_ORgene <- rbind(ORgene[, c("chr", "start", "end", "gene")], non_gene_intervals)
chr29_ORgene <- chr29_ORgene[order(chr29_ORgene$start, chr29_ORgene$end), ]
chr29_ORgene$gene[chr29_ORgene$gene!= 0] <- 1

# RNAseq
library(tidyverse)
rnaseq <- read.table('../20250718_rnaseq_reads_stat/merged_counts_final.txt',sep = '\t')
rnaseq <- rnaseq[,c(1,2,3,103)]
colnames(rnaseq) <- c('chr','start','end','rnaseq')
rnaseq1 <- rnaseq %>% separate(rnaseq, into = c("new_col1", "new_col2"), sep = " ")
rnaseq1 <- rnaseq1[,c(1,2,3,5)]
colnames(rnaseq1) <- c('chr','start','end','rnaseq')    
chr29_rnaseq <- rnaseq1[rnaseq1$chr == 'chr29',]

# pi
pi_d <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis/DP4_2allele_50miss_maf001/1.10kbwindow/pi/chicken_2036samples_domestic_10kbwindow_5kbstep.windowed.pi',sep = '\t',header = T)
pi_d <- pi_d[pi_d$CHROM == "chr29",]
pi_d$breed <- 'Domestic'
pi_w <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis/DP4_2allele_50miss_maf001/1.10kbwindow/pi/chicken_2036samples_wild_10kbwindow_5kbstep.windowed.pi',sep = '\t',header = T)
pi_w <- pi_w[pi_w$CHROM == "chr29",]
pi_w$breed <- "RJF" 
pi <- rbind(pi_d,pi_w)

#AF
af_d <- read.table('../20250721_freq/domestic_chr29_freq.frq',sep = '\t')
colnames(af_d) <- c('chr','pos','allele','CHR_NUMBER','domestic_REF','domestic_ALT')
af_w <- read.table('../20250721_freq/wild_chr29_freq.frq',sep = '\t')
colnames(af_w) <- c('chr','pos','allele','CHR_NUMBER','wild_REF','wild_ALT')
af <- cbind(af_d,af_w)
colnames(af) 
af$domestic_REF <- gsub("^[A-Z]:","",af$domestic_REF)
af$domestic_ALT <- gsub("^[A-Z]:","",af$domestic_ALT)
af$wild_REF <- gsub("^[A-Z]:","",af$wild_REF)
af$wild_ALT <- gsub("^[A-Z]:","",af$wild_ALT)
af$domestic_REF<- as.numeric(af$domestic_REF)
af$wild_REF<- as.numeric(af$wild_REF)
af$domestic_ALT<- as.numeric(af$domestic_ALT)
af$wild_ALT<- as.numeric(af$wild_ALT)
af$ALT_diff <- af$domestic_ALT-af$wild_ALT
af$REF_diff <- af$domestic_REF-af$wild_REF
af2 <- af[,c(1,2,5,6,11,12,13,14)]
af2$ALT_AF_ratio <- af2$domestic_ALT/af2$wild_ALT

#FST
fst <- read.table('~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis/DP4_2allele_50miss_maf001/1.10kbwindow/fst/chicken_2036samples_10kbwindow_5kbstep.windowed.weir.fst',sep = '\t',header = T)
fst <- fst[fst$CHROM == "chr29",]
#
p1 <- ggplot(chr29_repeat, aes(start, lowercase)) + geom_area(fill="#E64b35b2") + theme_classic()+xlab('')+ylab('Repeat/5kb')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p1
p2 <- ggplot(chr29_repeat, aes(start, gc)) + geom_area(fill="#4dbbd5b2") + theme_classic()+xlab('')+ylab('GC/5kb')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p2
p3 <- ggplot(chr29_snp, aes(start, snp)) + geom_area(fill="#00a087b2") + theme_classic()+xlab('')+ylab('SNPs/5kb')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p3
p4 <- ggplot(chr29_indel, aes(start, indel)) + geom_area(fill="#3c5488b2") + theme_classic()+xlab('')+ylab('Indels/5kb')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p4

p5 <- ggplot(chr29_rnaseq, aes(start, rnaseq)) + geom_area(fill="#8491b4b2") + theme_classic()+xlab('')+ylab('RNAseq')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text = element_blank(),axis.title.x = element_blank(),axis.ticks = element_blank(),
        axis.line.x = element_blank())
p5
p6 <- ggplot(pi)+geom_line(mapping = aes(BIN_START,PI,color=breed),linewidth=0.2)+ theme_classic()+xlab('')+ylab('PI')+
  scale_color_manual(values = c("Domestic" = "grey50", "RJF" = "#dc0000b2"), guide = "none")+
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),legend.position = 'none')
p6

p7 <- ggplot(af2,aes(x=pos,y=ALT_diff))+geom_point(color="#91d1c2b2",size=0.1)+theme_classic()+xlab('')+ylab('AF diff')+
  scale_y_continuous(limits = c(-0.7,0.7))+geom_hline(aes(yintercept=0),linetype="dashed",color='black')+
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p7
p8 <- ggplot(fst, aes(BIN_START, WEIGHTED_FST)) + geom_area(fill="#3c5488b2") + theme_classic()+xlab('')+ylab('Fst')+
  #geom_hline(aes(yintercept=repeat_mean),linetype="dashed")+ 
  theme(panel.grid =element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.ticks.x = element_blank(),
        axis.line.x = element_blank())
p8

p9 <- ggplot(chr29_ORgene, aes(x = start, xend = end, y = 0.5, yend = 0.5, color = factor(gene), fill = factor(gene))) +
  geom_segment(size = 10) +scale_color_manual(values = c("0" = "white", "1" = "#f39b7fb2"), guide = "none") +
  scale_fill_manual(values = c("0" = "white", "1" = "#f39b7fb2"), guide = "none") +labs(x = "Chromosome 29", y = "Olfactory gene") + 
  scale_x_continuous(breaks = c(0,1e6,2e6,3e6,4e6,5e6,6e6),labels = c('0','1Mb','2Mb','3Mb','4Mb','5Mb','6Mb'))+
  theme_minimal() +theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),panel.grid = element_blank())
p9 
library(patchwork)
final_plot <- p1 / p2 / p3 / p4 / p5 /p6/p7/p8/p9
final_plot
ggsave(final_plot,filename = 'chr29.pdf',dpi = 300,width = 12,height = 8)
library(ggsci)
library("scales")
mypal =pal_npg("nrc", alpha =0.7)(9)
mypal
show_col(mypal)
