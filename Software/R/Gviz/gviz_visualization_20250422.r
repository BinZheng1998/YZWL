setwd("E:/08-R测试/gviz/")
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
data <- read.table("E:/09-ref/output_for_gviz.txt",sep = '\t',header = T)

grtrack <- GeneRegionTrack(data, 
                           genome = 'chicken',
                           shape='smallArrow',
                           #shape = "arrow",
                           lwd = 1.5, 
                           #col.axis = "black", 
                           col.title = "black",
                           background.title = "white",
                           col=NULL,
                           col.line = NULL,
                           fontcolor='black',
                           fill = "#000080",
                           min.height = 10, #调整箭头大小
                           chromosome = 'chr2', 
                           name = "Gene" ) #左侧注释内容，可以不加


#fst
fst <- read.table("E:/02-群体进化/02-鸡/new_res/selection_analysisi_BW/2kb_window/fst/chicken_BW_2kbwindow_2kbstep.windowed.weir.fst",sep = '\t',header = T)
fst2 <- fst[,c(1,2,3,5)]
fst_gr <- GRanges(seqnames = fst2$CHROM,
              ranges = IRanges(start = fst2$BIN_START, end = fst2$BIN_END),
              FST = fst2$WEIGHTED_FST)
fstTrack <- DataTrack(range = fst_gr, data = "FST", 
                      name = "FST", type = "hist",
                      fill.histogram="darkblue",
                      col.histogram=NA,
                      col.axis = "black", 
                      col.title = "black",
                      background.title = "white",
                      yTicksAt=c(0.02,0.04,0.06,0.08),
                      ylim=c(0,0.1),
                      transformation = function(x) { x[x < 0] <- 0; x })

#pi
library(dplyr)
pi1 <- read.table("E:/02-群体进化/02-鸡/new_res/selection_analysisi_BW/2kb_window/pi/chicken_BW_high_2kbwindow_2kbstep.windowed.pi",sep = '\t',header = T)
colnames(pi1)[5] <- "PI_High"
pi2 <- read.table("E:/02-群体进化/02-鸡/new_res/selection_analysisi_BW/2kb_window/pi/chicken_BW_low_2kbwindow_2kbstep.windowed.pi",sep = '\t',header = T)
colnames(pi2)[5] <- "PI_Low"
pi <- pi1%>% left_join(pi2,by=c("CHROM","BIN_START","BIN_END"))
pi$PI_Low[is.na(pi$PI_Low)] <-pi$PI_High[is.na(pi$PI_Low)]
pi$PI <- pi$PI_Low/pi$PI_High

pi_2 <- pi[,c(1,2,3,8)]
pi_gr <- GRanges(seqnames = pi_2$CHROM,
                  ranges = IRanges(start = pi_2$BIN_START, end = pi_2$BIN_END),
                  PI = pi_2$PI)
piTrack <- DataTrack(range = pi_gr, data = "PI",name = "PIratio", 
                     col.axis = "black", 
                     col.title = "black",
                     background.title = "white",
                     type = "hist",
                     fill.histogram="darkblue",
                     col.histogram=NA,
                     yTicksAt=c(5,10,15,20,30,40),
                     ylim = c(0,50))

tajimaD <- read.table("E:/02-群体进化/02-鸡/new_res/selection_analysisi_BW/2kb_window/tajimaD/chicken_BW_high_2kb.Tajima.D",sep = '\t',header = T)
tajimaD$BIN_END <- tajimaD$BIN_START+2000
tajimaD2 <-tajimaD[,c(1,2,5,4)]
tajimaD_gr <- GRanges(seqnames = tajimaD2$CHROM,
                 ranges = IRanges(start = tajimaD2$BIN_START, end = tajimaD2$BIN_END),
                 TajimaD = tajimaD2$TajimaD)
tajimaTrack <- DataTrack(range = tajimaD_gr, data = "TajimaD", 
                         name = "Tajima",
                         col.axis = "black", 
                         col.title = "black",
                         background.title = "white",
                         type = "hist",
                         fill.histogram="darkblue",
                         col.histogram=NA,
                         yTicksAt=c(-2,-1,0,1,2),
                         ylim=c(-2.5,2.5))

genomeAxis <- GenomeAxisTrack(lwd=1,col='black',fontcolor='black')
tracks <- list(fstTrack,piTrack,tajimaTrack,genomeAxis,grtrack)
plotTracks(tracks, from =18600000,to = 19700000, chromosome = "chr11",
           #background.panel = "white",
           transcriptAnnotation = "symbol",
           labelPos = "above", #above below
           #cex.axis = 0.5,
           innerMargin = 0,
           # minDistance = 1,
           #add53 = TRUE, add35 = TRUE,
           title.width =0.4)
