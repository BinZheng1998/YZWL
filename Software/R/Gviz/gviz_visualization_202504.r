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
                           col="#000080",
                           fontcolor='black',
                           fill = "#000080",
                           min.height = 10, #调整箭头大小
                           chromosome = 'chr2', 
                           name = "Gene" ) #左侧注释内容，可以不加

#displayPars(grtrack) <- list(min.height = 10)
plotTracks(list(grtrack),from = 20000,to = 120000,
           transcriptAnnotation = "symbol",
           #transcriptAnnotation = "transcript",
           background.panel = "white", 
           #background.title = 'darkblue',
           #collapseTranscripts = "longest",
           title.width =0, #去掉左侧注释框
           col.line = NULL, col = NULL)

data <- read.table('E:/02-群体进化/02-鸡/new_res/selection_analysisi_BW/DCMS-fst-pi-tajimaD/regions/chicken_BW_neg_DCMS_regions.txt',sep = '\t',header = T)
df <- data[,c(1,2,3,4,5,6)]

gr <- GRanges(seqnames = df$CHROM,
              ranges = IRanges(start = df$BIN_START, end = df$BIN_END),
              FST = df$FST,
              PI = df$PI,
              TajimaD = df$tajimad_l)
dataTrack <- DataTrack(range = gr, name = "chicken BW")
fstTrack <- DataTrack(range = gr, data = "FST", 
                      name = "FST", type = "l",
                      col.axis = "black", 
                      col.title = "black",
                      background.title = "white",
                      yTicksAt=c(0.02,0.04,0.06,0.08),
                      ylim=c(0,0.1),
                      transformation = function(x) { x[x < 0] <- 0; x })
piTrack <- DataTrack(range = gr, data = "PI", 
                     col.axis = "black", 
                     col.title = "black",
                     background.title = "white",
                     name = "PI", type = "l",
                     yTicksAt=c(2,4,6,8),
                     ylim = c(0,10))
tajimaTrack <- DataTrack(range = gr, data = "TajimaD", 
                         col.axis = "black", 
                         col.title = "black",
                         background.title = "white",
                         name = "TajimaD",type = "hist",
                         fill.histogram="darkblue",
                         yTicksAt=c(-2,-1,0,1,2),
                         ylim=c(-2.5,2.5))

genomeAxis <- GenomeAxisTrack(lwd=1,col='black',fontcolor='black',name='chr2')
tracks <- list(fstTrack, piTrack, tajimaTrack, genomeAxis,grtrack)
plotTracks(tracks, from = 210000,to = 300000, chromosome = "chr2",
           background.panel = "white",
           transcriptAnnotation = "symbol",
           labelPos = "below", #above
           #cex.axis = 0.5,
           innerMargin = 0,
           #add53 = TRUE, add35 = TRUE,
           title.width =0.4)
