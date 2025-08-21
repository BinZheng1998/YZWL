setwd(dir = '~/project/01_evolution/convergent_evolution/07_result/20250805_annovar/domestication/')
library(dplyr)
library(tidyr)
library(stringr)
library(aplot)
library(funkyheatmap)
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

#longdata1 <- result_df %>%pivot_longer(cols = c('nonCDS_SNP','Exon_SNP','Intron_SNP'),names_to = 'Region',values_to = 'SNP_count')
longdata1 <- longdata1[!longdata1$Chromosome %in% c('chrZ','chrW'),]
chromosome <- paste0("chr",1:38)

longdata1$Chromosome <- factor(longdata1$Chromosome,levels = chromosome)
longdata1 <- longdata1[order(longdata1$Chromosome), ]
longdata1 <- longdata1[,c(1,6,7,8)]
head(longdata1)
colnames(longdata1)[1] <- 'id'

cinfo <- tibble(
  id = colnames(longdata1),
  group = c(NA, NA, NA, NA),
  options = lapply(seq(4), function(x) lst()))
cinfo$name <- c('Chromosome','NonCoding','Intron','Exon')
cinfo$palette <- c(NA, "NonCoding_palette","Intron_palette","Exon_palette" )
palettes <- list(NonCoding_palette = "Blues", Exon_palette = "Greens", Intron_palette = "Reds")
cinfo$geom <- c("text", "bar", "bar", "bar")

funky_heatmap(longdata1,column_info = cinfo,palettes = palettes)
#
data1 <- read.table('chicken/chicken.avinput.exonic_variant_function',sep = '\t')
longdata2 <- data1 %>%
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
#longdata2 <-longdata2 %>%pivot_longer(cols = c('Nonsynonymous','Synonymous','Unknown'),names_to = 'Type', values_to = 'SNP_count')

longdata2 <- longdata2[!longdata2$Chromosome %in% c('chrZ','chrW'),]
longdata2$Chromosome <- factor(longdata2$Chromosome,levels = chromosome)
longdata2 <- longdata2[order(longdata2$Chromosome), ]
longdata2 <- longdata2[,c(1,6,7,8)]
longdata2 <- longdata2 %>%
  mutate(across(.cols = c(Nonsynonymous_ratio, Synonymous_ratio, Unknown_ratio),
                .fns = ~ purrr::map(., list)))

longdata2 <- longdata2 %>%
  mutate(across(
    .cols = c(Nonsynonymous_ratio, Synonymous_ratio, Unknown_ratio),
    .fns = ~ purrr::map(., ~ c("value" = .x, "others" = 1 - .x))
  ))
longdata <- cbind(longdata1,longdata2)
longdata <- longdata[,-5]

#longdata <- longdata[,-5]
#longdata$pie_data <- longdata$exon_ratio
funky_heatmap(longdata)
cinfo <- tibble(
  id = colnames(longdata),
  group = c(NA, NA, NA, NA,NA,NA,NA),
  options = lapply(seq(7), function(x) lst()))
cinfo$name <- c('','NonCoding','Intron','Exon','Nonsynonymous','Synonymous','Unknown')
cinfo$palette <- c(NA, "NonCoding_palette","Intron_palette","Exon_palette" ,
                   "Nonsynonymous_palette","Synonymous_palette","Unknown_palette")
palettes <- list(
  NA,
  NonCoding_palette = "Blues", 
  Intron_palette = "Reds",
  Exon_palette = "Greens",
  Nonsynonymous_palette = c(value = "red", others = "white"),
  Synonymous_palette = c(value = "red", others = "white"),
  Unknown_palette = c(value = "red", others = "white")
)
cinfo$geom <- c("text", "bar", "bar", "bar","pie","pie","pie")
#图例

disabled_legends = list(
  list(palette = "Nonsynonymous_palette",enabled = FALSE),
  list(palette = "Synonymous_palette",enabled = FALSE))
legends <- c(legends, disabled_legends)

cinfo[[1, "options"]] <- list(list(width = 3))
cinfo[[2, "options"]] <- list(list(width = 6,draw_outline=F))
cinfo[[3, "options"]] <- list(list(width = 6,draw_outline=F))
cinfo[[4, "options"]] <- list(list(width = 6,draw_outline=F))
cinfo[[5, "options"]] <- list(list(width = 2,draw_outline=F))
cinfo[[6, "options"]] <- list(list(width = 2,draw_outline=F))
cinfo[[7, "options"]] <- list(list(width = 2))

funky_heatmap(longdata,column_info = cinfo,palettes = palettes)
