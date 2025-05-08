args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args || "-h" %in% args) {
 cat("Usage: Rscript plot.r --input <input.rds> --gene-list <genes> --ncol <ncol_value>\n")
 cat("\n")
 cat("Options:\n")
 cat(" --input Path to the input RDS file containing the single cell data.\n")
 cat(" --gene-list Comma-separated list of gene names to plot.\n")
 cat(" --ncol Number of columns for the Plot.\n")
 cat(" -h, --help Show this help message.\n")
 quit(status = 0)
} 

if (length(args) < 4) {
 stop("Usage: Rscript plot.r --input <input.rds> --gene-list <genes> --ncol <ncol_value>\n")
} 

input_file <- NULL
gene_list <- NULL
ncol <- NULL 

for (i in 1:length(args)) {
 if (args[i] == "--input") {
 input_file <- args[i + 1]
 } else if (args[i] == "--gene-list") {
 gene_list <- strsplit(args[i + 1], ",")[[1]] 
 } else if (args[i] == "--ncol") {
 ncol <- as.numeric(args[i + 1]) 
 }
} 

if (is.null(input_file) || is.null(gene_list) || is.null(ncol)) {
 stop("Arguments --input, --gene-list, and --ncol are required.\n")
} 

library(ggplot2)
library(Seurat) 

df <- readRDS(input_file) 
cat("Processing input:", input_file, "\n")
cat("Genes to plot:", paste(gene_list, collapse = ", "), "\n")
cat("Number of columns for plot:", ncol, "\n") 
p <- FeaturePlot(df,
 features = gene_list,
 min.cutoff = "q10",
 max.cutoff = "q90",
 cols = c("lightgrey", "red"),
 ncol=ncol) & 
NoLegend() & 
NoAxes() & 
theme( panel.border = element_rect(color = "black", linewidth = 1)) 

ggsave(p, filename = "marker_gene_featureplot.png", width = 8, height = 8) 

cat("Plot saved as", "marker_gene_featureplot.png", "\n")
