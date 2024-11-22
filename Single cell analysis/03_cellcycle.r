#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

print_help <- function() {
  cat("
Usage: cellcycle.r --input INPUT_FILE --out OUTPUT_FILE [options]

Options:
  --input   Path to the input RDS file
  --out     Path to the output RDS file
  -h, --help  Show this help message and exit
")
}

parse_args <- function(args) {
  if ("-h" %in% args || "--help" %in% args) {
    print_help()
    quit(status = 0)
  }

  parsed_args <- list(input = NULL, out = NULL)
  
  for (i in seq(1, length(args), by = 2)) {
    arg_name <- args[i]
    arg_value <- args[i + 1]
    
    if (arg_name == "--input") {
      parsed_args$input <- arg_value
    } else if (arg_name == "--out") {
      parsed_args$out <- arg_value
    } else {
      stop(paste("Unknown argument", arg_name))
    }
  }
  
  if (is.null(parsed_args$input)) {
    stop("Please provide an input RDS file using --input option.")
  }
  
  if (is.null(parsed_args$out)) {
    stop("Please provide an output file path using --out option.")
  }
  
  return(parsed_args)
}


opt <- parse_args(args)

library(Seurat)
cat("Reading input file:", opt$input, "\n")
seurat_obj <- readRDS(opt$input)

s.genes <- read.table("/home/bzheng/project/single_cell/new_result/00-script/chicken_s.gene.txt",sep='\t')
g2m.genes <- read.table("/home/bzheng/project/single_cell/new_result/00-script/chicken_g2m.gene.txt",sep='\t')

seurat_obj <- CellCycleScoring(seurat_obj, 
                               s.features = s.genes$V1, 
                               g2m.features = g2m.genes$V1, 
                               set.ident = TRUE)

cat("Regressing out cell cycle effects...\n")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"))

cat("Saving output file to:", opt$out, "\n")
saveRDS(seurat_obj, file = opt$out)

cat("Cell cycle analysis completed and saved to", opt$out, "\n")
