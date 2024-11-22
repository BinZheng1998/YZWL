#!/usr/bin/env Rscript

# 使用 commandArgs 手动解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 定义帮助信息函数
print_help <- function() {
  cat("
Usage: cellcycle.r --input INPUT_FILE --out OUTPUT_FILE [options]

Options:
  --input   Path to the input RDS file
  --out     Path to the output RDS file
  -h, --help  Show this help message and exit
")
}

# 定义帮助函数，检查并解析参数
parse_args <- function(args) {
  # 如果参数包含 -h 或 --help，显示帮助并退出
  if ("-h" %in% args || "--help" %in% args) {
    print_help()
    quit(status = 0)
  }

  # 初始化空列表以存储参数
  parsed_args <- list(input = NULL, out = NULL)
  
  # 遍历参数，逐个解析
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
  
  # 检查是否提供了必要的参数
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

# 4. 进行细胞周期评分
seurat_obj <- CellCycleScoring(seurat_obj, 
                               s.features = s.genes$V1, 
                               g2m.features = g2m.genes$V1, 
                               set.ident = TRUE)

# 5. 回归掉细胞周期的影响
cat("Regressing out cell cycle effects...\n")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score"))

# 6. 将结果保存为新的RDS文件
cat("Saving output file to:", opt$out, "\n")
saveRDS(seurat_obj, file = opt$out)

cat("Cell cycle analysis completed and saved to", opt$out, "\n")
