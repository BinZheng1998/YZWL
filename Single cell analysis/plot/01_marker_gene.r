# 使用 commandArgs 手动解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否有帮助参数（--help 或 -h）
if ("--help" %in% args || "-h" %in% args) {
  cat("Usage: Rscript plot.r --input <input.rds> --sample-name <sample_name>\n")
  cat("\n")
  cat("Options:\n")
  cat("  --input        Path to the input RDS file containing the single cell data.\n")
  cat("  --sample-name  The sample name used for the output plot file name.\n")
  cat("  -h, --help     Show this help message.\n")
  quit(status = 0)
}

# 检查参数个数
if (length(args) < 2) {
  stop("Usage: Rscript plot.r --input <input.rds> --sample-name <sample_name>\n")
}

# 解析命令行参数
input_file <- NULL
sample_name <- NULL

# 简单的参数解析
for (i in 1:length(args)) {
  if (args[i] == "--input") {
    input_file <- args[i + 1]
  } else if (args[i] == "--sample-name") {
    sample_name <- args[i + 1]
  }
}

# 确保两个参数都提供了
if (is.null(input_file) || is.null(sample_name)) {
  stop("Both --input and --sample-name arguments are required.\n")
}

# 加载所需的库
library(ggplot2)
library(Seurat)
library(scRNAtoolVis)
# 读取输入的RDS文件
df <- readRDS(input_file)

# 打印样本名称（如果需要，可以用于调试或日志记录）
cat("Processing sample:", sample_name, "\n")

# 基因列表
gene_list <- c(
  "gene-CYP17A1", "gene-STAR", "gene-CYP11A1", "gene-HSD3B1-gene-HAO2",
  "gene-HBAD", "STRG.62584-gene-HBZ", "gene-HBA1", "gene-PAX2", "gene-OSR1", 
  "gene-POSTN", "gene-COL5A2", "gene-DMRT1", "STRG.69136", "gene-COL3A1", 
  "gene-COL1A2", "gene-ACTA2", "gene-DCN", "gene-ACTG2", "gene-COL6A3", 
  "gene-TCF21", "gene-MYH11", "gene-ALDH1A3", "gene-CYP19A1", "gene-FSHR", 
  "STRG.77116-gene-KRT7", "gene-KRT18-STRG.77064", "gene-KRT24", "gene-CHGB", 
  "gene-DSP", "gene-EMB", "gene-MPEG1", "STRG.20429", "gene-CSF1R", "gene-SPI1", 
  "gene-RUNX1", "gene-SOX18", "gene-CDH5", "gene-APLNR", "gene-KDR", "gene-FGFR3", 
  "gene-TOX3", "gene-AMH", "gene-NR5A1", "gene-ESRRG", "gene-HEMGN", 
  "gene-PITX2-gene-LOC124417965", "gene-ESR1", "gene-EMX2", "gene-LHX9", "gene-GJA1", 
  "STRG.56044-STRG.56046-gene-MASP1", "gene-SOX10", "gene-MAP2", "gene-PMP22", 
  "gene-GAL", "2-MCQ4078410.1-gene-CHGA", "gene-MOXD1", "gene-SOX10", "gene-TAGLN3", 
  "gene-TFAP2B", "gene-SLC12A1", "gene-CA2", "gene-GATA3", "gene-MIOX", "gene-CLDN3", 
  "gene-ELF3", "gene-SULT1C3", "gene-MECOM", "gene-EZR", "gene-HNF1B", "gene-MAP3K13", 
  "gene-DAZL", "STRG.81315", "gene-CLDN1", "gene-PIWIL1", "gene-ALCAM", "gene-NEGR1", 
  "gene-SLC34A2", "gene-KIT", "gene-ESR1", "gene-EYA1", "gene-CCL4", "gene-IL8L2", 
  "gene-NPHS2", "gene-PODXL", "gene-PTPRO", "gene-RARRES1", "gene-UNC5C", 
  "gene-TEX2-gene-PECAM1", "gene-SMOC2", "gene-GATA2", "gene-MID1"
)

# 绘制DotPlot
p <- jjDotPlot(
  object = df,
  gene = gene_list,
  cluster.order = 40:0, 
  ytree = FALSE, 
  rescale = TRUE, 
  rescale.min = 0, 
  x.text.vjust = 0.5
)

# 保存图像
output_file <- paste0(sample_name, "_marker.png")
ggsave(p, filename = output_file, width = 18, height = 12)

# 输出完成消息
cat("Plot saved as", output_file, "\n")
