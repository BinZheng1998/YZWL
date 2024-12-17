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
# 读取输入的RDS文件
df <- readRDS(input_file)

# 打印样本名称（如果需要，可以用于调试或日志记录）
cat("Processing sample:", sample_name, "\n")

mycolors <- c('#7FC97F','#BEAED4','#FDC086','#FFFF99','#386CB0','#F0027F','#BF5B17',
            '#636363','#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02',
            '#A6761D','#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#E31A1C',
            '#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#B15928','#FBB4AE','#B3CDE3',
            '#CCEBC5','#DECBE4','#FED9A6','#FFFFCC','#E5D8BD','#FDDAEC','#F2F2F2',
            '#B3E2CD','#FDCDAC')

violin_data <- df@meta.data[, c("seurat_clusters", "nFeature_RNA")]
p<-ggplot(violin_data, aes(x = seurat_clusters, y = nFeature_RNA,fill=seurat_clusters)) +
  geom_violin() + scale_fill_manual(values=mycolors)+
  theme_classic() +
  labs(title="nFeature_RNA",
       x = " ",
       y = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="None",plot.title = element_text(hjust = 0.5))

output_file <- paste0(sample_name, "_cluster_gene_number.png")
ggsave(p, filename = output_file, width = 12, height = 4)

cat("Plot saved as", output_file, "\n")
