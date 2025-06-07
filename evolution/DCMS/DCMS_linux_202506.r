#!/usr/bin/env Rscript

# 加载必要的 R 包
library(tidyverse)
library(MINOTAUR)
library(rrcovNA)
library(MASS)
library(qvalue)
library(optparse)

# 设置选项以避免科学计数法
options(scipen = 999)

# 定义命令行选项 (放在脚本最前面)
option_list <- list(
  make_option(c("--input-fst"), type = "character", default = NULL, 
              help = "Path to the input FST file (required)", metavar = "character"),
  make_option(c("--input-pi1"), type = "character", default = NULL, 
              help = "Path to the input pi1 file for high group (required)", metavar = "character"),
  make_option(c("--input-pi2"), type = "character", default = NULL, 
              help = "Path to the input pi2 file for low group (required)", metavar = "character"),
  make_option(c("--species"), type = "character", default = "chicken", 
              help = "Species to analyze: chicken, pig, cattle, sheep, dog [default: %default]", 
              metavar = "character"),
  make_option(c("--out-prefix"), type = "character", default = "output", 
              help = "Prefix for output files [default: %default]", metavar = "character"),
  make_option(c("--output-folder"), type = "character", default = ".", 
              help = "Folder where output files will be saved [default: %default]", 
              metavar = "character")
)

# 创建 OptionParser 对象，添加 usage 和脚本描述
opt_parser <- OptionParser(option_list = option_list, 
                           usage = "usage: %prog [options]\nThis script performs positive and negative selection analysis using the DCMS method.")

# 解析命令行参数
opt <- parse_args(opt_parser, args = commandArgs(trailingOnly = TRUE))

# 调试：打印解析的参数
cat("Parsed arguments:\n")
print(opt)

# 检查必需的参数是否提供
if (is.null(opt$`input-fst`) || is.null(opt$`input-pi1`) || is.null(opt$`input-pi2`)) {
  print_help(opt_parser)
  stop("Error: Required arguments --input-fst, --input-pi1, and --input-pi2 must be provided.", call. = FALSE)
}

# 扩展 ~ 为完整路径
opt$`input-fst` <- path.expand(opt$`input-fst`)
opt$`input-pi1` <- path.expand(opt$`input-pi1`)
opt$`input-pi2` <- path.expand(opt$`input-pi2`)
opt$`output-folder` <- path.expand(opt$`output-folder`)

# 调试：打印扩展后的路径
cat("Expanded file paths:\n")
cat("input-fst:", opt$`input-fst`, "\n")
cat("input-pi1:", opt$`input-pi1`, "\n")
cat("input-pi2:", opt$`input-pi2`, "\n")
cat("output-folder:", opt$`output-folder`, "\n")

# 检查输入文件是否存在
if (!file.exists(opt$`input-fst`)) stop(paste("Input FST file does not exist:", opt$`input-fst`))
if (!file.exists(opt$`input-pi1`)) stop(paste("Input pi1 file does not exist:", opt$`input-pi1`))
if (!file.exists(opt$`input-pi2`)) stop(paste("Input pi2 file does not exist:", opt$`input-pi2`))

# 检查输出文件夹是否存在，如果不存在则创建
if (!dir.exists(opt$`output-folder`)) {
  dir.create(opt$`output-folder`, recursive = TRUE)
  cat("Created output folder:", opt$`output-folder`, "\n")
}

# 定义自定义函数

add_chr_prefix <- function(chrom) {
  ifelse(grepl("^chr", chrom, ignore.case = TRUE), chrom, paste0("chr", chrom))
}

## 数据检查函数
data_checks2 <- function(dfv, column.nums, subset, S, M, check.na = TRUE, check.S = TRUE, check.M = FALSE) {
  if (!is.matrix(dfv) & !is.data.frame(dfv)) stop("dfv must be a matrix or data frame")
  if (inherits(try(dfv[, column.nums], silent = TRUE), "try-error")) stop("column.nums must contain valid indexes for choosing columns in dfv")
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  if (any(!apply(df.vars, 2, is.numeric))) stop("all selected columns of dfv must be numeric")
  if (check.na & any(is.na(df.vars))) stop("dfv cannot contain NA values")
  if (nrow(df.vars) < 2) stop("dfv must contain at least two rows")
  if (inherits(try(df.vars[subset, ], silent = TRUE), "try-error")) stop("subset must contain valid indexes for choosing rows in dfv")
  df.vars_subset <- as.matrix(df.vars[subset, , drop = FALSE])
  if (nrow(df.vars_subset) < 2) stop("subset must index at least two rows in dfv")
  if (check.S) {
    if (is.null(S)) S <- stats::cov(df.vars_subset, use = "pairwise.complete.obs")
    if (!is.matrix(S)) stop("S must be a matrix")
    if (nrow(S) != ncol(df.vars) | ncol(S) != ncol(df.vars)) stop("S must contain the same number of rows and columns as there are selected variables in dfv")
    if (any(is.na(S))) stop("covariance matrix S contains NA values")
    if (inherits(try(solve(S), silent = TRUE), "try-error")) stop("covariance matrix S is exactly singular")
    S_inv <- solve(S)
  } else {
    S <- NULL
    S_inv <- NULL
  }
  if (check.M) {
    if (is.null(M)) M <- colMeans(df.vars_subset, na.rm = TRUE)
    M <- as.vector(unlist(M))
    if (length(M) != ncol(df.vars)) stop("M must contain one value per selected column of dfv")
  } else {
    M <- NULL
  }
  return(list(S = S, S_inv = S_inv, M = M))
}

## DCMS 计算函数
DCMS2 <- function(dfv, column.nums = 1:ncol(dfv), subset = 1:nrow(dfv), S = NULL, dfp, column.nums.p = 1:ncol(dfp)) {
  if (length(column.nums) != length(column.nums.p)) stop("column.nums must contain same number of values as column.nums.p")
  dfv_check <- data_checks2(dfv, column.nums, subset, S, M = NULL, check.na = TRUE, check.M = FALSE)
  dfp_check <- data_checks2(dfp, column.nums.p, subset, S = NULL, M = NULL, check.na = TRUE, check.S = FALSE, check.M = FALSE)
  if (nrow(dfv) != nrow(dfp)) stop("dfv and dfp must contain the same number of entries")
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  n <- nrow(df.vars)
  d <- ncol(df.vars)
  df.p <- as.matrix(dfp[, column.nums.p, drop = FALSE])
  S <- dfv_check$S
  corrMat <- S / sqrt(outer(diag(S), diag(S)))
  DCMS <- 0
  for (i in 1:d) {
    DCMS <- DCMS + (log(1 - df.p[, i]) - log(df.p[, i])) / sum(abs(corrMat[i, ]))
  }
  return(DCMS)
}

## 统计量转 p 值函数
stat_to_pvalue2 <- function(dfv, column.nums = 1:ncol(dfv), subset = 1:nrow(dfv), 
                            two.tailed = rep(TRUE, length(column.nums)), right.tailed = rep(FALSE, length(column.nums))) {
  dfv_check <- data_checks2(dfv, column.nums, subset, S = NULL, M = NULL, check.na = TRUE, check.S = FALSE, check.M = FALSE)
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  n <- nrow(df.vars)
  d <- ncol(df.vars)
  df.p <- as.data.frame(matrix(0, n, d))
  if (length(two.tailed) != d) stop("two.tailed must be a vector of same length as column.nums")
  if (length(right.tailed) != d) stop("right.tailed must be a vector of same length as column.nums")
  noSubset <- (length(subset) == nrow(dfv))
  if (noSubset) {
    for (i in 1:d) {
      df.p[, i] <- (rank(df.vars[, i]) - 1) / (n - 1)
      if (two.tailed[i]) {
        df.p[, i] <- 1 - 2 * abs(df.p[, i] - 0.5)
      } else if (right.tailed[i]) {
        df.p[, i] <- 1 - df.p[, i]
      }
      df.p[, i] <- (df.p[, i] * n + 1) / (n + 2)
    }
  } else {
    df.vars_subset <- as.matrix(df.vars[subset, , drop = FALSE])
    n2 <- nrow(df.vars_subset)
    for (i in 1:d) {
      df.p[, i] <- findInterval(df.vars[, i], sort(df.vars_subset[, i])) / n2
      if (two.tailed[i]) {
        df.p[, i] <- 1 - 2 * abs(df.p[, i] - 0.5)
      } else if (right.tailed[i]) {
        df.p[, i] <- 1 - df.p[, i]
      }
      df.p[, i] <- (df.p[, i] * n2 + 1) / (n2 + 2)
    }
  }
  return(df.p)
}

## 合并区间函数
merge_intervals <- function(df, threshold) {
  colnames(df) <- c("chromosome", "start", "end")
  df <- df[order(df$chromosome, df$start), ]
  merged <- data.frame(chromosome = character(), start = integer(), end = integer(), stringsAsFactors = FALSE)
  current_row <- df[1, ]
  for (i in 2:nrow(df)) {
    if (df$chromosome[i] == current_row$chromosome && df$start[i] <= current_row$end + threshold) {
      current_row$end <- max(current_row$end, df$end[i])
    } else {
      merged <- rbind(merged, current_row)
      current_row <- df[i, ]
    }
  }
  merged <- rbind(merged, current_row)
  return(merged)
}

# 根据物种定义要排除的染色体
exclude_chrom <- switch(opt$species,
                        "chicken" = c("Z", "W", "MT","chrZ","chrW","chrMT"),
                        "pig" = c("19", "20", "22","chr19","chr20","chr22"),
                        "sheep" = c("27", "28","chr27","chr28"),
                        "cattle" = c("X", "Y","chrX","chrY"),
                        "dog" = c("chrX","X"),
                        stop("Invalid species. Choose from chicken, pig, cattle, sheep, dog."))

# 读取输入文件
fst <- read.table(opt$`input-fst`, sep = "\t", header = TRUE)
fst$BIN_END <- fst$BIN_END + 1
pi_1 <- read.table(opt$`input-pi1`, sep = "\t", header = TRUE)
pi_1$BIN_END <- pi_1$BIN_END + 1 
pi_2 <- read.table(opt$`input-pi2`, sep = "\t", header = TRUE)
pi_2$BIN_END <- pi_2$BIN_END + 1

# 数据预处理
fst <- fst[, c(1, 2, 3, 5)]  # 选择 CHROM, BIN_START, BIN_END, FST 列
fst <- fst[!fst$CHROM %in% exclude_chrom, ]  # 排除指定染色体
fst$CHROM <- add_chr_prefix(fst$CHROM)
#fst$CHROM <- paste0("chr", fst$CHROM)  # 添加 "chr" 前缀

pi_1 <- pi_1[, c(1, 2, 3, 5)]  # 选择 CHROM, BIN_START, BIN_END, PI 列
pi_2 <- pi_2[, c(1, 2, 3, 5)]
pi <- pi_1 %>%
  dplyr::left_join(pi_2, by = c("CHROM", "BIN_START", "BIN_END")) %>%
  dplyr::rename(pi_1 = PI.x, pi_2 = PI.y) %>%
  dplyr::mutate(pi2VSpi1 = pi_2 / pi_1)
pi <- pi[!pi$CHROM %in% exclude_chrom, ]
#pi$CHROM <- paste0("chr", pi$CHROM)
pi$CHROM <- add_chr_prefix(pi$CHROM)
pi <- pi[, c(1, 2, 3, 6)]  # 保留 CHROM, BIN_START, BIN_END, pi2VSpi1

df <- fst %>%
  dplyr::left_join(pi, by = c("CHROM", "BIN_START", "BIN_END"))
colnames(df) <- c("CHROM", "BIN_START", "BIN_END", "FST", "PI")
df$FST[is.na(df$FST)] <- 0
df$PI[is.na(df$PI)] <- 1
df2 <- as.data.frame(df[, c(4, 5)])  # 提取 FST 和 PI (pi2VSpi1)

# 正选择分析
p <- stat_to_pvalue2(dfv = df2, column.nums = 1:2, two.tailed = c(FALSE, FALSE), right.tailed = c(TRUE, TRUE))
colnames(p) <- c("fst", "pi")
mcd <- CovNAMcd(p, alpha = 0.75, nsamp = 30000)
dfv <- as.matrix(df2)
dfp <- as.matrix(p)
dcms <- DCMS2(dfv = dfv, S = mcd@cov, dfp = dfp)
rlm_model <- rlm(dcms ~ 1)
mu <- coef(rlm_model)[1]
sigma <- summary(rlm_model)$sigma
p_values <- pnorm(dcms, mean = mu, sd = sigma, lower.tail = FALSE)
qobj <- qvalue(p_values)
q_values <- qobj$qvalues
q_values_BH <- p.adjust(p = p_values, method = "BH")
final <- cbind(df[, c(1, 2, 3)], dfv, dcms, p_values, q_values, q_values_BH)

# 选择 top threshold% 的区域
p_001 <- quantile(p_values, 0.01)
data1 <- final[p_values <= p_001, c(1, 2, 3)]
p_005 <- quantile(p_values, 0.05)
data2 <- final[p_values <= p_005, c(1, 2, 3)]

# 写入正选择结果
out_prefix <- file.path(opt$`output-folder`, opt$`out-prefix`)
write.table(final, file = paste0(out_prefix, "_pos_DCMS_regions.txt"), row.names = FALSE, sep = "\t",quote = FALSE)
write.table(data1, file = paste0(out_prefix, "_pos_DCMS_top_001_raw_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
write.table(data2, file = paste0(out_prefix, "_pos_DCMS_top_005_raw_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
result1 <- merge_intervals(data1, 20000)
write.table(result1, file = paste0(out_prefix, "_pos_DCMS_top_001_merged_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
result2 <- merge_intervals(data2, 20000)
write.table(result2, file = paste0(out_prefix, "_pos_DCMS_top_005_merged_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)


# 负选择分析
p <- stat_to_pvalue2(dfv = df2, column.nums = 1:2, two.tailed = c(FALSE, FALSE), right.tailed = c(TRUE, FALSE))
colnames(p) <- c("fst", "pi")
mcd <- CovNAMcd(p, alpha = 0.75, nsamp = 300000)
dfv <- as.matrix(df2)
dfp <- as.matrix(p)
dcms <- DCMS2(dfv = dfv, S = mcd@cov, dfp = dfp)
rlm_model <- rlm(dcms ~ 1)
mu <- coef(rlm_model)[1]
sigma <- summary(rlm_model)$sigma
p_values <- pnorm(dcms, mean = mu, sd = sigma, lower.tail = FALSE)
qobj <- qvalue(p_values)
q_values <- qobj$qvalues
q_values_BH <- p.adjust(p = p_values, method = "BH")
final2 <- cbind(df[, c(1, 2, 3)], dfv, dcms, p_values, q_values, q_values_BH)

# 选择 top threshold% 的区域
p_001 <- quantile(p_values, 0.01)
data3 <- final2[p_values <= p_001, c(1, 2, 3)]
p_005 <- quantile(p_values, 0.05)
data4 <- final2[p_values <= p_005, c(1, 2, 3)]

# 写入负选择结果
write.table(final2, file = paste0(out_prefix, "_neg_DCMS_regions.txt"), row.names = FALSE, sep = "\t",quote = FALSE)
write.table(data3, file = paste0(out_prefix, "_neg_DCMS_top_001_raw_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
result3 <- merge_intervals(data3, 20000)
write.table(result3, file = paste0(out_prefix, "_neg_DCMS_top_001_merged_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
result4 <- merge_intervals(data4, 20000)
write.table(result4, file = paste0(out_prefix, "_neg_DCMS_top_005_merged_regions.txt"), row.names = FALSE, col.names = FALSE, sep = "\t",quote = FALSE)
