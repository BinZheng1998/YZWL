library(tidyverse)
library(MINOTAUR)
options(scipen = 999) 
data_checks2 <- function (dfv, column.nums, subset, S, M, check.na = TRUE, check.S = TRUE, 
                          check.M = FALSE) 
{
  if (!is.matrix(dfv) & !is.data.frame(dfv)) 
    stop("dfv must be a matrix or data frame")
  if (inherits(try(dfv[, column.nums], silent = TRUE), "try-error")) 
    stop("column.nums must contain valid indexes for choosing columns in dfv")
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  if (any(!(apply(df.vars, 2, is.numeric)))) 
    stop("all selected columns of dfv must be numeric")
  if (check.na & any(is.na(df.vars))) 
    stop("dfv cannot contain NA values")
  if (nrow(df.vars) < 2) 
    stop("dfv must contain at least two rows")
  if (inherits(try(df.vars[subset, ], silent = TRUE), "try-error")) 
    stop("subset must contain valid indexes for choosing rows in dfv")
  df.vars_subset <- as.matrix(df.vars[subset, , drop = FALSE])
  if (nrow(df.vars_subset) < 2) 
    stop("subset must index at least two rows in dfv")
  if (check.S) {
    if (is.null(S)) 
      S <- stats::cov(df.vars_subset, use = "pairwise.complete.obs")
    if (!is.matrix(S)) 
      stop("S must be a matrix")
    if (nrow(S) != ncol(df.vars) | ncol(S) != ncol(df.vars)) 
      stop("S must contain the same number of rows and columns as there are selected variables in dfv")
    if (any(is.na(S))) 
      stop("covariance matrix S contains NA values")
    if (inherits(try(solve(S), silent = TRUE), "try-error")) 
      stop("covariance matrix S is exactly singular")
    S_inv <- solve(S)
  }
  else {
    S <- NULL
    S_inv <- NULL
  }
  if (check.M) {
    if (is.null(M)) 
      M <- colMeans(df.vars_subset, na.rm = TRUE)
    M <- as.vector(unlist(M))
    if (length(M) != ncol(df.vars)) 
      stop("M must contain one value per selected column of dfv")
  }
  else {
    M <- NULL
  }
  output <- list(S = S, S_inv = S_inv, M = M)
  return(output)
}

DCMS2 <-function (dfv, column.nums = 1:ncol(dfv), subset = 1:nrow(dfv), 
                  S = NULL, dfp, column.nums.p = 1:ncol(dfp)) 
{
  if (length(column.nums) != length(column.nums.p)) 
    stop("column.nums must contain same number of values as column.nums.p")
  dfv_check <- data_checks2(dfv, column.nums, subset, S, M = NULL, 
                            check.na = TRUE, check.M = FALSE)
  dfp_check <- data_checks2(dfp, column.nums.p, subset, S = NULL, 
                            M = NULL, check.na = TRUE, check.S = FALSE, check.M = FALSE)
  if (nrow(dfv) != nrow(dfp)) 
    stop("dfv and dfp must contain the same number of entries")
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  n <- nrow(df.vars)
  d <- ncol(df.vars)
  df.p <- as.matrix(dfp[, column.nums.p, drop = FALSE])
  S <- dfv_check$S
  corrMat <- S/sqrt(outer(diag(S), diag(S)))
  DCMS <- 0
  for (i in 1:d) {
    DCMS <- DCMS + (log(1 - df.p[, i]) - log(df.p[, i]))/sum(abs(corrMat[i, 
    ]))
  }
  return(DCMS)
}

stat_to_pvalue2 <- function (dfv, column.nums = 1:ncol(dfv), subset = 1:nrow(dfv), 
                             two.tailed = rep(TRUE, length(column.nums)), right.tailed = rep(FALSE, 
                                                                                             length(column.nums))) 
{
  dfv_check <- data_checks2(dfv, column.nums, subset, S = NULL, 
                            M = NULL, check.na = TRUE, check.S = FALSE, check.M = FALSE)
  df.vars <- as.matrix(dfv[, column.nums, drop = FALSE])
  n <- nrow(df.vars)
  d <- ncol(df.vars)
  df.p <- as.data.frame(matrix(0, n, d))
  if (length(two.tailed) != d) 
    stop("two.tailed must be a vector of same length as column.nums")
  if (length(right.tailed) != d) 
    stop("right.tailed must be a vector of same length as column.nums")
  noSubset <- (length(subset) == nrow(dfv))
  if (noSubset) {
    for (i in 1:d) {
      df.p[, i] <- (rank(df.vars[, i]) - 1)/(n - 1)
      if (two.tailed[i]) {
        df.p[, i] <- 1 - 2 * abs(df.p[, i] - 0.5)
      }
      else {
        if (right.tailed[i]) 
          df.p[, i] <- 1 - df.p[, i]
      }
      df.p[, i] <- (df.p[, i] * n + 1)/(n + 2)
    }
  }
  if (!noSubset) {
    df.vars_subset <- as.matrix(df.vars[subset, , drop = FALSE])
    n2 <- nrow(df.vars_subset)
    for (i in 1:d) {
      df.p[, i] <- findInterval(df.vars[, i], sort(df.vars_subset[, 
                                                                  i]))/n2
      if (two.tailed[i]) {
        df.p[, i] <- 1 - 2 * abs(df.p[, i] - 0.5)
      }
      else {
        if (right.tailed[i]) 
          df.p[, i] <- 1 - df.p[, i]
      }
      df.p[, i] <- (df.p[, i] * n2 + 1)/(n2 + 2)
    }
  }
  return(df.p)
}

#合并附近的区间 比如10kb
merge_intervals <- function(df, threshold) {
  colnames(df) <- c("chromosome", "start", "end")
  
  df <- df[order(df$chromosome, df$start), ]
  merged <- data.frame(chromosome = character(), start = integer(), end = integer())
  
  current_row <- df[1, ]
  
  for (i in 2:nrow(df)) {
    if (df$chromosome[i] == current_row$chromosome &&
        df$start[i] <= current_row$end + threshold) {
      current_row$end <- max(current_row$end, df$end[i])
    } else {
      merged <- rbind(merged, current_row)
      current_row <- df[i, ]
    }
  }
  merged <- rbind(merged, current_row)
  
  return(merged)
}

fst <- read.table('E:/02-群体进化/03-猪/selection_body_weight/10kb_window_5kb_step/top5/fst/pig_bw_10kbwindow_10kbstep.windowed.weir.fst',sep = '\t',header = T)
fst <- fst[,c(1,2,3,5)]
fst <- fst[fst$CHROM != "19",]
fst <- fst[fst$CHROM != "20",]
fst <- fst[fst$CHROM != "22",]
colnames(fst)
fst$CHROM <- paste0("chr","",fst$CHROM)

# BW high
pi_1 <- read.table('E:/02-群体进化/03-猪/selection_body_weight/10kb_window_5kb_step/top5/pi/pig_bw_high_10kbwindow_10kbstep.windowed.pi',sep = '\t',header = T)
# wild
pi_2 <- read.table('E:/02-群体进化/03-猪/selection_body_weight/10kb_window_5kb_step/top5/pi/pig_bw_low_10kbwindow_10kbstep.windowed.pi',sep = '\t',header = T)
pi <- pi_1%>% 
  left_join(pi_2,by=c("CHROM","BIN_START","BIN_END")) %>% 
  dplyr::rename(pi_1=PI.x, pi_2=PI.y)%>%
  mutate(pi2VSpi1=pi_2/pi_1)
pi <- pi[pi$CHROM != "19",]
pi <- pi[pi$CHROM != "20",]
pi <- pi[pi$CHROM != "22",]
pi$CHROM <- paste0("chr","",pi$CHROM)
pi <- pi[,c(1,2,3,8)]

df<- fst%>% 
  left_join(pi,by=c("CHROM","BIN_START","BIN_END")) 
colnames(df)<-c("CHROM","BIN_START","BIN_END","FST","PI")

df$FST[is.na(df$FST)] <- 0
df$PI[is.na(df$PI)] <- 1
df2 <- as.data.frame(df[,c(4,5)])

#pos selection

#step1
p <- stat_to_pvalue2(
  dfv = df2,                   
  column.nums = 1:2,           
  two.tailed=c(FALSE,FALSE),    # 单尾检验
  right.tailed = c(TRUE,TRUE) # 右尾检验
)

head(p)

colnames(p) <- c("fst","pi")

#step2
library(rrcovNA )
mcd <- CovNAMcd(p,alpha = 0.75,nsamp = 30000)
mcd@cov
dfv <- as.matrix(df2)
dfp <- as.matrix(p)
dcms <- DCMS2(dfv = dfv,S=mcd@cov,dfp = dfp)

#step3
library(MASS)
rlm_model <- rlm(dcms~1)
mu <- coef(rlm_model)[1]
sigma <- summary(rlm_model)$sigma  
p_values <- pnorm(dcms, mean = mu, sd = sigma, lower.tail = FALSE)

#step4
library(qvalue)
qobj <- qvalue(p_values)
q_values <- qobj$qvalues
q_values_BH <- p.adjust(p = p_values, method = "BH")
final <- cbind(df[,c(1,2,3)],dfv,dcms,p_values,q_values,q_values_BH)
write.table(final,file = "pig_BW_pos_DCMS_regions.txt",row.names = F,sep = '\t')

data1 <- final[final$p_values < 0.01,]
data1 <- data1[,c(1,2,3)]
write.table(data1,file = "pig_BW_pos_DCMS_p001_raw_regions.txt",row.names = F,col.names = F,sep = '\t')

result1 <- merge_intervals(data1, 20000)
write.table(result1,file = "pig_BW_pos_DCMS_p001_merged_regions.txt",row.names = F,col.names = F,sep = '\t')

#neg selection
#step1
p <- stat_to_pvalue2(
  dfv = df2,                   
  column.nums = 1:2,           
  two.tailed=c(FALSE,FALSE),    # 单尾检验
  right.tailed = c(TRUE,FALSE) # 右尾检验
)
head(p)
colnames(p) <- c("fst","pi")

#step2
library(rrcovNA )
mcd <- CovNAMcd(p,alpha = 0.75,nsamp = 30000)
mcd@cov
dfv <- as.matrix(df2)
dfp <- as.matrix(p)
dcms <- DCMS2(dfv = dfv,S=mcd@cov,dfp = dfp)

#step3
library(MASS)
rlm_model <- rlm(dcms~1)
mu <- coef(rlm_model)[1]
sigma <- summary(rlm_model)$sigma 
p_values <- pnorm(dcms, mean = mu, sd = sigma, lower.tail = FALSE)

#step4
library(qvalue)
qobj <- qvalue(p_values)
q_values <- qobj$qvalues
q_values_BH <- p.adjust(p = p_values, method = "BH")
final2 <- cbind(df[,c(1,2,3)],dfv,dcms,p_values,q_values,q_values_BH)
write.table(final2,file = "pig_BW_neg_DCMS_regions.txt",row.names = F,sep = '\t')

data2 <- final2[final2$p_values<0.01,]
data2 <- data2[,c(1,2,3)]
write.table(data2,file = "pig_BW_neg_DCMS_p001_raw_regions.txt",row.names = F,col.names = F,sep = '\t')

result2 <- merge_intervals(data2, 20000)
write.table(result2,file = "pig_BW_neg_DCMS_p001_merged_regions.txt",row.names = F,col.names = F,sep = '\t')
