 加载必要的包
library(data.table)

# 函数：检查单个 FASTA 文件的序列长度
check_fasta_lengths <- function(file_path) {
  # 读取 FASTA 文件
  tryCatch({
    raw_seqs <- read.FASTA(file_path, type = "AA")
    lengths <- nchar(raw_seqs)
    
    # 统计长度分布
    length_table <- table(lengths)
    
    # 判断序列长度是否一致
    is_consistent <- length(unique(lengths)) == 1
    
    # 返回结果
    list(
      file = basename(file_path),
      length_distribution = length_table,
      is_consistent = is_consistent,
      num_sequences = length(lengths)
    )
  }, error = function(e) {
    # 错误处理
    message(sprintf("错误: 无法读取文件 %s: %s", file_path, e$message))
    return(NULL)
  })
}

# 函数：处理文件夹中的所有 FASTA 文件
process_fasta_files <- function(folder_path) {
  # 获取文件夹中所有 .fa 文件
  fasta_files <- list.files(path = folder_path, pattern = "trimmed\\.fa$", full.names = TRUE)
  
  if (length(fasta_files) == 0) {
    message("文件夹中没有找到 .fa 文件！")
    return(NULL)
  }
  
  # 存储结果
  results <- list()
  
  # 遍历每个 FASTA 文件
  for (file in fasta_files) {
    result <- check_fasta_lengths(file)
    if (!is.null(result)) {
      results[[length(results) + 1]] <- result
    }
  }
  
  # 创建表格
  output_table <- data.table(
    File = character(),
    Num_Sequences = integer(),
    Length_Distribution = character(),
    Is_Consistent = logical()
  )
  
  # 填充表格
  for (res in results) {
    # 将长度分布转换为字符串
    len_dist <- paste(names(res$length_distribution), 
                      res$length_distribution, 
                      sep = ":", 
                      collapse = ", ")
    output_table <- rbind(output_table, data.table(
      File = res$file,
      Num_Sequences = res$num_sequences,
      Length_Distribution = len_dist,
      Is_Consistent = res$is_consistent
    ))
  }
  
  # 输出表格
  print("序列长度分析结果：")
  print(output_table)
  
  # 保存表格到 CSV 文件
  output_file <- file.path(folder_path, "fasta_length_analysis.csv")
  fwrite(output_table, output_file)
  message(sprintf("结果已保存到: %s", output_file))
  
  return(output_table)
}

# 设置文件夹路径（替换为你实际的文件夹路径）
folder_path <- "alignment_mafft"  # 例如 "test" 文件夹

# 运行分析
result_table <- process_fasta_files(folder_path)
