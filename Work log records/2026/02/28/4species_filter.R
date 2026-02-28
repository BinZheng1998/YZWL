setwd('~/project/01_evolution/convergent_evolution/07_result/20260227_BW_new_filter_analysis/')

species_list <- c("chicken", "pig", "cattle", "sheep")
directions <- c("pos", "neg")

for (sp in species_list) {
  message(paste(">>> Processing species:", sp))
  
  for (dir in directions) {
    # 1. 构建文件路径
    input_file  <- sprintf('%s/DCMS/%s_BW_10kwindow_5kstep_%s_DCMS_regions.txt', sp, sp, dir)
    gene_file   <- sprintf('%s/DCMS/%s_BW_10kwindow_5kstep_%s_DCMS_top_005_raw_regions_50kb_genes.txt', sp, sp, dir)
    output_file <- sprintf('%s/DCMS/%s_BW_10kwindow_5kstep_%s_DCMS_top005_filter_regions_50kb_genes.txt', sp, sp, dir)
    
    # 检查输入文件是否存在
    if (!file.exists(input_file)) {
      warning(paste("File missing, skipping:", input_file))
      next
    }
    
    # 2. 读取原始数据
    df <- read.table(input_file, sep = '\t', header = TRUE)
    
    # 3. 核心逻辑：根据方向(pos/neg)预过滤 PI
    if (dir == "pos") {
      df_filtered <- df[df$PI > 1, ]
    } else {
      df_filtered <- df[df$PI < 1, ]
    }
    
    # 检查过滤后是否有数据
    if (nrow(df_filtered) == 0) {
      message(paste("    - No data after PI filter for", sp, dir))
      next
    }
    
    # 4. 计算 FST 阈值并标记 top5
    # 注意：这里是在 PI 过滤后的基础上计算前 5% 的 FST
    threshold_val <- quantile(df_filtered$FST, probs = 0.95, na.rm = TRUE)
    
    df_filtered$fst_top5 <- threshold_val
    df_filtered$top5 <- 'NS'
    df_filtered$top5[df_filtered$FST > df_filtered$fst_top5 & 
                     df_filtered$p_values < df_filtered$threshold_005] <- 'top5'
    
    # 打印统计结果
    message(paste("    - Direction:", dir, "| FST Threshold:", round(threshold_val, 4)))
    print(table(df_filtered$top5))
    
    # 5. 提取 top5 区域并关联基因
    df1 <- df_filtered[df_filtered$top5 == 'top5', ]
    
    if (nrow(df1) > 0) {
      df_coords <- df1[, c(1:3)]
      colnames(df_coords) <- c('chromosome', 'start', 'end')
      
      # 读取对应的基因文件
      if (file.exists(gene_file)) {
        genes <- read.table(gene_file, sep = '\t', header = TRUE)
        # 合并
        df2 <- merge(df_coords, genes, by = c('chromosome', 'start', 'end'))
        # 输出
        write.table(df2, output_file, sep = '\t', row.names = FALSE, quote = FALSE)
      } else {
        warning(paste("      Gene file missing:", gene_file))
      }
    } else {
      message("    - No top5 regions found.")
    }
  }
}
