#!/bin/bash

# 读取 1.txt 每一行
while IFS=$'\t' read -r fst pi1 pi2 species prefix outdir
do
  echo "处理物种: $species | 输出目录: $outdir"

  # 确保输出目录存在
  mkdir -p "$outdir"

  # 步骤 1: 执行 DCMS.r
  echo "运行 DCMS.r ..."
  /usr/bin/Rscript /home/bzheng/project/01_evolution/script/DCMS_202506.r \
    --input-fst "$fst" \
    --input-pi1 "$pi1" \
    --input-pi2 "$pi2" \
    --species "$species" \
    --out-prefix "$prefix" \
    --output-folder "$outdir"

  # 步骤 2: 对 merged_regions.txt 文件执行 gene.sh
  echo "查找并处理 merged_regions.txt 文件 ..."
  for merged_file in "$outdir"/*merged_regions.txt; do
    [ -e "$merged_file" ] || continue  # 跳过不存在的情况
    echo "处理 $merged_file ..."
    bash /home/bzheng/project/01_evolution/script/get_genes_202506.sh \
      --input "$merged_file" \
      --species "$species" \
      --output "$outdir"
  done

  # 步骤 3: 对 genes.txt 文件执行 GO 富集分析
  echo "运行 GO 富集分析 ..."
  /usr/bin/Rscript /home/bzheng/project/01_evolution/script/GO_20250607.r* \
    --input-path "$outdir" \
    --species "$species" \
    --output-folder "$outdir"

  echo "完成: $prefix"
  echo "-----------------------------"

done < DCMS2GO_input.txt
