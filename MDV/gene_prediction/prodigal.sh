#!/bin/bash

genome_dir="/home/binz/MDV/genome/fasta"
output_dir="/home/binz/MDV/result/Gene_prediction_Prodigal/result"

# 遍历基因组文件
for genome_file in "$genome_dir"/*.fna; do
    # 提取基因组文件名（不含路径和扩展名）
    genome_name=$(basename "$genome_file" .fna)
    
    # 运行 Prodigal 进行基因预测，输出到指定目录
    prodigal -i "$genome_file" -a "$output_dir/$genome_name.faa" -o "$output_dir/$genome_name.gff" -p single
     
    echo "Processed $genome_name"
done
