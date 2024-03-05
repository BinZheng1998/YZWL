#!/bin/bash

# 设置输入基因组文件夹路径
input_folder="/home/binz/MDV/genome/fasta"

# 设置输出文件夹路径
output_folder="/home/binz/MDV/result/pan-genome/prokka"

# 遍历输入文件夹中的所有基因组文件
for genome_file in ${input_folder}/*.fna; do
    # 提取基因组文件的文件名（不含路径和扩展名）
    genome_name=$(basename ${genome_file} .fna)

    # 执行 Prokka 进行注释
    prokka ${genome_file} \
           --outdir ${output_folder}/${genome_name}_annotation \
           --prefix ${genome_name} \
           --kingdom Viruses
done
