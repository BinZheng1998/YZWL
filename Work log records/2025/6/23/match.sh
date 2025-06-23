#!/bin/bash

# 设置序列数据文件数组
species=("pig" "dog" "sheep" "cattle" "chicken")

# 设置 BLAST 命令的基本路径和参数
blast_cmd="/home/bzheng/software/ncbi-blast-2.15.0+/bin/blastn"
evalue="1e-6"
outfmt="6"
threads="24"
max_target_seqs="1"

# 循环遍历数组中的物种文件
for ref in "${species[@]}"; do
    for query in "${species[@]}"; do
        # 跳过相同的物种比对（如：pig 对 pig）
        if [ "$ref" != "$query" ]; then
            # 生成以 ref 为参考数据库的比对
            $blast_cmd -query "../data/${query}/${query}_BW_LOCgene.fa" -db "../data/${ref}/${ref}_BW_LOCgene" \
            -evalue $evalue -outfmt $outfmt -num_threads $threads \
            -max_target_seqs $max_target_seqs \
            -out "../result/${ref}_ref_${query}_query.txt"

            # 生成以 query 为参考数据库的比对
            $blast_cmd -query "../data/${ref}/${ref}_BW_LOCgene.fa" -db "../data/${query}/${query}_BW_LOCgene" \
            -evalue $evalue -outfmt $outfmt -num_threads $threads \
            -max_target_seqs $max_target_seqs \
            -out "../result/${query}_ref_${ref}_query.txt"
        fi
    done
done
