#!/bin/bash

# 设置物种列表（查询物种）
ref=("pig" "dog" "sheep" "cattle" "chicken")

# 固定参考物种为human
species="human"

# BLAST参数配置
blast_cmd="/home/bzheng/software/ncbi-blast-2.15.0+/bin/blastp"
evalue="1e-6"
outfmt="6"              # 输出为制表符分隔格式
threads="48"            # 并行线程数
max_target_seqs="1"     # 每个查询序列仅保留最佳匹配

# 遍历所有查询物种
for ref in "${ref[@]}"; do
    # 执行BLAST比对：human为参考，当前物种为查询
    $blast_cmd -query "../${species}/${species}.fa" \
               -db "../${ref}/${ref}" \
               -evalue $evalue \
               -outfmt $outfmt \
               -num_threads $threads \
               -max_target_seqs $max_target_seqs \
               -out "${ref}_ref_human_query.txt"
done
