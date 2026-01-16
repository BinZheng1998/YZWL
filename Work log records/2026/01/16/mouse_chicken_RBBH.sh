#!/bin/bash

# ================= 配置区域 =================

# 定义两个物种的名称 (请确保文件名与此处一致，例如 chicken.fa 和 mouse.fa)
sp1="chicken"
sp2="mouse"

# 设置 BLAST 命令路径和参数
blast_cmd="/home/bzheng/software/ncbi-blast-2.15.0+/bin/blastp"
evalue="1e-4"
outfmt="6"
threads="48"
# max_target_seqs=1 是找 RBH 的关键，只保留最佳的一个
max_target_seqs="1"

# 数据文件目录
data_dir="/home/bzheng/project/10_RER/data"
# 结果输出目录
result_dir="/home/bzheng/project/10_RER/result/260116_mouse_chicken_RBBH"

# 创建结果目录
mkdir -p "$result_dir"

# ================= BLAST 比对阶段 =================

echo "开始 BLAST 比对..."

# 1. Chicken 查询 -> Mouse 数据库 (sp1 -> sp2)
# 输出文件: chicken_vs_mouse.txt
if [ ! -f "${result_dir}/${sp1}_vs_${sp2}.txt" ]; then
    echo "Running BLAST: ${sp1} query against ${sp2} db..."
    $blast_cmd -query "${data_dir}/${sp1}.fa" -db "${data_dir}/${sp2}_db" \
        -evalue $evalue -outfmt $outfmt -num_threads $threads -max_target_seqs $max_target_seqs \
        -out "${result_dir}/${sp1}_vs_${sp2}.txt" || {
        echo "Error running BLAST for ${sp1} -> ${sp2}"
        exit 1
    }
else
    echo "BLAST result for ${sp1} -> ${sp2} already exists, skipping..."
fi

# 2. Mouse 查询 -> Chicken 数据库 (sp2 -> sp1)
# 输出文件: mouse_vs_chicken.txt
if [ ! -f "${result_dir}/${sp2}_vs_${sp1}.txt" ]; then
    echo "Running BLAST: ${sp2} query against ${sp1} db..."
    $blast_cmd -query "${data_dir}/${sp2}.fa" -db "${data_dir}/${sp1}_db" \
        -evalue $evalue -outfmt $outfmt -num_threads $threads -max_target_seqs $max_target_seqs \
        -out "${result_dir}/${sp2}_vs_${sp1}.txt" || {
        echo "Error running BLAST for ${sp2} -> ${sp1}"
        exit 1
    }
else
    echo "BLAST result for ${sp2} -> ${sp1} already exists, skipping..."
fi

# ================= RBH 提取阶段 =================

echo "开始提取双向最佳比对 (RBH)..."

# 定义最终输出文件
output_file="${result_dir}/${sp1}_${sp2}_orthologs.txt"

# 写入表头
echo -e "${sp1}_protein_id\t${sp2}_protein_id" > "$output_file"

# 定义关联数组
declare -A forward_map

# 1. 读取 Forward 结果 (Chicken -> Mouse)
# 格式: query(chicken) subject(mouse) ...
# 存储: forward_map[chicken_id] = mouse_id
# 注意：sort -u -k1,1 确保每个 query 只取第一行（虽然 max_target_seqs=1 通常只给一个，但防万一）
echo "Loading forward results..."
while read -r q_id s_id; do
    forward_map["$q_id"]=$s_id
done < <(awk '{print $1 "\t" $2}' "${result_dir}/${sp1}_vs_${sp2}.txt")

# 2. 读取 Reverse 结果 (Mouse -> Chicken) 并验证 RBH
# 格式: query(mouse) subject(chicken) ...
# 逻辑: 如果 forward_map[chicken_id] == mouse_id，则是 RBH
echo "Checking reciprocal hits and writing output..."
count=0
while read -r q_id s_id; do
    # 这里 q_id 是 mouse (sp2), s_id 是 chicken (sp1)
    
    # 检查 Chicken -> Mouse 的结果是否也是 Mouse
    expected_mouse="${forward_map[$s_id]}"
    
    if [[ "$expected_mouse" == "$q_id" ]]; then
        # 匹配成功，写入文件 (输出顺序: Chicken ID, Mouse ID)
        echo -e "$s_id\t$q_id" >> "$output_file"
        ((count++))
    fi
done < <(awk '{print $1 "\t" $2}' "${result_dir}/${sp2}_vs_${sp1}.txt")

echo "完成！共找到 $count 对直系同源基因。"
echo "结果保存在: $output_file"
