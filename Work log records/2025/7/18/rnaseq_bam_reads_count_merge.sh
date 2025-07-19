#!/bin/bash

# 创建临时文件保存所有输入文件列表
file_list=$(mktemp)
find . -maxdepth 1 -type f -name "*_reads_counts.txt" > "$file_list"

# 检查是否找到文件
if [ ! -s "$file_list" ]; then
    echo "Error: No *_reads_counts.txt files found in current directory"
    exit 1
fi

# 获取所有样本名称
sample_names=$(sed 's|^\./||; s/_reads_counts\.txt$//' "$file_list")

# 获取第一个文件的前三列（染色体、开始、结束位置）
first_file=$(head -1 "$file_list" | tr -d '\n')
awk '{print $1 "\t" $2 "\t" $3}' "$first_file" > merged_counts.txt

# 添加每个样本的reads列
while IFS= read -r file; do
    # 提取样本名称
    sample_name=$(basename "$file" | sed 's/_reads_counts\.txt$//')
    
    # 添加当前文件的reads列
    paste merged_counts.txt <(awk '{print $4}' "$file") > merged_counts.tmp
    mv merged_counts.tmp merged_counts.txt
    
    # 添加列名到标题行
    sed -i "1s/$/\t$sample_name/" merged_counts.txt
done < "$file_list"

# 计算总和并添加到最后一列
awk '
NR == 1 {
    # 打印标题行并添加"Total"列
    print $0, "Total"
    next
}
{
    total = 0
    # 从第4列开始求和
    for (i=4; i<=NF; i++) {
        total += $i
    }
    # 打印整行并添加总和
    print $0, total
}' merged_counts.txt > merged_counts_final.txt

# 清理临时文件
rm "$file_list"

echo "Combined file created: merged_counts_final.txt"
