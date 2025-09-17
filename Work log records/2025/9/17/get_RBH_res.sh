#!/bin/bash

# 设置序列数据文件数组
species=("japanese_quail" "water_buffalo" "european_hare" "rabbit" "fishing_cat"
         "clouded_leopard" "donkey" "equus_quagga" "raccoon_dog" "goat"
         "tufted_duck" "ruddy_duck" "goose" "northern_pintail" "duck"
         "ring_necked_pheasant" "rock_ptarmigan" "sumatran_serow" "red_fox"
         "common_warthog" "takin" "scimitar_horned_oryx" "bison_bison_bison" "cat"
         "horse" "turkey" "pig" "dog" "sheep" "cattle" "chicken")
# 输入的人类蛋白 ID 文件路径
human_id_file="/home/bzheng/project/11_multic_species_homologs/protein_index_202509/human_proteinID.txt"

# 结果目录
result_dir="/home/bzheng/project/11_multic_species_homologs/RBBH_res_202509"

# 创建最终输出文件，包含列名
output_file="${result_dir}/rbh_combined.txt"
echo -e "human_protein_id\t${species[*]// /_protein_id\t}_protein_id" > "$output_file"

# 读取人类蛋白 ID 文件并初始化关联数组
declare -A human_to_species
declare -A species_to_human

# 循环处理每个物种的 RBH
echo "开始处理BLAST结果文件..."
for sp in "${species[@]}"; do
    echo "处理物种: $sp"
    
    # 两个方向的 BLAST 结果文件
    human_to_sp="${result_dir}/human_ref_${sp}_query.txt"
    sp_to_human="${result_dir}/${sp}_ref_human_query.txt"
    
    # 检查文件是否存在
    if [[ ! -f "$human_to_sp" ]]; then
        echo "警告: 文件不存在 - $human_to_sp"
        continue
    fi
    if [[ ! -f "$sp_to_human" ]]; then
        echo "警告: 文件不存在 - $sp_to_human"
        continue
    fi
    
    # 使用 awk 提取 RBH 并存储到关联数组
    # 人类→物种的比对结果（BLAST输出格式6的第1列是query，第2列是target）
    while IFS=$'\t' read -r human_id sp_id rest; do
        human_to_species["$human_id:$sp"]=$sp_id
    done < "$human_to_sp"
    
    # 物种→人类的比对结果（BLAST输出格式6的第1列是query，第2列是target）
    while IFS=$'\t' read -r sp_id human_id rest; do
        species_to_human["$sp_id:$sp"]=$human_id
    done < "$sp_to_human"
done

# 处理人类蛋白 ID 文件并生成最终输出
echo "生成最终RBH结果..."

# 创建临时文件来存储提取的蛋白ID
temp_id_file="${result_dir}/temp_human_ids.txt"

# 提取human_proteinID.txt的第一列（蛋白ID）
awk '{print $1}' "$human_id_file" > "$temp_id_file"

# 读取提取的蛋白ID文件
while IFS= read -r human_id; do
    # 跳过空行和注释行
    [[ -z "$human_id" ]] && continue
    [[ "$human_id" =~ ^# ]] && continue
    
    # 初始化输出行，包含人类蛋白 ID
    output_line="$human_id"
    
    # 检查每个物种的 RBH
    for sp in "${species[@]}"; do
        sp_id="${human_to_species[$human_id:$sp]}"
        
        # 验证是否为 RBH (双向最佳匹配)
        if [[ -n "$sp_id" && "${species_to_human[$sp_id:$sp]}" == "$human_id" ]]; then
            output_line="$output_line\t$sp_id"
        else
            output_line="$output_line\t-"
        fi
    done
    
    # 将结果写入输出文件
    echo -e "$output_line" >> "$output_file"
done < "$temp_id_file"

# 清理临时文件
rm -f "$temp_id_file"

echo "RBH分析完成! 结果保存在: $output_file"
echo "总共处理了 $(wc -l < "$output_file") 行数据（包括表头）"
