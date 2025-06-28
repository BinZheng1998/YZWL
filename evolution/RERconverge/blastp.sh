#!/bin/bash

# 设置序列数据文件数组（不包括 human，因为只需要其他物种）
species=("pig" "dog" "sheep" "cattle" "chicken" "aardvark" "armadillo" "big_brown_bat" "black_flying_fox" "cape_elephant_shrew" "cape_golden_mole" "chinese_tree_shrew" "coelacanth" "collared_flycatcher" "european_shrew" "green_anole" "green_monkey" "large_flying_fox" "lesser_Egyptian_jerboa" "medium_ground_finch" "opossum" "small_eared_galago" "star_nosed_mole" "tasmanian_devil" "tibetan_antelope" "tibetan_ground_jay" "tropical_clawed_frog" "white_throated_sparrow") # 可添加任意数量的物种，例如 ("pig" "dog" "sheep" "cattle" "chicken" "horse" "cat")

# 设置 BLAST 命令的基本路径和参数
blast_cmd="/home/bzheng/software/ncbi-blast-2.15.0+/bin/blastp"
evalue="1e-4"
outfmt="6"
threads="128"
max_target_seqs="1"

# 输入的人类蛋白 ID 文件路径
human_id_file="/home/bzheng/project/10_RER/data/human_protein_id.txt" # 请替换为实际的人类蛋白 ID 文件路径

# 创建结果目录
mkdir -p /home/bzheng/project/10_RER/result

# 进行 BLAST 比对，仅针对 human 和其他物种
for sp in "${species[@]}"; do
    # human 查询，物种作为参考数据库
    $blast_cmd -query "/home/bzheng/project/10_RER/data/human.fa" -db "/home/bzheng/project/10_RER/data/${sp}_db" \
        -evalue $evalue -outfmt $outfmt -num_threads $threads -max_target_seqs $max_target_seqs \
        -out "/home/bzheng/project/10_RER/result/human_ref_${sp}_query.txt" || {
        echo "Error running BLAST for human query, ${sp} database"
        exit 1
    }

    # 物种查询，human 作为参考数据库
    $blast_cmd -query "/home/bzheng/project/10_RER/data/${sp}.fa" -db "/home/bzheng/project/10_RER/data/human_db" \
        -evalue $evalue -outfmt $outfmt -num_threads $threads -max_target_seqs $max_target_seqs \
        -out "/home/bzheng/project/10_RER/result/${sp}_ref_human_query.txt" || {
        echo "Error running BLAST for ${sp} query, human database"
        exit 1
    }
done

# 创建最终输出文件，包含列名
output_file="/home/bzheng/project/10_RER/result/rbh_combined.txt"
echo -e "human_protein_id\t${species[*]// /_protein_id\t}_protein_id" > "$output_file"

# 读取人类蛋白 ID 文件并初始化关联数组
declare -A human_to_species
declare -A species_to_human

# 循环处理每个物种的 RBH
for sp in "${species[@]}"; do
    # 读取两个方向的 BLAST 结果
    human_to_sp="/home/bzheng/project/10_RER/result/human_ref_${sp}_query.txt"
    sp_to_human="/home/bzheng/project/10_RER/result/${sp}_ref_human_query.txt"

    # 使用 awk 提取 RBH 并存储到关联数组
    while read -r human_id sp_id; do
        human_to_species["$human_id:$sp"]=$sp_id
    done < <(awk '{print $1 "\t" $2}' "$human_to_sp")

    while read -r sp_id human_id; do
        species_to_human["$sp_id:$sp"]=$human_id
    done < <(awk '{print $1 "\t" $2}' "$sp_to_human")
done

# 处理人类蛋白 ID 文件并生成最终输出
while IFS= read -r human_id; do
    # 跳过空行
    [[ -z "$human_id" ]] && continue

    # 初始化输出行，包含人类蛋白 ID
    output_line="$human_id"

    # 检查每个物种的 RBH
    for sp in "${species[@]}"; do
        sp_id="${human_to_species[$human_id:$sp]}"
        # 验证是否为 RBH
        if [[ -n "$sp_id" && "${species_to_human[$sp_id:$sp]}" == "$human_id" ]]; then
            output_line="$output_line\t$sp_id"
        else
            output_line="$output_line\t-"
        fi
    done

    # 将结果写入输出文件
    echo -e "$output_line" >> "$output_file"
done < "$human_id_file"
