#!/bin/bash

# 检查是否提供了输入文件
if [ $# -ne 1 ]; then
    echo "用法: $0 <input_list.txt>"
    echo "输入文件应包含三列（制表符分隔）："
    echo "1. 特征表文件绝对路径"
    echo "2. 物种名称"
    echo "3. 蛋白序列文件绝对路径"
    exit 1
fi

INPUT_FILE="$1"

# 检查输入文件是否存在
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 不存在"
    exit 1
fi

# 检查 parallel 命令是否可用
if ! command -v parallel &> /dev/null; then
    echo "警告: GNU parallel 未安装，将使用后台进程并行运行（可能不稳定）"
    USE_PARALLEL=false
else
    USE_PARALLEL=true
fi

# 检查 makeblastdb 命令是否可用
if ! command -v makeblastdb &> /dev/null; then
    echo "错误: makeblastdb 未安装，无法构建BLASTP索引"
    exit 1
fi

# 设置最大并行任务数（可根据机器性能调整）
MAX_JOBS=4

# 函数：执行单个任务
run_job() {
    feature_txt="$1"
    species="$2"
    protein_faa="$3"
    
    # 检查文件是否存在
    if [ ! -f "$feature_txt" ]; then
        echo "错误: 特征表文件 $feature_txt 不存在，跳过"
        return 1
    fi
    if [ ! -f "$protein_faa" ]; then
        echo "错误: 蛋白序列文件 $protein_faa 不存在，跳过"
        return 1
    fi

    # 运行 Python 脚本
    echo "运行: python /home/bzheng/project/11_multic_species_homologs/script/get_protein.py --input-feature-txt $feature_txt --species $species --input-protein-faa $protein_faa"
    python /home/bzheng/project/11_multic_species_homologs/script/get_allProteins_index.py --input-feature-txt "$feature_txt" --species "$species" --input-protein-faa "$protein_faa"
    if [ $? -ne 0 ]; then
        echo "错误: $species 的蛋白提取任务失败"
        return 1
    fi

    # 构建BLASTP索引
    output_fasta="${species}.fa"
    if [ -f "$output_fasta" ]; then
        echo "为 $species 构建BLASTP索引..."
        makeblastdb -in "$output_fasta" -dbtype prot -out "${species}_blastdb" -parse_seqids -title "${species} protein database"
        if [ $? -eq 0 ]; then
            echo "成功为 $species 构建BLASTP索引"
        else
            echo "错误: $species 的BLASTP索引构建失败"
            return 1
        fi
    else
        echo "错误: 输出文件 $output_fasta 未生成，无法构建索引"
        return 1
    fi
}

# 导出函数以便 parallel 使用
export -f run_job

# 主逻辑
if [ "$USE_PARALLEL" = true ]; then
    # 使用 GNU parallel 并行执行
    cat "$INPUT_FILE" | parallel --colsep '\t' --jobs "$MAX_JOBS" run_job {1} {2} {3}
else
    # 使用后台进程并行执行
    while IFS=$'\t' read -r feature_txt species protein_faa; do
        # 跳过空行
        [ -z "$feature_txt" ] && continue
        run_job "$feature_txt" "$species" "$protein_faa" &
        # 控制并行任务数量
        while [ $(jobs -r | wc -l) -ge "$MAX_JOBS" ]; do
            sleep 1
        done
    done < "$INPUT_FILE"
    # 等待所有后台任务完成
    wait
fi

echo "所有任务已完成！"
