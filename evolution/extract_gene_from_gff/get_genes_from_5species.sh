#!/bin/bash

# 定义物种GFF路径（关联数组）
declare -A GFF_PATHS
GFF_PATHS["chicken"]="/home/bzheng/ref/huxu_v23_ref/20250514_final_ref/20250529_chicken.gtf"
GFF_PATHS["dog"]="/home/bzheng/others_species_analysis/dog/ref/dog.gff"
GFF_PATHS["cattle"]="/home/bzheng/others_species_analysis/cattle/ref/Btau_5.0.1_new.gff"
GFF_PATHS["pig"]="/home/bzheng/others_species_analysis/pig/ref/pig.gff"
GFF_PATHS["sheep"]="/home/bzheng/others_species_analysis/sheep/ref/sheep.gff"

# 默认参数
INPUT_FILE=""
SPECIES="cattle"
OUTPUT_DIR="./"  # 默认输出到当前目录
SHOW_HELP=false

# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input) INPUT_FILE="$2"; shift ;;
        --species) SPECIES="$2"; shift ;;
        --output) OUTPUT_DIR="$2"; shift ;;
        -h|--help) SHOW_HELP=true ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
    shift
done

# 帮助信息
if $SHOW_HELP; then
    echo "用法: $0 --input <输入文件> [选项]"
    echo "选项:"
    echo "  --input    输入文件路径（必需）"
    echo "  --species  物种（cattle/dog/pig/sheep/chicken）[默认: cattle]"
    echo "  --output   输出目录路径 [默认: 当前目录]"
    exit 0
fi

# 输入校验
if [ -z "$INPUT_FILE" ]; then
    echo "错误: 必须指定输入文件" >&2
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件不存在: $INPUT_FILE" >&2
    exit 1
fi

if [ -z "${GFF_PATHS[$SPECIES]}" ]; then
    echo "错误: 不支持的物种，可用选项: ${!GFF_PATHS[@]}" >&2
    exit 1
fi

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 获取输入文件基名（不含扩展名）
INPUT_BASENAME=$(basename "$INPUT_FILE" | cut -d. -f1)

# 定义窗口大小
WINDOWS=("0kb" "50kb")

# 处理每个窗口
for WINDOW in "${WINDOWS[@]}"; do
    OUTPUT_FILE="${OUTPUT_DIR}/${INPUT_BASENAME}_${WINDOW}_genes.txt"
    
    # 调用Python脚本
    python Get_gene_from_GTF-or-GFF_202506.py \
        "$INPUT_FILE" \
        "${GFF_PATHS[$SPECIES]}" \
        "$WINDOW" \
        "$OUTPUT_FILE"
    
    echo "已生成: $OUTPUT_FILE"
done
