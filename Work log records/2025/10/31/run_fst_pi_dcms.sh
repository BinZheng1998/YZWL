#!/bin/bash

# 默认参数
VCF=""
GROUP_FILE=""
RSCRIPT=""
SPECIES="chicken"
OUTPUT_FOLDER="./results"
WINDOW_SIZE=10000
STEP_SIZE=5000
JOBS=$(nproc)  # 默认使用所有CPU核心

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --vcf)
            VCF="$2"
            shift 2
            ;;
        --group)
            GROUP_FILE="$2"
            shift 2
            ;;
        --rscript)
            RSCRIPT="$2"
            shift 2
            ;;
        --species)
            SPECIES="$2"
            shift 2
            ;;
        --output-folder)
            OUTPUT_FOLDER="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# 检查必要参数
if [[ -z "$VCF" || -z "$GROUP_FILE" || -z "$RSCRIPT" ]]; then
    echo "Usage: $0 --vcf <vcf_file> --group <group_file> --rscript <r_script> [--species <species>] [--output-folder <folder>] [--jobs <num_jobs>]"
    echo ""
    echo "Required parameters:"
    echo "  --vcf          Input VCF file (.vcf or .vcf.gz)"
    echo "  --group        Sample information file (tab-separated with header)"
    echo "  --rscript      R script for analysis"
    echo ""
    echo "Optional parameters:"
    echo "  --species      Species name (default: chicken)"
    echo "  --output-folder Output directory (default: ./results)"
    echo "  --jobs         Number of parallel jobs (default: all CPU cores)"
    exit 1
fi

# 检查输入文件是否存在
if [[ ! -f "$VCF" ]]; then
    echo "Error: VCF file '$VCF' not found!"
    exit 1
fi

if [[ ! -f "$GROUP_FILE" ]]; then
    echo "Error: Group file '$GROUP_FILE' not found!"
    exit 1
fi

if [[ ! -f "$RSCRIPT" ]]; then
    echo "Error: R script '$RSCRIPT' not found!"
    exit 1
fi

# 检查是否安装了parallel
if ! command -v parallel &> /dev/null; then
    echo "Error: GNU Parallel is not installed. Please install it first."
    echo "On Ubuntu/Debian: sudo apt-get install parallel"
    echo "On CentOS/RHEL: sudo yum install parallel"
    exit 1
fi

echo "Starting PARALLEL Fst and Pi analysis with sliding windows..."
echo "VCF file: $VCF"
echo "Group file: $GROUP_FILE"
echo "R script: $RSCRIPT"
echo "Species: $SPECIES"
echo "Output folder: $OUTPUT_FOLDER"
echo "WindowOUTPUT_FOLDER"
echo "Window size: $WINDOW_SIZE bp"
echo "Step size: $STEP_SIZE bp"
echo "Parallel jobs: $JOBS"

# 创建输出目录
mkdir -p "$OUTPUT_FOLDER"

# 创建临时工作目录
TEMP_DIR=$(mktemp -d -p "$OUTPUT_FOLDER" temp.XXXXXX)
echo "Temporary directory: $TEMP_DIR"

# 函数：清理临时文件
cleanup() {
    echo "Cleaning up temporary files..."
    rm -rf "$TEMP_DIR"
}
trap cleanup EXIT

# 解析样本信息文件
process_sample_info() {
    local group_file="$1"
    
    # 清空现有文件
    > "$TEMP_DIR/high_varieties.txt"
    > "$TEMP_DIR/low_varieties.txt"
    
    # 读取样本信息（跳过第一行列名）
    {
        read header;
        while IFS=$'\t' read -r sample variety phenotype; do
            # 去除可能的回车符
            sample=$(echo "$sample" | tr -d '\r')
            variety=$(echo "$variety" | tr -d '\r') 
            phenotype=$(echo "$phenotype" | tr -d '\r')
            
            if [[ "$phenotype" == "high" ]]; then
                echo -e "$sample\t$variety" >> "$TEMP_DIR/high_varieties.txt"
            elif [[ "$phenotype" == "low" ]]; then
                echo -e "$sample\t$variety" >> "$TEMP_DIR/low_varieties.txt"
            else
                echo "Warning: Unknown phenotype '$phenotype' for sample $sample"
            fi
        done
    } < "$group_file"
}

# 获取唯一品种列表
get_unique_varieties() {
    local file="$1"
    awk -F'\t' '{print $2}' "$file" | sort | uniq
}

# 为特定品种提取样本列表
extract_samples_by_variety() {
    local varieties_file="$1"
    local target_variety="$2"
    awk -F'\t' -v var="$target_variety" '$2 == var {print $1}' "$varieties_file"
}

# 单次比较的函数 - 将被并行调用
single_comparison() {
    local high_var="$1"
    local low_var="$2"
    local vcf="$3"
    local rscript="$4"
    local species="$5"
    local output_folder="$6"
    local window_size="$7"
    local step_size="$8"
    
    # 创建当前比较的前缀
    safe_high_var=$(echo "$high_var" | sed 's/[^a-zA-Z0-9]/_/g')
    safe_low_var=$(echo "$low_var" | sed 's/[^a-zA-Z0-9]/_/g')
    current_prefix="${safe_high_var}_vs_${safe_low_var}"
    
    # 提取样本
    extract_samples_by_variety "$TEMP_DIR/high_varieties.txt" "$high_var" > "${TEMP_DIR}/high_samples_${safe_high_var}_${safe_low_var}.txt"
    extract_samples_by_variety "$TEMP_DIR/low_varieties.txt" "$low_var" > "${TEMP_DIR}/low_samples_${safe_high_var}_${safe_low_var}.txt"
    
    # 检查是否有足够的样本
    high_count=$(wc -l < "${TEMP_DIR}/high_samples_${safe_high_var}_${safe_low_var}.txt")
    low_count=$(wc -l < "${TEMP_DIR}/low_samples_${safe_high_var}_${safe_low_var}.txt")
    
    if [[ $high_count -eq 0 ]] || [[ $low_count -eq 0 ]]; then
        echo "  Warning: Insufficient samples for $high_var vs $low_var, skipping..."
        return 1
    fi
    
    echo "  [$BASHPID] Processing $high_var ($high_count samples) vs $low_var ($low_count samples)"
    
    # 计算窗口Fst
    vcftools --gzvcf "$vcf" \
             --weir-fst-pop "${TEMP_DIR}/high_samples_${safe_high_var}_${safe_low_var}.txt" \
             --weir-fst-pop "${TEMP_DIR}/low_samples_${safe_high_var}_${safe_low_var}.txt" \
             --fst-window-size "$window_size" \
             --fst-window-step "$step_size" \
             --out "$current_prefix" 2>/dev/null
    
    # 计算高表达组的窗口Pi
    vcftools --gzvcf "$vcf" \
             --keep "${TEMP_DIR}/high_samples_${safe_high_var}_${safe_low_var}.txt" \
             --window-pi "$window_size" \
             --window-pi-step "$step_size" \
             --out "${current_prefix}_high" 2>/dev/null
    
    # 计算低表达组的窗口Pi
    vcftools --gzvcf "$vcf" \
             --keep "${TEMP_DIR}/low_samples_${safe_high_var}_${safe_low_var}.txt" \
             --window-pi "$window_size" \
             --window-pi-step "$step_size" \
             --out "${current_prefix}_low" 2>/dev/null
    
    # 设置输出文件名
    fst_file="${current_prefix}.windowed.weir.fst"
    pi1_file="${current_prefix}_high.windowed.pi"
    pi2_file="${current_prefix}_low.windowed.pi"
    
    # 如果所有输出文件都存在，调用R脚本
    if [[ -f "$fst_file" ]] && [[ -f "$pi1_file" ]] && [[ -f "$pi2_file" ]]; then
        /usr/bin/Rscript "$rscript" \
            --input-fst="$fst_file" \
            --input-pi1="$pi1_file" \
            --input-pi2="$pi2_file" \
            --species="$species" \
            --out-prefix="$current_prefix" \
            --output-folder="$output_folder"
        
        # 清理本次比较的临时文件
        rm "${TEMP_DIR}/high_samples_${safe_high_var}_${safe_low_var}.txt" \
           "${TEMP_DIR}/low_samples_${safe_high_var}_${safe_low_var}.txt"
        
        echo "  [$BASHPID] Completed $high_var vs $low_var"
    else
        echo "  [$BASHPID] Failed to generate all output files for $high_var vs $low_var"
        return 1
    fi
    
    return 0
}

export -f single_comparison
export -f extract_samples_by_variety

# 处理样本信息
echo "Processing sample information..."
process_sample_info "$GROUP_FILE"

# 获取唯一的品种列表
readarray -t HIGH_VARIETIES < <(get_unique_varieties "$TEMP_DIR/high_varieties.txt")
readarray -t LOW_VARIETIES < <(get_unique_varieties "$TEMP_DIR/low_varieties.txt")

echo "Found ${#HIGH_VARIETIES[@]} varieties in high group: ${HIGH_VARIETIES[*]}"
echo "Found ${#LOW_VARIETIES[@]} varieties in low group: ${LOW_VARIETIES[*]}"

# 生成所有比较对的列表
COMPARISON_LIST="$TEMP_DIR/comparison_list.txt"
> "$COMPARISON_LIST"

for high_var in "${HIGH_VARIETIES[@]}"; do
    for low_var in "${LOW_VARIETIES[@]}"; do
        echo "$high_var"$'\t'"$low_var" >> "$COMPARISON_LIST"
    done
done

total_comparisons=$(wc -l < "$COMPARISON_LIST")
echo "Total comparisons to perform: $total_comparisons"

# 使用parallel并行执行
echo "Starting parallel execution with $JOBS jobs..."

# 导出必要的变量给子进程
export VCF
export RSCRIPT
export SPECIES
export OUTPUT_FOLDER
export WINDOW_SIZE
export STEP_SIZE
export TEMP_DIR

# 执行并行计算
cat "$COMPARISON_LIST" | parallel --progress --bar --jobs "$JOBS" \
    --colsep '\t' \
    single_comparison "{1}" "{2}" "$VCF" "$RSCRIPT" "$SPECIES" "$OUTPUT_FOLDER" "$WINDOW_SIZE" "$STEP_SIZE"

echo ""
echo "=========================================="
echo "PARALLEL analysis completed!"
echo "Total comparisons: $total_comparisons"
echo "Results are saved in: $OUTPUT_FOLDER"
echo "=========================================="
