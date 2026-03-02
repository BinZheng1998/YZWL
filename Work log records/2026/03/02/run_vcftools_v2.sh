#!/bin/bash

# 初始化变量
VCF=""
POP1=""
POP2=""
SPECIES=""

# 解析参数的函数
usage() {
    echo "使用方法: $0 --vcf <file.vcf.gz> --pop1_highBW <pop1.txt> --pop2_lowBW <pop2.txt> --species <name>"
    echo "参数说明:"
    echo "  --vcf      输入的压缩 VCF 文件 (.vcf.gz)"
    echo "  --pop1_highBW     群体 1 的样本列表文件 (每行一个 ID)"
    echo "  --pop2_lowBW     群体 2 的样本列表文件 (每行一个 ID)"
    echo "  --species  物种名称 (如 chicken, pig)，将作为输出文件唯一前缀"
    exit 1
}

# 循环解析参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --vcf) VCF="$2"; shift ;;
        --pop1_highBW) POP1="$2"; shift ;;
        --pop2_lowBW) POP2="$2"; shift ;;
        --species) SPECIES="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "未知参数: $1"; usage ;;
    esac
    shift
done

# 检查必要参数
if [[ -z "$VCF" || -z "$POP1" || -z "$POP2" || -z "$SPECIES" ]]; then
    echo "错误: 必须提供 --vcf, --pop1_higbBW, --pop2_lowBW 和 --species 参数。"
    usage
fi

# 设定输出前缀仅为物种名
OUT_PREFIX="${SPECIES}"

# 定义窗口大小组合: "窗口大小:步长"
WINDOWS=("5000:2500" "10000:5000" "20000:10000" "50000:25000")

for WIN_SET in "${WINDOWS[@]}"; do
    SIZE=$(echo $WIN_SET | cut -d: -f1)
    STEP=$(echo $WIN_SET | cut -d: -f2)
    TAG_SIZE=$(awk "BEGIN {print $SIZE/1000}")
    TAG_STEP=$(awk "BEGIN {print $STEP/1000}")
    TAG="${TAG_SIZE}kwindow_${TAG_STEP}kstep"
   #TAG="$((SIZE/1000))kwindow_$((STEP/1000))kstep"

    echo "===================================================="
    echo "物种: $SPECIES | 正在计算: [$TAG] 窗口"
    echo "===================================================="

    # 1. 计算 Fst
    vcftools --gzvcf "$VCF" --weir-fst-pop "$POP1" --weir-fst-pop "$POP2" \
        --fst-window-size "$SIZE" --fst-window-step "$STEP" \
        --out "./${OUT_PREFIX}_${TAG}_Fst"

    # 2. 计算 Pi
    vcftools --gzvcf "$VCF" --keep "$POP1" --window-pi "$SIZE" --window-pi-step "$STEP" \
        --out "./${OUT_PREFIX}_${TAG}_Pop1_highBW_Pi"
    vcftools --gzvcf "$VCF" --keep "$POP2" --window-pi "$SIZE" --window-pi-step "$STEP" \
        --out "./${OUT_PREFIX}_${TAG}_Pop2_lowBW_Pi"

    # 3. 计算 Tajima's D
    vcftools --gzvcf "$VCF" --keep "$POP1" --TajimaD "$SIZE" \
        --out "./${OUT_PREFIX}_${TAG}_Pop1_highBW_TajimaD"
    vcftools --gzvcf "$VCF" --keep "$POP2" --TajimaD "$SIZE" \
        --out "./${OUT_PREFIX}_${TAG}_Pop2_lowBW_TajimaD"
done

# 4. 计算全位点等位基因频率
echo ">>> 正在计算 ${SPECIES} 的全位点等位基因频率..."
vcftools --gzvcf "$VCF" --freq --keep "$POP1" --out "./${OUT_PREFIX}_Pop1_highBW_AlleleFreq"
vcftools --gzvcf "$VCF" --freq --keep "$POP2" --out "./${OUT_PREFIX}_Pop2_lowBW_AlleleFreq"

echo "===================================================="
echo "任务完成！"
echo "生成文件示例: ${OUT_PREFIX}_20k_10k_Fst.windowed.weir.fst"
