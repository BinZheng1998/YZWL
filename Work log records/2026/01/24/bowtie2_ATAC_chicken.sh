#!/bin/bash

# ================= 配置区域 =================

# 1. 软件路径
PICARD_JAR="/home/bzheng/software/picard.jar"
# 这里填写您的 Batu5.0.1 索引前缀
REF_INDEX="/home/bzheng/project/05_epigenetics/ref/chicken_index" 

# 2. 核心参数
THREADS=8
MAX_FRAGMENT_SIZE=2000

# 3. 线粒体过滤名单 (正则匹配)
MITO_REGEX="^(chrM|chrMT|MT|mt|M)$"

# ===========================================

# 检查软件
if [ ! -f "$PICARD_JAR" ]; then echo "Error: 找不到 Picard"; exit 1; fi
mkdir -p results/logs results/bams results/metrics

echo ">>> [ATAC-seq Mode] 开始扫描 _*.clean.fastq.gz 文件..."

# 识别 _1.clean.fastq.gz 结尾的文件并提取样本名
SAMPLES=$(ls *_1.clean.fastq.gz 2>/dev/null | sed 's/_1.clean.fastq.gz$//' | sort | uniq)

if [ -z "$SAMPLES" ]; then
    echo "当前目录下未找到 *_1.clean.fastq.gz 文件！"
    exit 1
fi

for SAMPLE in $SAMPLES; do
    echo "========================================================"
    echo "正在处理样本: $SAMPLE"
    echo "开始时间: $(date)"

    # 构造文件名
    R1="${SAMPLE}_1.clean.fastq.gz"
    R2="${SAMPLE}_2.clean.fastq.gz"
    
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "警告: 找不到样本 $SAMPLE 的配对文件 ($R2)，跳过。"
        continue
    fi

    # 定义输出文件路径
    SORTED_BAM="results/bams/${SAMPLE}.sorted.bam"       # 中间文件 (很大, 稍后删除)
    METRICS_FILE="results/metrics/${SAMPLE}.marked_dup_metrics.txt"
    FINAL_BAM="results/bams/${SAMPLE}.atac_final.bam"    # 最终结果
    LOG_FILE="results/logs/${SAMPLE}.log"

    # ================= 核心处理管道 =================
    
    echo ">>> Step 1: Bowtie2 Alignment..."
    # Bowtie2 -> samtools view (丢弃 unmapped) -> samtools sort
    /usr/bin/time -v bowtie2 --threads $THREADS \
        -X $MAX_FRAGMENT_SIZE \
        --very-sensitive --no-mixed --no-discordant \
        -x "$REF_INDEX" -1 "$R1" -2 "$R2" 2>> "$LOG_FILE" | \
    samtools view -h -b -F 4 | \
    samtools sort --threads $(($THREADS - 1)) -o "$SORTED_BAM"

    if [ -f "$SORTED_BAM" ]; then
        echo ">>> Step 2: MarkDuplicates & Filtering..."
        
        # Picard -> samtools filter -> awk mito remove -> final bam
        /usr/bin/time -v /home/bzheng/software/jdk-22.0.1/bin/java -jar "$PICARD_JAR" MarkDuplicates \
            INPUT="$SORTED_BAM" \
            OUTPUT=/dev/stdout \
            METRICS_FILE="$METRICS_FILE" \
            VALIDATION_STRINGENCY=LENIENT \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=false \
            QUIET=true \
            COMPRESSION_LEVEL=0 2>> "$LOG_FILE" | \
        samtools view -h -F 1804 -f 2 -q 30 /dev/stdin | \
        awk -v mito="$MITO_REGEX" 'BEGIN {OFS="\t"} /^@/ {print $0; next} $3 !~ mito {print $0}' | \
        samtools view -b - > "$FINAL_BAM"
        
        # 建立索引
        samtools index "$FINAL_BAM"
        
        # --- 【关键修改】安全删除中间文件 ---
        # 逻辑：检查最终文件是否存在(-f) 且 大小不为0(-s)
        if [ -s "$FINAL_BAM" ]; then
            echo ">>> 最终文件生成成功，正在删除巨大的中间文件: $SORTED_BAM"
            rm "$SORTED_BAM"
        else
            echo ">>> Error: 最终 BAM 文件似乎有问题，保留中间文件以供检查！"
        fi
        
        echo ">>> $SAMPLE 处理完成！"
    else
        echo ">>> Error: Step 1 比对失败，请查看 $LOG_FILE"
    fi

done

echo "========================================================"
echo "所有任务完成。空间已释放。"
