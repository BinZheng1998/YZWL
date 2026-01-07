#!/bin/bash
# =============================================================================
# Modern Fast ATAC-seq Pipeline (Split Version)
# Logic: Chromap -> .sam file -> Samtools Filter -> .bam file
# =============================================================================

set -e

# ===================== 1. 配置区域 =====================

PROJECT_DIR="$(pwd)"
RAW_DATA_DIR="${PROJECT_DIR}"
RESULTS_DIR="${PROJECT_DIR}/results_fast_encode"

# 请确认这里是修正后的绝对路径
REF_FASTA="/home/bzheng/others_species_analysis/pig/ref/pig.fa"
CHROMAP_INDEX="/home/bzheng/project/05_epigenetics/ATAC/pig/0.ref/chromap/pig_chromap_index"

GENOME_SIZE="2.5e9"
THREADS=12

echo ">>> Scanning for samples in current directory..."
SAMPLES=($(find . -maxdepth 1 -name "*_1.fastq.gz" | sed 's/^\.\///' | sed 's/_1.fastq.gz//' | sort -u))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "ERROR: No *_1.fastq.gz files found!"
    exit 1
fi

# ===================== 2. 目录准备 =====================

DIR_CLEAN="${RESULTS_DIR}/01_clean"
DIR_ALIGN="${RESULTS_DIR}/02_alignment"
DIR_PEAKS="${RESULTS_DIR}/03_peaks"
DIR_LOGS="${RESULTS_DIR}/logs"

mkdir -p ${DIR_CLEAN} ${DIR_ALIGN} ${DIR_PEAKS} ${DIR_LOGS}

# ===================== 3. 主循环 =====================

for SAMPLE in "${SAMPLES[@]}"; do
    echo "======================================================="
    echo "Processing Sample: ${SAMPLE}"
    echo "======================================================="

    R1="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"
    R2="${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz"
    
    CLEAN_R1="${DIR_CLEAN}/${SAMPLE}_1.clean.fq.gz"
    CLEAN_R2="${DIR_CLEAN}/${SAMPLE}_2.clean.fq.gz"
    
    # 定义中间 SAM 文件和最终 BAM 文件
    SAM_RAW="${DIR_ALIGN}/${SAMPLE}.raw.sam"
    BAM_FINAL="${DIR_ALIGN}/${SAMPLE}.final.bam"

    # --- Step 1: fastp ---
    if [ ! -f "${CLEAN_R1}" ]; then
        echo "[1/4] QC & Trimming with fastp..."
        fastp -i ${R1} -I ${R2} \
              -o ${CLEAN_R1} -O ${CLEAN_R2} \
              --detect_adapter_for_pe \
              --thread 12 \
              -h ${DIR_LOGS}/${SAMPLE}.fastp.html \
              -j ${DIR_LOGS}/${SAMPLE}.fastp.json \
              2> ${DIR_LOGS}/${SAMPLE}.fastp.log
    else
        echo "[1/4] Fastp results exist. Skipping."
    fi

    # --- Step 2a: Chromap Align (只比对，输出 SAM) ---
    # 如果最终 BAM 不存在，且中间 SAM 也不存在，就开始比对
    if [ ! -f "${BAM_FINAL}" ]; then
        if [ ! -f "${SAM_RAW}" ]; then
            echo "[2a/4] Chromap Mapping to SAM..."
            
            # 注意：这里去掉了 | 管道，直接输出到 -o ${SAM_RAW}
            /home/bzheng/software/chromap/chromap --preset atac \
                    -r ${REF_FASTA} \
                    -x ${CHROMAP_INDEX} \
                    -1 ${CLEAN_R1} -2 ${CLEAN_R2} \
                    --remove-pcr-duplicates \
                    -q 25 \
                    -t ${THREADS} \
                    --SAM \
                    -o ${SAM_RAW} \
                    2> ${DIR_LOGS}/${SAMPLE}.chromap.log
        else
            echo "[2a/4] Raw SAM exists. Skipping mapping."
        fi

        # --- Step 2b: Samtools Filtering (过滤 + 排序) ---
        echo "[2b/4] Filtering (MAPQ>25, no chrM) & Sorting..."
        
        # 现在的输入是 ${SAM_RAW}
        samtools view -h -q 25 -F 1804 ${SAM_RAW} | \
        grep -v "chrM" | grep -v "MT" | \
        samtools sort -@ 4 -o ${BAM_FINAL} -
        
        samtools index ${BAM_FINAL}
        
        # --- Step 2c: 清理巨大的 SAM 文件 ---
        # 这一步很重要！SAM 文件比 BAM 大 4-5 倍，不删会撑爆硬盘
        echo "[Cleanup] Removing intermediate SAM file..."
        rm -f ${SAM_RAW}
        
    else
        echo "[2/4] Final BAM exists. Skipping alignment."
    fi

    # --- Step 3: MACS3 Peak Calling ---
    if [ ! -f "${DIR_PEAKS}/${SAMPLE}_peaks.narrowPeak" ]; then
        echo "[3/4] Calling peaks with MACS3..."
        macs3 callpeak \
            -t ${BAM_FINAL} \
            -n ${SAMPLE} \
            -g ${GENOME_SIZE} \
            --outdir ${DIR_PEAKS} \
            -p 0.01 \
            --nomodel \
            --shift -75 \
            --extsize 150 \
            -B --SPMR \
            --keep-dup all \
            --call-summits \
            2> ${DIR_LOGS}/${SAMPLE}.macs3.log
    else
        echo "[3/4] Peaks exist. Skipping."
    fi

done

echo ">>>> Analysis Pipeline Finished! <<<<"
