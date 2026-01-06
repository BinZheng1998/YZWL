#!/bin/bash
# =============================================================================
# RNA-seq Analysis Pipeline (Pig Data) - Auto Cleanup Version
# 流程: Trimmomatic -> STAR -> StringTie -> featureCounts
# 特性: 自动识别SE/PE、自动断点续跑、即时清理中间文件节省空间
# =============================================================================

set -e

# ===================== 1. 核心配置 =====================

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DATA_DIR="${PROJECT_DIR}"
RESULTS_DIR="${PROJECT_DIR}/results"

# --- 参考基因组 ---
REF_DIR="/home/bzheng/project/04_rnaseq/pig/0.ref/reference_genome"
REF_GENOME_FA="${REF_DIR}/GCF_000003025.6_Sscrofa11.1_genomic.renamed.fna"
GTF_FILE="${REF_DIR}/Sus_scrofa.Sscrofa11.1.100.gtf"
STAR_INDEX_DIR="${REF_DIR}/star_index"

# --- 软件路径 ---
TRIMMOMATIC_JAR="/home/bzheng/.conda/envs/pigGTEx/share/trimmomatic-0.39-2/trimmomatic.jar"
ADAPTER_DIR="/home/bzheng/.conda/envs/pigGTEx/share/trimmomatic-0.39-2/adapters"

# --- 并行策略 ---
MAX_PARALLEL_SAMPLES=8  # 建议保持为4，防止内存溢出
THREADS_PER_JOB=10
THREADS_FC=24

# --- 输出目录 ---
TRIMMED_DIR="${RESULTS_DIR}/trimmed_data"
MAPPED_DIR="${RESULTS_DIR}/mapped_data"
COUNTS_DIR="${RESULTS_DIR}/gene_counts"
STRINGTIE_DIR="${RESULTS_DIR}/stringtie_abundance"

mkdir -p ${TRIMMED_DIR}/logs ${MAPPED_DIR}/logs ${COUNTS_DIR} ${STRINGTIE_DIR}/logs

# ===================== 2. 样本扫描 =====================

cd ${RAW_DATA_DIR}
# 使用 find 查找，支持 .fastq.gz 结尾的所有文件，自动识别名字
SAMPLES=($(find . -maxdepth 1 -name "*.fastq.gz" | \
    sed 's|^\./||' | \
    sed -E 's/(_1|_2)?\.fastq\.gz//' | \
    sort -u))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: No .fastq.gz files found."
    exit 1
fi

echo ">>> Found ${#SAMPLES[@]} samples."

# ===================== 3. 处理函数 (含即时清理) =====================

process_sample_parallel() {
    local SAMPLE=$1
    local SAMPLE_MAPPED_DIR="${MAPPED_DIR}/${SAMPLE}"
    local BAM_FILE="${SAMPLE_MAPPED_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    local STRINGTIE_ABUND="${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt"

    # --- 1. 自动判断 SE/PE ---
    if [ -f "${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz" ]; then
        MODE="PE"
        R1="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"
        R2="${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz"
        TRIM_R1="${TRIMMED_DIR}/${SAMPLE}_1.clean.fq.gz"
        TRIM_R2="${TRIMMED_DIR}/${SAMPLE}_2.clean.fq.gz"
        # PE Trimmomatic & STAR Args
        TRIM_CMD="PE -phred33 ${R1} ${R2} ${TRIM_R1} ${TRIMMED_DIR}/${SAMPLE}_1.unpaired.fq.gz ${TRIM_R2} ${TRIMMED_DIR}/${SAMPLE}_2.unpaired.fq.gz ILLUMINACLIP:${ADAPTER_DIR}/TruSeq3-PE-2.fa:2:30:10:2:TRUE"
        STAR_INPUT="${TRIM_R1} ${TRIM_R2}"
    else
        MODE="SE"
        if [ -f "${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz" ]; then R1="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"; else R1="${RAW_DATA_DIR}/${SAMPLE}.fastq.gz"; fi        TRIM_R1="${TRIMMED_DIR}/${SAMPLE}.clean.fq.gz"
        # SE Trimmomatic & STAR Args
        TRIM_CMD="SE -phred33 ${R1} ${TRIM_R1} ILLUMINACLIP:${ADAPTER_DIR}/TruSeq3-SE.fa:2:30:10:2:TRUE"
        STAR_INPUT="${TRIM_R1}"
    fi

    echo "[${SAMPLE}] Mode: ${MODE}. Job Started..."

    # --- 2. Trimmomatic ---
    if [ ! -f "${TRIM_R1}" ]; then
        java -jar ${TRIMMOMATIC_JAR} ${TRIM_CMD} \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
            >> ${TRIMMED_DIR}/logs/${SAMPLE}.trim.log 2>&1
    fi

    # --- 3. STAR Alignment & 即时清理 ---
    if [ ! -f "${BAM_FILE}" ]; then
        mkdir -p ${SAMPLE_MAPPED_DIR}
        STAR --runThreadN ${THREADS_PER_JOB} \
            --genomeDir ${STAR_INDEX_DIR} --sjdbGTFfile ${GTF_FILE} \
            --readFilesIn ${STAR_INPUT} --readFilesCommand zcat \
            --outFileNamePrefix ${SAMPLE_MAPPED_DIR}/${SAMPLE}_ \
            --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
            --outFilterMismatchNmax 3 --quantMode TranscriptomeSAM \
            --limitBAMsortRAM 2000000000 \
            >> ${MAPPED_DIR}/logs/${SAMPLE}.star.log 2>&1
        
        # 【核心修改：即时清理】
        # STAR 一旦成功生成 BAM，立刻删除巨大的 clean fastq 文件
        echo "[${SAMPLE}] Cleaning up intermediate fastq files..."
        if [ "$MODE" == "PE" ]; then
            rm -f ${TRIM_R1} ${TRIM_R2} ${TRIMMED_DIR}/${SAMPLE}_*.unpaired.fq.gz
        else
            rm -f ${TRIM_R1}
        fi
    fi

    # --- 4. StringTie ---
    if [ ! -f "${STRINGTIE_ABUND}" ]; then
        stringtie -e -B -p ${THREADS_PER_JOB} -G ${GTF_FILE} \
            -o ${STRINGTIE_DIR}/${SAMPLE}.stringtie.gtf \
            -A ${STRINGTIE_ABUND} ${BAM_FILE} \
            >> ${STRINGTIE_DIR}/logs/${SAMPLE}.stringtie.log 2>&1
    fi

    echo "[${SAMPLE}] Finished."
}
export -f process_sample_parallel
export PROJECT_DIR RAW_DATA_DIR TRIMMED_DIR MAPPED_DIR STRINGTIE_DIR
export TRIMMOMATIC_JAR ADAPTER_DIR STAR_INDEX_DIR GTF_FILE THREADS_PER_JOB

# ===================== 4. 并行执行 =====================

echo ">>> Starting Loop..."
RUNNING_JOBS=0
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample_parallel "${SAMPLE}" &
    RUNNING_JOBS=$((RUNNING_JOBS + 1))
    if [ ${RUNNING_JOBS} -ge ${MAX_PARALLEL_SAMPLES} ]; then
        wait -n || true
        RUNNING_JOBS=$((RUNNING_JOBS - 1))
    fi
done
wait || echo "Warning: Some jobs failed."

# ===================== 5. 汇总与最终清理 =====================

echo ">>> Generating Final Matrices..."
ALL_BAMS=$(find ${MAPPED_DIR} -name "*_Aligned.sortedByCoord.out.bam" | sort)

if [ -n "$ALL_BAMS" ]; then
    # 智能 FeatureCounts (自动探测 PE/SE)
    if ls ${RAW_DATA_DIR}/*_2.fastq.gz 1> /dev/null 2>&1; then FC_OPTS="-p"; else FC_OPTS=""; fi
    
    featureCounts -T ${THREADS_FC} ${FC_OPTS} -t exon -g gene_id \
        -a ${GTF_FILE} -o ${COUNTS_DIR}/all_samples_raw_counts.txt \
        ${ALL_BAMS} > ${COUNTS_DIR}/featureCounts.log 2>&1
        
    cut -f1,7- "${COUNTS_DIR}/all_samples_raw_counts.txt" | \
        sed 's|'${MAPPED_DIR}'/[^/]*/||g' | sed 's|_Aligned.sortedByCoord.out.bam||g' | \
        grep -v "^#" > "${COUNTS_DIR}/all_samples_counts_matrix.txt"
fi

# 合并 TPM (Pythonic Awk)
TMP_MERGE="${STRINGTIE_DIR}/tmp_merge"
mkdir -p ${TMP_MERGE}
FIRST_SAMPLE="${SAMPLES[0]}"
awk -F'\t' 'NR>1 {print $1"\t"$2}' "${STRINGTIE_DIR}/${FIRST_SAMPLE}.gene_abundance.txt" > ${TMP_MERGE}/gene_info.txt
HEADER="Gene_ID\tGene_Name"
for SAMPLE in "${SAMPLES[@]}"; do
    if [ -f "${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt" ]; then
        HEADER="${HEADER}\t${SAMPLE}"
        awk -F'\t' 'NR>1 {print $9}' "${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt" > "${TMP_MERGE}/${SAMPLE}.tpm"
        awk -F'\t' 'NR>1 {print $8}' "${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt" > "${TMP_MERGE}/${SAMPLE}.fpkm"
    fi
done
echo -e "${HEADER}" > "${STRINGTIE_DIR}/all_samples_TPM_matrix.txt"
paste ${TMP_MERGE}/gene_info.txt ${TMP_MERGE}/*.tpm >> "${STRINGTIE_DIR}/all_samples_TPM_matrix.txt"
echo -e "${HEADER}" > "${STRINGTIE_DIR}/all_samples_FPKM_matrix.txt"
paste ${TMP_MERGE}/gene_info.txt ${TMP_MERGE}/*.fpkm >> "${STRINGTIE_DIR}/all_samples_FPKM_matrix.txt"
rm -rf ${TMP_MERGE}

# ===================== 6. 最终大扫除 (Step 6) =====================

echo ">>>> 开始最终清理..."

# 1. 确保删除所有剩余的 clean fastq (以防万一函数里没删干净)
if [ -d "${TRIMMED_DIR}" ]; then
    echo "清理残留的 Fastq 文件..."
    # 删除所有 .fq.gz 文件
    rm -f ${TRIMMED_DIR}/*.fq.gz
fi

# 2. 删除 StringTie 的中间文件 (.ctab)
# StringTie 会生成一堆 t_data.ctab, e_data.ctab 等，我们只需要 .gtf 和 .abundance
if [ -d "${STRINGTIE_DIR}" ]; then
    echo "清理 StringTie ctab 文件..."
    find ${STRINGTIE_DIR} -name "*.ctab" -type f -delete
fi

echo "======================================================="
echo "Pipeline Completed Successfully!"
echo "Matrices location: ${COUNTS_DIR} and ${STRINGTIE_DIR}"
echo "======================================================="
