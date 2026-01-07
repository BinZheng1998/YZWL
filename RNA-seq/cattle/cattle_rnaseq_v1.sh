#!/bin/bash
# =============================================================================
# RNA-seq Analysis Pipeline (Cattle Data) - Final Integrated Version
# 功能: 自动处理 SE/PE 数据，支持 .fastq.gz 和 .fq.gz，包含即时清理
# =============================================================================

set -e

# ===================== 1. 核心配置 =====================

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DATA_DIR="${PROJECT_DIR}"
RESULTS_DIR="${PROJECT_DIR}/results"

# --- 参考基因组 ---
REF_DIR="/home/bzheng/project/04_rnaseq/cattle/0.ref/genome_reference"
REF_GENOME_FA="${REF_DIR}/GCF_002263795.1_ARS-UCD1.2_genomic.fna"
GTF_FILE="${REF_DIR}/GCF_002263795.1_ARS-UCD1.2_genomic.gtf"
STAR_INDEX_DIR="/home/bzheng/project/04_rnaseq/cattle/0.ref/star_index"

# --- 软件路径 ---
TRIMMOMATIC_JAR="/home/bzheng/.conda/envs/cattle/share/trimmomatic-0.39-2/trimmomatic.jar"
ADAPTER_DIR="/home/bzheng/.conda/envs/cattle/share/trimmomatic-0.39-2/adapters"

# --- 并行策略 ---
MAX_PARALLEL_SAMPLES=4  # 保持为4以防 STAR 内存溢出
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
echo ">>> Scanning for samples..."

# 查找所有 .fastq.gz 和 .fq.gz 文件，去除 _1/_2 后缀
SAMPLES=($(find . -maxdepth 1 -name "*.fastq.gz" -o -name "*.fq.gz" | \
    sed 's|^\./||' | \
    sed -E 's/(_1|_2)?\.(fastq|fq)\.gz//' | \
    sort -u))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "Error: No .fastq.gz or .fq.gz files found."
    exit 1
fi

echo ">>> Found ${#SAMPLES[@]} samples."

# ===================== 3. 处理函数 (含即时清理) =====================

process_sample_parallel() {
    local SAMPLE=$1
    local SAMPLE_MAPPED_DIR="${MAPPED_DIR}/${SAMPLE}"
    local BAM_FILE="${SAMPLE_MAPPED_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    local STRINGTIE_ABUND="${STRINGTIE_DIR}/${SAMPLE}.gene_abundance.txt"

    # --- 1. 自动判断 SE/PE 及 后缀 (修正版) ---
    MODE=""
    R1=""
    R2=""
    SUFFIX=""

    # 优先检测 PE (R2存在即为PE)
    if [ -f "${RAW_DATA_DIR}/${SAMPLE}_2.fastq.gz" ]; then
        MODE="PE"
        SUFFIX=".fastq.gz"
    elif [ -f "${RAW_DATA_DIR}/${SAMPLE}_2.fq.gz" ]; then
        MODE="PE"
        SUFFIX=".fq.gz"
    else
        MODE="SE"
    fi

    # 根据模式设置路径和命令
    if [ "$MODE" == "PE" ]; then
        R1="${RAW_DATA_DIR}/${SAMPLE}_1${SUFFIX}"
        R2="${RAW_DATA_DIR}/${SAMPLE}_2${SUFFIX}"
        
        TRIM_R1="${TRIMMED_DIR}/${SAMPLE}_1.clean.fq.gz"
        TRIM_R2="${TRIMMED_DIR}/${SAMPLE}_2.clean.fq.gz"
        
        TRIM_CMD="PE -phred33 ${R1} ${R2} ${TRIM_R1} ${TRIMMED_DIR}/${SAMPLE}_1.unpaired.fq.gz ${TRIM_R2} ${TRIMMED_DIR}/${SAMPLE}_2.unpaired.fq.gz ILLUMINACLIP:${ADAPTER_DIR}/TruSeq3-PE-2.fa:2:30:10:2:TRUE"
        STAR_INPUT="${TRIM_R1} ${TRIM_R2}"
    else
        # SE 模式：探测 R1 具体名称
        if [ -f "${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz" ]; then
            R1="${RAW_DATA_DIR}/${SAMPLE}_1.fastq.gz"
        elif [ -f "${RAW_DATA_DIR}/${SAMPLE}_1.fq.gz" ]; then
            R1="${RAW_DATA_DIR}/${SAMPLE}_1.fq.gz"
        elif [ -f "${RAW_DATA_DIR}/${SAMPLE}.fastq.gz" ]; then
            R1="${RAW_DATA_DIR}/${SAMPLE}.fastq.gz"
        elif [ -f "${RAW_DATA_DIR}/${SAMPLE}.fq.gz" ]; then
            R1="${RAW_DATA_DIR}/${SAMPLE}.fq.gz"
        else
            echo "Error: Cannot find input file for ${SAMPLE}"
            return 1
        fi

        TRIM_R1="${TRIMMED_DIR}/${SAMPLE}.clean.fq.gz"
        TRIM_CMD="SE -phred33 ${R1} ${TRIM_R1} ILLUMINACLIP:${ADAPTER_DIR}/TruSeq3-SE.fa:2:30:10:2:TRUE"
        STAR_INPUT="${TRIM_R1}"
    fi

    echo "[${SAMPLE}] Mode: ${MODE}. Input: $(basename ${R1})..."

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
            --limitBAMsortRAM 4000000000 --twopassMode Basic \
            >> ${MAPPED_DIR}/logs/${SAMPLE}.star.log 2>&1
        
        # 【即时清理】STAR 成功后删除 clean fastq
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
    # 智能探测 FeatureCounts 参数 (兼容 fastq.gz 和 fq.gz)
    if ls ${RAW_DATA_DIR}/*_2.fastq.gz 1> /dev/null 2>&1 || ls ${RAW_DATA_DIR}/*_2.fq.gz 1> /dev/null 2>&1; then 
        FC_OPTS="-p"
    else 
        FC_OPTS=""
    fi
    
    echo "FeatureCounts Opts: ${FC_OPTS}"
    featureCounts -T ${THREADS_FC} ${FC_OPTS} -t exon -g gene_id \
        -a ${GTF_FILE} -o ${COUNTS_DIR}/all_samples_raw_counts.txt \
        ${ALL_BAMS} > ${COUNTS_DIR}/featureCounts.log 2>&1
        
    cut -f1,7- "${COUNTS_DIR}/all_samples_raw_counts.txt" | \
        sed 's|'${MAPPED_DIR}'/[^/]*/||g' | sed 's|_Aligned.sortedByCoord.out.bam||g' | \
        grep -v "^#" > "${COUNTS_DIR}/all_samples_counts_matrix.txt"
fi

# 合并 TPM
TMP_MERGE="${STRINGTIE_DIR}/tmp_merge"
mkdir -p ${TMP_MERGE}
FIRST_SAMPLE="${SAMPLES[0]}"
# 获取基因信息头
if [ -f "${STRINGTIE_DIR}/${FIRST_SAMPLE}.gene_abundance.txt" ]; then
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
fi
rm -rf ${TMP_MERGE}

# ===================== 6. 最终大扫除 =====================

echo ">>>> 开始最终清理..."
if [ -d "${TRIMMED_DIR}" ]; then
    rm -f ${TRIMMED_DIR}/*.fq.gz
fi
if [ -d "${STRINGTIE_DIR}" ]; then
    find ${STRINGTIE_DIR} -name "*.ctab" -type f -delete
fi

echo "======================================================="
echo "Pipeline Completed Successfully!"
echo "Matrices location: ${COUNTS_DIR} and ${STRINGTIE_DIR}"
echo "======================================================="
