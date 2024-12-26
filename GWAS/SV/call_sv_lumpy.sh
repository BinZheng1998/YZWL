#!/bin/bash

# 定义输入参数
fastq_1="$1"
fastq_2="$2"
reference_genome="/share/home/yzwl_zhengbin/ref/huxu_v23_ref/chicken.v23.fa"

# 创建输出目录
sample=$(basename "${fastq_1}" | cut -d '_' -f 1)
output_dir="/share/home/yzwl_zhengbin/GWAS/call_sv/03_vcf/lumpy/${sample}_output"
mkdir -p $output_dir

# 激活 Conda 环境
source /share/home/yzwl_zhengbin/miniconda3/bin/activate py2

# 运行 BWA
bwa mem -t 12 -R "@RG\tID:${sample}\tSM:${sample}\tLB:lib\tPL:illumina" "${reference_genome}" "$fastq_1" "$fastq_2" \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -bS - \
    > "${output_dir}/${sample}.bam"

# 提取不匹配读取
samtools view -@ 8 -b -F 1294 "${output_dir}/${sample}.bam" > "${output_dir}/${sample}_unsorted.discordants.bam"

# 提取分割读取
samtools view -@ 8 -h "${output_dir}/${sample}.bam" \
    | /share/home/yzwl_zhengbin/miniconda3/envs/py2/bin/extractSplitReads_BwaMem -i stdin \
    | samtools view -bS - \
    > "${output_dir}/${sample}_unsorted.splitters.bam"

# 排序 BAM 文件
samtools sort -@ 8 "${output_dir}/${sample}_unsorted.discordants.bam" -o "${output_dir}/${sample}.discordants.bam"
samtools sort -@ 8 "${output_dir}/${sample}_unsorted.splitters.bam" -o "${output_dir}/${sample}.splitters.bam"
# 运行 Lumpy
lumpyexpress \
    -B "${output_dir}/${sample}.bam" \
    -S "${output_dir}/${sample}.splitters.bam" \
    -D "${output_dir}/${sample}.discordants.bam" \
    -o "${output_dir}/${sample}.lumpy.vcf"

# 删除中间文件
rm -rf "${output_dir}/${sample}.discordants.bam"
rm -rf "${output_dir}/${sample}_unsorted.discordants.bam"
rm -rf "${output_dir}/${sample}.splitters.bam"
rm -rf "${output_dir}/${sample}_unsorted.splitters.bam"
