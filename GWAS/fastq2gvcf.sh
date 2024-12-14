#!/bin/bash

#source /share/org/YZWL/yzwl_zhengbin/software/miniconda3/bin/activate gwas
# 输入文件路径
input_file_1="$1"
input_file_2="$2"

# 检查输入参数是否为空
if [ -z "$input_file_1" ] || [ -z "$input_file_2" ]; then
    echo "Error: Please provide two input files."
    exit 1
fi

# 检查输入文件是否存在
if [ ! -f "$input_file_1" ] || [ ! -f "$input_file_2" ]; then
    echo "Error: Input file does not exist."
    exit 1
fi

# 参考基因组路径
reference_genome="/share/org/YZWL/yzwl_zhengbin/project/00-ref/huxu.fa"
picard_jar="/share/org/YZWL/yzwl_zhengbin/software/picard.jar"
gatk_jar="/share/org/YZWL/yzwl_zhengbin/software/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"
gatk="/share/org/YZWL/yzwl_zhengbin/software/gatk-4.5.0.0/gatk"
# 提取样本名称
sample=$(basename "$input_file_1" | sed 's/_[12].fq.gz//')
output_dir="/share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/${sample}_output"
mkdir -p $output_dir

fastp -i "$input_file_1" -I "$input_file_2" -q 20 -u 50 -l 50 -e 20 -o "${output_dir}/${sample}_clean_R1.fastq.gz" -O "${output_dir}/${sample}_clean_R2.fastq.gz"

# Step 2: 从FASTQ数据比对到BAM文件
bwa mem -t 24 -aM "$reference_genome" "${output_dir}/${sample}_clean_R1.fastq.gz" "${output_dir}/${sample}_clean_R2.fastq.gz" > "${output_dir}/${sample}_aligned.sam"

# Step 3: 将SAM文件转换为BAM文件
samtools view -bS "${output_dir}/${sample}_aligned.sam" > "${output_dir}/${sample}_aligned.bam"

# Step 4: 对BAM文件进行排序
samtools sort -@ 12 -m 8G -T "${output_dir}/${sample}_sorted" -o "${output_dir}/${sample}_sorted.bam" "${output_dir}/${sample}_aligned.bam"

# Step 5: 为BAM文件建立索引
samtools index "${output_dir}/${sample}_sorted.bam"

# Step 6: 添加或替换读组信息
java -Xmx240g -jar "$picard_jar" AddOrReplaceReadGroups \
      -I "${output_dir}/${sample}_sorted.bam" \
      -O "${output_dir}/${sample}_sorted_rg.bam" \
      --RGLB library \
      --RGPL illumina \
      --RGPU "$sample" \
      --RGSM "$sample"

# Step 7: 标记重复序列
java -Xmx240g -jar "$picard_jar" MarkDuplicates \
      -I "${output_dir}/${sample}_sorted_rg.bam" \
      -O "${output_dir}/${sample}_dups_marked.bam" \
      --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
      --REMOVE_DUPLICATES false \
      -M "${output_dir}/${sample}_dups_metrics.txt"

# Step 8: 为BAM文件建立索引
java -Xmx240g -jar "$picard_jar" BuildBamIndex \
      -I "${output_dir}/${sample}_dups_marked.bam" \
      -O "${output_dir}/${sample}_dups_marked.bai"

# Step 9: 使用GATK进行变异检测
"$gatk" --java-options "-Xmx240G" HaplotypeCaller \
     -R "$reference_genome" \
     -I "${output_dir}/${sample}_dups_marked.bam" \
     -ERC GVCF \
     -mbq 20 \
     --output-mode "EMIT_ALL_CONFIDENT_SITES" \
     -O "${output_dir}/${sample}.raw.snps.indels.g.vcf"

# 删除中间文件
rm -f "${output_dir}/${sample}_clean_1.fq.gz" \
      "${output_dir}/${sample}_clean_2.fq.gz" \
      "${output_dir}/${sample}_aligned.sam" \
      #"${output_dir}/${sample}_aligned.bam" \
      "${output_dir}/${sample}_sorted.bam" \
      "${output_dir}/${sample}_sorted.bam.bai" \
      "${output_dir}/${sample}_sorted_rg.bam" \
      "${output_dir}/${sample}_dups_metrics.txt" 
      #"${output_dir}/${sample}_dups_marked.bam" \
      #"${output_dir}/${sample}_dups_marked.bai"
