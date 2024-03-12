#!/bin/bash
#CSUB -J zhengbin
#CSUB -q cpu
#CSUB -o fastq2gvcf.out
#CSUB -e fastq2gvcf.error
#CSUB -n 12
#CSUB -cwd /share/home/yzwl_zhengbin/GWAS/call_snp_indel

source /share/home/yzwl_zhengbin/miniconda3/bin/activate py3

# 记录开始时间
start_time=$(date +%s)

# 创建用于存储时间数据的文件
time_file="time_data.txt"
echo "Sample,Step,Time (seconds)" > "$time_file"

# 输入FASTQ文件夹路径
fastq_dir="/home/binz/GWAS_test/test_fastq"
# 参考基因组路径(这里用的是胡须鸡的参考基因组)
reference_genome="/home/binz/ref/chicken/chicken.v23.fa"
# Picard.jar 路径
picard_jar="/home/binz/software/picard.jar"
# GATK.jar 路径
gatk_jar="/home/binz/software/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar"

# 获取输入FASTQ文件夹中的所有FASTQ文件
fastq_files=("$fastq_dir"/*_R1.fastq.gz)

# 处理每个样本
for fastq_file in "${fastq_files[@]}"; do
    # 提取样本名称（假设文件名格式为sample_name_R1.fastq.gz）
    sample=$(basename "$fastq_file" | sed 's/_R[12].fastq.gz//')
    echo "Processing sample: $sample"

    # 创建输出文件夹
    output_dir="${sample}_output"
    mkdir -p "$output_dir"

    # 提取配对的R2文件名
    r2_file="${fastq_dir}/${sample}_R2.fastq.gz"

    # Step 1: 解压缩FASTQ文件
    echo "Step 1: Decompressing FASTQ files..."
    gunzip -c "$fastq_file" > "${output_dir}/${sample}_R1.fastq"
    gunzip -c "$r2_file" > "${output_dir}/${sample}_R2.fastq"

# Step 2: 从FASTQ数据比对到BAM文件
    echo "Step 2: Aligning reads using BWA mem..."
    bwa_mem_start_time=$(date +%s)
    bwa mem -aM  "$reference_genome" "${output_dir}/${sample}_R1.fastq" "${output_dir}/${sample}_R2.fastq" > "${output_dir}/${sample}_aligned.sam"
    bwa_mem_end_time=$(date +%s)
    bwa_mem_time=$((bwa_mem_end_time - bwa_mem_start_time))
    echo "$sample,bwa mem,$bwa_mem_time" >> "$time_file"

    # Step 3: 将SAM文件转换为BAM文件
    echo "Step 3: Converting SAM to BAM..."
    samtools_view_start_time=$(date +%s)
    samtools view -bS "${output_dir}/${sample}_aligned.sam" > "${output_dir}/${sample}_aligned.bam"
    samtools_view_end_time=$(date +%s)
    samtools_view_time=$((samtools_view_end_time - samtools_view_start_time))
    echo "$sample,samtools view,$samtools_view_time" >> "$time_file"

    # Step 4: 对BAM文件进行排序
    echo "Step 4: Sorting BAM..."
    samtools_sort_start_time=$(date +%s)
    samtools sort "${output_dir}/${sample}_aligned.bam" -f "${output_dir}/${sample}_sorted.bam"
    samtools_sort_end_time=$(date +%s)
    samtools_sort_time=$((samtools_sort_end_time - samtools_sort_start_time))
    echo "$sample,samtools sort,$samtools_sort_time" >> "$time_file"

    # Step 5: 为BAM文件建立索引
    echo "Step 5: Building BAM index..."
        samtools_index_start_time=$(date +%s)
    samtools index "${output_dir}/${sample}_sorted.bam"
    samtools_index_end_time=$(date +%s)
    samtools_index_time=$((samtools_index_end_time - samtools_index_start_time))
    echo "$sample,samtools index,$samtools_index_time" >> "$time_file"

    # Step 6: 添加或替换读组信息
    echo "Step 6: Adding or replacing read groups..."
    java -Xmx60g -jar "$picard_jar" AddOrReplaceReadGroups \
          -I "${output_dir}/${sample}_sorted.bam" \
          -O "${output_dir}/${sample}_sorted_rg.bam" \
          --RGLB library \
          --RGPL illumina \
          --RGPU $sample \
          --RGSM $sample
    picard_addrg_time=$((end_time - start_time))
    echo "$sample,picard AddOrReplaceReadGroups,$picard_addrg_time" >> "$time_file"

    # Step 7: 标记重复序列
    echo "Step 7: Marking duplicates..."
    java -Xmx60g -jar "$picard_jar" MarkDuplicates \
          -I "${output_dir}/${sample}_sorted_rg.bam" \
          -O "${output_dir}/${sample}_dups_marked.bam" \
          --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
          --REMOVE_DUPLICATES false \
          -M "${output_dir}/${sample}_dups_metrics.txt"
    picard_markdup_time=$((end_time - start_time))
    echo "$sample,picard MarkDuplicates,$picard_markdup_time" >> "$time_file"

    # Step 8: 为BAM文件建立索引
    echo "Step 8: Building BAM index..."
    java -Xmx60g -jar "$picard_jar" BuildBamIndex \
          -I "${output_dir}/${sample}_dups_marked.bam" \
          -O "${output_dir}/${sample}_dups_marked.bai"
    picard_buildbamindex_time=$((end_time - start_time))
    echo "$sample,picard BuildBamIndex,$picard_buildbamindex_time" >> "$time_file"

    # Step 9: 使用GATK进行变异检测
    echo "Step 9: Variant calling with GATK..."
    gatk_haplotypecaller_start_time=$(date +%s)
    gatk --java-options "-Xmx60G" HaplotypeCaller \
         -R "$reference_genome" \
         -I "${output_dir}/${sample}_dups_marked.bam" \
         -ERC GVCF \
         -mbq 20 \
         --output-mode "EMIT_ALL_CONFIDENT_SITES" \
         -O "${output_dir}/${sample}.raw.snps.indels.g.vcf"
    gatk_haplotypecaller_time=$((end_time - start_time))
    echo "$sample,GATK HaplotypeCaller,$gatk_haplotypecaller_time" >> "$time_file"

    # 删除中间文件
    #rm -f "${output_dir}/${sample}_R1.fastq" \
          #"${output_dir}/${sample}_R2.fastq" \
          #"${output_dir}/${sample}_aligned.sam" \
          #"${output_dir}/${sample}_aligned.bam" \
          #"${output_dir}/${sample}_sorted.bam" \
          #"${output_dir}/${sample}_sorted.bam.bai" \
          #"${output_dir}/${sample}_sorted_rg.bam" \
          #"${output_dir}/${sample}_deduplicated.bam.bai" 
done

# 记录结束时间
end_time=$(date +%s)

# 计算总时间
total_time=$((end_time - start_time))
echo "Total processing time: $total_time seconds" >> "$time_file"

  
