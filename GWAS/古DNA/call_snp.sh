#/bin/bash

# 输入FASTQ文件夹路径
fastq_dir="/home/bzheng/others_species_analysis/cattle/fastq/Bos_primigenius33"
# 参考基因组路径(这里用的是胡须鸡的参考基因组)
reference_genome="/home/bzheng/others_species_analysis/cattle/ref/Btau_5.0.1_genome.fa"
# Picard.jar 路径
picard_jar="/home/bzheng/software/picard.jar"
gatk="/home/bzheng/software/gatk-4.5.0.0/gatk"
java="/home/bzheng/software/jdk-22.0.1/bin/java"

fastq_files=("$fastq_dir"/*.fastq.gz)

# 处理每个样本
for fastq_file in "${fastq_files[@]}"; do
    # 提取样本名称（假设文件名格式为sample_name_R1.fastq.gz或sample_name.fastq.gz）
    sample=$(basename "$fastq_file" | sed 's/_[12].fastq.gz//;s/.fastq.gz//')
    echo "Processing sample: $sample"

    output_dir="${sample}_output"
    mkdir -p "$output_dir"

    r2_file="${fastq_dir}/${sample}_2.fastq.gz"
    if [[ -f "$r2_file" ]]; then
        echo "Detected paired-end data for $sample"

        # Step 1: 使用fastp进行质量控制和去除低质量序列
        echo "Step 1: Running fastp for quality control..."
        fastp -i "$fastq_file" -I "$r2_file" \
              -q 20 -u 50 -l 30  -e 20 \
            -o "${output_dir}/${sample}_R1_clean.fastq.gz" \
            -O "${output_dir}/${sample}_R2_clean.fastq.gz" \
            --thread 12

        # Step 2: 从清理后的压缩FASTQ数据比对到排序的BAM文件
        echo "Step 2: Aligning reads and sorting BAM..."
        bwa aln -t 16 -l 1024 -n 0.01 "$reference_genome" "${output_dir}/${sample}_R1_clean.fastq.gz" > "${output_dir}/${sample}_R1.sai"
        bwa aln -t 16 -l 1024 -n 0.01 "$reference_genome" "${output_dir}/${sample}_R2_clean.fastq.gz" > "${output_dir}/${sample}_R2.sai"

# 使用 bwa sampe 合并两个sai文件并生成SAM文件
        bwa sampe "$reference_genome" \
                "${output_dir}/${sample}_R1.sai" "${output_dir}/${sample}_R2.sai" \
                "${output_dir}/${sample}_R1_clean.fastq.gz" "${output_dir}/${sample}_R2_clean.fastq.gz" | \
        samtools view -bS > "${output_dir}/${sample}.bam" 
        samtools index "${output_dir}/${sample}.bam"
        samtools sort -@ 16 -m 10G "${output_dir}/${sample}.bam" -T "${output_dir}/${sample}_sorted" -o "${output_dir}/${sample}_sorted.bam"
    else
        echo "Detected single-end data for $sample"
        # 单端测序的处理流程

        # Step 1: 使用fastp进行质量控制和去除低质量序列
        echo "Step 1: Running fastp for quality control..."
        fastp -i "$fastq_file" \
              -q 20 -u 50 -l 30  -e 20 \
            -o "${output_dir}/${sample}_clean.fastq.gz" \
            --thread 16

        # Step 2: 从清理后的压缩FASTQ数据比对到排序的BAM文件
        echo "Step 2: Aligning reads and sorting BAM..."
        bwa aln -t 16 -l 1024 -n 0.01 "$reference_genome" "${output_dir}/${sample}_clean.fastq.gz" > "${output_dir}/${sample}.sai"
        bwa samse "$reference_genome" "${output_dir}/${sample}.sai" "${output_dir}/${sample}_clean.fastq.gz" | \
        samtools view -bS > "${output_dir}/${sample}.bam" 
        samtools index "${output_dir}/${sample}.bam"
        samtools sort -@ 16 -m 10G "${output_dir}/${sample}.bam" -T "${output_dir}/${sample}_sorted" -o "${output_dir}/${sample}_sorted.bam"
    fi

    # Step 3: 为BAM文件建立索引
    echo "Step 3: Building BAM index..."
    samtools index "${output_dir}/${sample}_sorted.bam"
    samtools view -bF 4 "${output_dir}/${sample}_sorted.bam" > "${output_dir}/${sample}_filter.bam"
    samtools index "${output_dir}/${sample}_filter.bam"

    # Step 4: 添加或替换读组信息
    echo "Step 4: Adding or replacing read groups..."
    $java -Xmx260g -jar "$picard_jar" AddOrReplaceReadGroups \
          -I "${output_dir}/${sample}_filter.bam" \
          -O "${output_dir}/${sample}_sorted_rg.bam" \
          --RGLB library \
          --RGPL illumina \
          --RGPU $sample \
          --RGSM $sample

    # Step 5: 标记重复序列
    echo "Step 5: Marking duplicates..."
    $java -Xmx260g -jar "$picard_jar" MarkDuplicates \
          -I "${output_dir}/${sample}_sorted_rg.bam" \
          -O "${output_dir}/${sample}_dups_marked.bam" \
          --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
          --REMOVE_DUPLICATES true \
          -M "${output_dir}/${sample}_dups_metrics.txt"

    # Step 6: 为BAM文件建立索引
    echo "Step 6: Building BAM index..."
    $java -Xmx260g -jar "$picard_jar" BuildBamIndex \
          -I "${output_dir}/${sample}_dups_marked.bam" \
          -O "${output_dir}/${sample}_dups_marked.bai"

    # Step 7: 使用GATK进行变异检测
    echo "Step 7: Variant calling with GATK..."
    $gatk --java-options "-Xmx260G" HaplotypeCaller \
         -R "$reference_genome" \
         -I "${output_dir}/${sample}_dups_marked.bam" \
         -ERC GVCF \
         -mbq 20 \
         --output-mode "EMIT_ALL_CONFIDENT_SITES" \
         -O "${output_dir}/${sample}.raw.snps.indels.g.vcf"

    # 删除中间文件（如需要，可取消注释）
    rm -f "${output_dir}/${sample}_R1_clean.fastq.gz" \
          "${output_dir}/${sample}_R2_clean.fastq.gz" \
          "${output_dir}/${sample}_clean.fastq.gz" \
          "${output_dir}/${sample}.sai"
          "${output_dir}/${sample}_sorted.bam" \
          "${output_dir}/${sample}_sorted.bam.bai" \
          "${output_dir}/${sample}_filter.bam"
          "${output_dir}/${sample}_filter.bam.bai"
          "${output_dir}/${sample}_sorted_rg.bam" \
          "${output_dir}/${sample}_dups_marked.bam" \
          "${output_dir}/${sample}_dups_marked.bam.bai"
done
