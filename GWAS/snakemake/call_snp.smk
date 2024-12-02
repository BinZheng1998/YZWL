# 设置全局变量
FASTQ_DIR = "/home/bzheng/others_species_analysis/cattle/fastq/Bos_primigenius33"
REFERENCE_GENOME = "/home/bzheng/others_species_analysis/cattle/ref/Btau_5.0.1_genome.fa"
PICARD_JAR = "/home/bzheng/software/picard.jar"
GATK = "/home/bzheng/software/gatk-4.5.0.0/gatk"
JAVA = "/home/bzheng/software/jdk-22.0.1/bin/java"
THREADS = 16  # 线程数，视情况调整

# 从txt文件读取样本前缀名
SAMPLES_FILE = "/home/bzheng/others_species_analysis/cattle/vcf/Bos_primigenius33/sample_name.txt"  
with open(SAMPLES_FILE, 'r') as f:
    FASTQ_FILES_R1 = [line.strip() for line in f if line.strip()]

# 默认目标
rule all:
    input:
        expand("output/{sample}.raw.snps.indels.g.vcf", sample=FASTQ_FILES_R1),
        expand("output/{sample}_dups_marked.bam", sample=FASTQ_FILES_R1),
        expand("output/{sample}_dups_marked.bai", sample=FASTQ_FILES_R1)

# 规则：Quality Control with fastp
rule fastp:
    input:
        r1 = lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}_1.fastq.gz",
        r2 = lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}_2.fastq.gz"
    output:
        r1_clean = "output/{sample}_R1_clean.fastq.gz",
        r2_clean = "output/{sample}_R2_clean.fastq.gz"
    threads: THREADS
    shell:
        "fastp -i {input.r1} -I {input.r2} "
        "-q 20 -u 50 -l 30 -e 20 "
        "-o {output.r1_clean} -O {output.r2_clean} "
        "--thread {threads}"

# 规则：单端FASTQ质量控制
rule fastp_single:
    input:
        r1 = lambda wildcards: f"{FASTQ_DIR}/{wildcards.sample}.fastq.gz"
    output:
        r1_clean = "output/{sample}_clean.fastq.gz"
    threads: THREADS
    shell:
        "fastp -i {input.r1} "
        "-q 20 -u 50 -l 30 -e 20 "
        "-o {output.r1_clean} "
        "--thread {threads}"

# 规则：BWA 比对并生成 BAM 文件（配对端）
rule bwa_align_paired:
    input:
        r1_clean = "output/{sample}_R1_clean.fastq.gz",
        r2_clean = "output/{sample}_R2_clean.fastq.gz"
    output:
        bam = "output/{sample}.bam",
        sai_r1 = "output/{sample}_R1.sai",
        sai_r2 = "output/{sample}_R2.sai"
    threads: THREADS
    shell:
        "bwa aln -t {threads} -l 1024 -n 0.01 {REFERENCE_GENOME} {input.r1_clean} > {output.sai_r1} && "
        "bwa aln -t {threads} -l 1024 -n 0.01 {REFERENCE_GENOME} {input.r2_clean} > {output.sai_r2} && "
        "bwa sampe {REFERENCE_GENOME} {output.sai_r1} {output.sai_r2} "
        "{input.r1_clean} {input.r2_clean} | "
        "samtools view -bS > {output.bam}"

# 规则：BWA 比对并生成 BAM 文件（单端）
rule bwa_align_single:
    input:
        r1_clean = "output/{sample}_clean.fastq.gz"
    output:
        bam = "output/{sample}.bam",
        sai_r1 = "output/{sample}.sai"
    threads: THREADS
    shell:
        "bwa aln -t {threads} -l 1024 -n 0.01 {REFERENCE_GENOME} {input.r1_clean} > {output.sai_r1} && "
        "bwa samse {REFERENCE_GENOME} {output.sai_r1} {input.r1_clean} | "
        "samtools view -bS > {output.bam}"

# 规则：排序并去重
rule sort_and_dedup:
    input:
        bam = "output/{sample}.bam"
    output:
        sorted_bam = "output/{sample}_sorted.bam",
        sorted_bai = "output/{sample}_sorted.bam.bai",
        dedup_bam = "output/{sample}_dups_marked.bam",
        dedup_bai = "output/{sample}_dups_marked.bai",
        metrics = "output/{sample}_dups_metrics.txt"
    threads: THREADS
    shell:
        "{JAVA} -Xmx260g -jar {PICARD_JAR} SortSam "
        "I={input.bam} O={output.sorted_bam} SO=coordinate && "
        "samtools index {output.sorted_bam} && "
        "{JAVA} -Xmx260g -jar {PICARD_JAR} MarkDuplicates "
        "I={output.sorted_bam} O={output.dedup_bam} M={output.metrics} REMOVE_DUPLICATES=true && "
        "samtools index {output.dedup_bam}"

# 规则：GATK 变异检测
rule haplotype_caller:
    input:
        bam = "output/{sample}_dups_marked.bam"
    output:
        gvcf = "output/{sample}.raw.snps.indels.g.vcf"
    shell:
        "{GATK} --java-options '-Xmx260G' HaplotypeCaller "
        "-R {REFERENCE_GENOME} -I {input.bam} -ERC GVCF "
        "-mbq 20 --output-mode 'EMIT_ALL_CONFIDENT_SITES' "
        "-O {output.gvcf}"
