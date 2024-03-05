#!/bin/bash

# 设置参考基因组文件名
reference_genome="/home/binz/MDV/ref/MDV_ref.fna"

# 设置基因组目录
genome_directory="/home/binz/MDV/genome/fasta/"

# 遍历基因组目录中的所有文件
for genome_file in ${genome_directory}*.fna; do
    # 获取基因组文件名（不含路径和扩展名）
    genome_name=$(basename "${genome_file}" .fna)

    # 运行 NUCmer 比对基因组到参考基因组，并使用基因组名称作为前缀
    nucmer -p "${genome_name}" -t 24 "${reference_genome}" "${genome_file}"

    # 使用delta-filter 过滤比对结果
    delta-filter -1 "${genome_name}.delta" > "${genome_name}.filtered.delta"

    # 提取比对结果的摘要
    show-coords -T -H -r -d "${genome_name}.filtered.delta" > "${genome_name}.coords"

    # 提取snp和InDel
    show-snps -C -l -r -x 1 -T "${genome_name}.filtered.delta" > "${genome_name}.filter.snps"

    # 将snps文件转换为vcf文件
    /home/binz/MDV/result/call_SNP_INDEL_CNV/MUMmerSNPs2VCF.py "${genome_name}.filter.snps" "${genome_name}.filter.snps.vcf"
done
