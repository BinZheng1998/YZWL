#!/bin/bash
#CSUB -J call_sv1
#CSUB -q cpu
#CSUB -o bam2gvcf1.out
#CSUB -e bam2gvcf1.error
#CSUB -n 12
#CSUB -cwd /share/home/yzwl_zhengbin/GWAS/call_sv/delly

# 设置输入和输出目录
input_parent_dir="/share/home/yzwl_zhengbin/GWAS/call_snp_indel"
reference_genome="/share/home/yzwl_zhengbin/ref/huxu_v23_ref/chicken.v23.fa"
# 遍历每个样本文件夹
for sample_dir in ${input_parent_dir}/*_output/; do
    sample_name=$(basename "${sample_dir}" | sed 's/_output//')  # 提取样本名

    # 使用delly进行SV分析 (在Python 3环境中运行)
    source /share/home/yzwl_zhengbin/miniconda3/bin/activate py3
    delly call \
        -o ${sample_name}.delly.bcf \
        -g $reference_genome \
        ${sample_dir}/${sample_name}_dups_marked.bam
    #将bcf文件转换为vcf文件
    bcftools view ${sample_name}.delly.bcf > ${sample_name}.delly.vcf

    echo "样本 ${sample_name} 的SV分析完成！"
done
