#!/bin/bash

# 激活指定的 Conda 环境
source /share/home/yzwl_zhengbin/miniconda3/bin/activate manta

input_bam=$1

ref=/share/home/yzwl_zhengbin/ref/huxu_v23_ref/chicken.v23.fa

input_bam_basename=$(basename $input_bam)

output_prefix=${input_bam_basename%.bam}

output_dir=/share/home/yzwl_zhengbin/GWAS/call_sv/03_vcf/manta/${output_prefix}_manta
mkdir $output_dir

samtools sort -@8 $input_bam -o ${output_dir}/${output_prefix}_sorted.bam

samtools index ${output_dir}/${output_prefix}_sorted.bam

python2 /share/home/yzwl_zhengbin/miniconda3/envs/manta/share/manta-1.6.0-2/bin/configManta.py \
        --bam ${output_dir}/${output_prefix}_sorted.bam \
        --runDir ${output_dir} \
        --referenceFasta /share/home/yzwl_zhengbin/ref/huxu_v23_ref/chicken.v23.fa

python2 ${output_dir}/runWorkflow.py -j 4 -m local

mv ${output_dir}/results/variants/diploidSV.vcf.gz ${output_dir}/results/variants/${output_prefix}_diploidSV.vcf.gz
mv ${output_dir}/results/variants/${output_prefix}_diploidSV.vcf.gz ${output_dir}

source /share/home/yzwl_zhengbin/miniconda3/bin/activate py3

vcftools --gzvcf  ${output_dir}/${output_prefix}_diploidSV.vcf.gz --minQ 20 --recode --recode-INFO-all --out ${output_dir}/${output_prefix}_clean

bgzip ${output_dir}/${output_prefix}_clean.recode.vcf

tabix -p vcf ${output_dir}/${output_prefix}_clean.recode.vcf.gz

rm -r ${output_dir}/${output_prefix}_sorted.bam* \
        ${output_dir}/runWorkflow.py* \
        ${output_dir}/work*
