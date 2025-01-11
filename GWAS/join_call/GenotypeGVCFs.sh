#!/bin/bash
#CSUB -J GWAS       # job name
#CSUB -q c01             # queue name
#CSUB -o %J.out          # output file
#CSUB -e %J.error        # error file
#CSUB -n 16              # 使用16个CPU核
#CSUB -cwd /share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/vcf

Ref_genome="/share/org/YZWL/yzwl_zhengbin/project/00-ref/huxu.fa"

for chr in {1..38} Z W MT; do 
    chr_script="chr${chr}_genotype.sh"
    echo '#!/bin/bash' >> "$chr_script"
    echo '#CSUB -J 'gwas_chr${chr}'' >> "$chr_script"  # 作业名
    echo '#CSUB -q c01' >> "$chr_script"  # 指定队列
    echo '#CSUB -o 'chr${chr}_genotype.out'' >> "$chr_script"  # 输出文件
    echo '#CSUB -e 'chr${chr}_genotype.err'' >> "$chr_script"  # 错误文件
    echo '#CSUB -n 32' >> "$chr_script"  # 使用16个线程
    echo '#CSUB -R span[hosts=1]' >> "$chr_script"
    echo '#CSUB -cwd /share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/vcf' >> "$chr_script"  # 指定作业运行目录
    echo '' >> "$chr_script"
    echo 'source /share/org/YZWL/yzwl_zhengbin/software/miniconda3/bin/activate' >> "$chr_script"
    echo 'conda activate gwas' >> "$chr_script"

    echo 'gatk --java-options "-Xmx400g" GenotypeGVCFs \' >> "$chr_script"
    echo '-R' ${Ref_genome} '\' >> "$chr_script"
    echo '-V 'chr${chr}.g.vcf.gz' \' >> "$chr_script"
    echo '-O 'chr${chr}.vcf.gz'' >> "$chr_script"
done
