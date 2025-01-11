#!/bin/bash
#CSUB -J GWAS    
#CSUB -q c01           
#CSUB -o %J.out         
#CSUB -e %J.error       
#CSUB -n 16             
#CSUB -cwd vcf

Ref_genome="huxu.fa"

for chr in {1..38} Z W MT; do 
    chr_script="chr${chr}_genotype.sh"
    echo '#!/bin/bash' >> "$chr_script"
    echo '#CSUB -J 'gwas_chr${chr}'' >> "$chr_script" 
    echo '#CSUB -q c01' >> "$chr_script"  
    echo '#CSUB -o 'chr${chr}_genotype.out'' >> "$chr_script"  
    echo '#CSUB -e 'chr${chr}_genotype.err'' >> "$chr_script"  
    echo '#CSUB -n 32' >> "$chr_script"  
    echo '#CSUB -R span[hosts=1]' >> "$chr_script"
    echo '#CSUB -cwd /share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/vcf' >> "$chr_script" 
    echo '' >> "$chr_script"
    echo 'source /share/org/YZWL/yzwl_zhengbin/software/miniconda3/bin/activate' >> "$chr_script"
    echo 'conda activate gwas' >> "$chr_script"

    echo 'gatk --java-options "-Xmx400g" GenotypeGVCFs \' >> "$chr_script"
    echo '-R' ${Ref_genome} '\' >> "$chr_script"
    echo '-V 'chr${chr}.g.vcf.gz' \' >> "$chr_script"
    echo '-O 'chr${chr}.vcf.gz'' >> "$chr_script"
done
