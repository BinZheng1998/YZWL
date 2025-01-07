#!/bin/bash
#CSUB -J GWAS       
#CSUB -q c01             
#CSUB -o %J.out          
#CSUB -e %J.error        
#CSUB -n 16              
#CSUB -cwd /share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/vcf

Ref_genome="/share/org/YZWL/yzwl_zhengbin/project/00-ref/huxu.fa"

for chr in {1..38} Z W MT; do 
    chr_script="chr${chr}.sh"
    echo '#!/bin/bash' >> "$chr_script"
    echo '#CSUB -J 'gwas_chr${chr}'' >> "$chr_script"  
    echo '#CSUB -q c01' >> "$chr_script" 
    echo '#CSUB -o 'chr${chr}.out'' >> "$chr_script"  
    echo '#CSUB -e 'chr${chr}.err'' >> "$chr_script"  
    echo '#CSUB -n 32' >> "$chr_script" 
    echo '#CSUB -R span[hosts=1]' >> "$chr_script"
    echo '#CSUB -cwd /share/org/YZWL/yzwl_zhengbin/project/01-gwas/01-result/vcf' >> "$chr_script" 
    echo '' >> "$chr_script"
    echo 'source /share/org/YZWL/yzwl_zhengbin/software/miniconda3/bin/activate' >> "$chr_script"
    echo 'conda activate gwas' >> "$chr_script"

    echo 'gatk --java-options "-Xmx300g" CombineGVCFs \' >> "$chr_script"
    echo '--tmp-dir ./tmp \' >> "$chr_script" 
    echo '-R' ${Ref_genome} '\' >> "$chr_script"
    for i in ../gvcf/*.g.vcf; do
        if [ -f "$i" ]; then  
            echo '--variant '${i}' \' >> "$chr_script"
        else
            echo "Warning: No GVCF files found for chromosome ${chr}."
        fi
    done

    echo '-L 'chr${chr}' \' >> "$chr_script"
    echo '-O 'chr${chr}.g.vcf.gz'' >> "$chr_script"
done
