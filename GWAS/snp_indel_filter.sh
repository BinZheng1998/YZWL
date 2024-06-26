#举个例子
#snp
vcftools --gzvcf chr37.vcf.gz --max-missing 0.5 --mac 3 --minQ 20 --recode --recode-INFO-all --stdout | bgzip -c > 240524_chr37.raw.g5mac3.vcf.gz 
vcftools --gzvcf 240524_chr37.raw.g5mac3.vcf.gz --minDP 6 --recode --recode-INFO-all --stdout | bgzip -c > 240524_chr37.raw.g5mac3dp6.vcf.gz
vcftools --gzvcf 240524_chr37.raw.g5mac3dp6.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | bgzip -c > 240524_chr37.snps.vcf.gz
vcftools --gzvcf 240524_chr37.snps.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout --recode-INFO-all | bgzip -c > 240524_chr37.SNP.2allele.vcf.gz 
vcftools --gzvcf 240524_chr37.SNP.2allele.vcf.gz --max-missing 0.90 --maf 0.05 --hwe 0.000001 --recode --stdout | bgzip -c > 240524_chr37_clean_snps.vcf.gz 
#indel
vcftools --gzvcf 240524_chr37.raw.g5mac3dp6.vcf.gz --keep-only-indels --recode --recode-INFO-all --stdout | bgzip -c > 240524_chr37.indels_all.vcf.gz 
vcftools --gzvcf 240524_chr37.indels_all.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout --recode-INFO-all | bgzip -c > 240524_chr37.INDEL.2allele.vcf.gz 
vcftools --gzvcf 240524_chr37.INDEL.2allele.vcf.gz --max-missing 0.90 --maf 0.05 --hwe 0.000001 --recode --stdout | bgzip -c > 240524_chr37_clean_indels.vcf.gz
