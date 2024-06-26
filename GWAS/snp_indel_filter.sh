#!/bin/bash

input_vcf="$1"

filename=$(basename -- "$input_vcf")
sample_name="${filename%%_*}"  # Extracts 'sample' from sample_6x.vcf.gz

#snp
# Filter SNPs based on missing data, MAC, and quality
vcftools --gzvcf "$input_vcf" --max-missing 0.5 --mac 3 --minQ 20 --recode --recode-INFO-all --stdout | bgzip -c > "${sample_name}_raw.g5mac3.vcf.gz"

# Filter based on minimum read depth
vcftools --gzvcf "${sample_name}_raw.g5mac3.vcf.gz" --minDP 6 --recode --recode-INFO-all --stdout | bgzip -c > "${sample_name}_raw.g5mac3dp6.vcf.gz"

# Remove indels
vcftools --gzvcf "${sample_name}_raw.g5mac3dp6.vcf.gz" --remove-indels --recode --recode-INFO-all --stdout | bgzip -c > "${sample_name}_snps.vcf.gz"

# Filter biallelic SNPs
vcftools --gzvcf "${sample_name}_snps.vcf.gz" --min-alleles 2 --max-alleles 2 --recode --stdout --recode-INFO-all | bgzip -c > "${sample_name}_SNP.2allele.vcf.gz"

# Final cleaning of SNPs
vcftools --gzvcf "${sample_name}_SNP.2allele.vcf.gz" --max-missing 0.50 --maf 0.05 --hwe 0.000001 --recode --stdout | bgzip -c > "${sample_name}_clean_snps.vcf.gz"

#indel
# Filter indels
vcftools --gzvcf "${sample_name}_raw.g5mac3dp6.vcf.gz" --keep-only-indels --recode --recode-INFO-all --stdout | bgzip -c > "${sample_name}_indels_all.vcf.gz"

# Filter biallelic indels
vcftools --gzvcf "${sample_name}_indels_all.vcf.gz" --min-alleles 2 --max-alleles 2 --recode --stdout --recode-INFO-all | bgzip -c > "${sample_name}_INDEL.2allele.vcf.gz"

# Final cleaning of indels
vcftools --gzvcf "${sample_name}_INDEL.2allele.vcf.gz" --max-missing 0.50 --maf 0.05 --hwe 0.000001 --recode --stdout | bgzip -c > "${sample_name}_clean_indels.vcf.gz"

rm "${sample_name}_raw.g5mac3.vcf.gz" \
        "${sample_name}_raw.g5mac3dp6.vcf.gz" \
        "${sample_name}_snps.vcf.gz" \
        "${sample_name}_SNP.2allele.vcf.gz" \
        "${sample_name}_indels_all.vcf.gz" \
        "${sample_name}_INDEL.2allele.vcf.gz"
