#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.vcf.gz input_sample.txt"
    exit 1
fi

input_vcf="$1"
pop_samples_file="$2"   

java=/home/bzheng/software/jdk-22.0.1/bin/java
beagle_jar="/home/bzheng/software/beagle.27May24.118.jar"

vcf_files=()

while IFS=" " read -r population samples; do
    echo "Processing population: $population with samples: $samples"

    output_vcf="${population}_subset"

    sample_file="./${population}_samples.txt"
    echo "$samples" | tr ',' '\n' > "$sample_file"
    
    vcftools --gzvcf "$input_vcf" --keep "$sample_file" --recode --out "$output_vcf" 
    wait
    
    bgzip "$output_vcf.recode.vcf"          
    wait

    tabix -p vcf "$output_vcf.recode.vcf.gz"
    wait
    
    beagle_output="./${population}_beagle"
    $java -Xmx280g -jar $beagle_jar gt="${output_vcf}.recode.vcf.gz" impute=true nthreads=16 out=$beagle_output
    wait
    
    tabix -p vcf "${beagle_output}.vcf.gz"
    wait
    
    vcf_files+=("${beagle_output}.vcf.gz")

    rm "$sample_file"

    echo "Beagle imputation and indexing completed for $population. Output saved as ${beagle_output}.vcf.gz"

done < "$pop_samples_file"

merged_vcf="./merged_beagle.vcf.gz"
echo "Concatenating VCF files: ${vcf_files[*]}"

bcftools concat -Oz -o "$merged_vcf" "${vcf_files[@]}"

tabix -p vcf "$merged_vcf"

echo "All populations processed, concatenated, and indexed into ${merged_vcf}"
echo "Concatenation and indexing completed."
