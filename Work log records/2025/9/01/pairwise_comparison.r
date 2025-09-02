#!/bin/bash

BWhigh="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250901_breeds_pairwise_comparisons/0.metadata/body_weight_high"
BWlow="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250901_breeds_pairwise_comparisons/0.metadata/body_weight_low"
vcffile="/home/bzheng/Pop/dna/snp/05_clean_vcf_phased/BW_2036samples/merged_beagle.vcf.gz"
res="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250901_breeds_pairwise_comparisons/2.res"

cd "$BWhigh"
txt1=(`find $(pwd) -type f -name "*.txt"`)
cd "$BWlow"
txt2=(`find $(pwd) -type f -name "*.txt"`)

cd "$res"
#pi
for pop1 in "${txt1[@]}"
do
                prefix1=`basename $pop1 | cut -d "." -f 1`
                vcftools --gzvcf $vcffile --keep $pop1 --window-pi 10000 --window-pi-step 5000 --out "$res"/"$prefix1"_high_pi_10kW_5kS_0901
done

for pop2 in "${txt2[@]}"
do
                prefix2=`basename $pop1 | cut -d "." -f 1`
                vcftools --gzvcf $vcffile --keep $pop2 --window-pi 10000 --window-pi-step 5000 --out "$res"/"$prefix2"_low_pi_10kW_5kS_0901
done

#fst
for pop1 in "${txt1[@]}"; do
    for pop2 in "${txt2[@]}"; do
        prefix1=$(basename "$pop1" .txt)
        prefix2=$(basename "$pop2" .txt)

        echo "Calculating FST between $prefix1 (BWhigh) and $prefix2 (BWlow)"
        vcftools --gzvcf "$vcffile" \
                 --weir-fst-pop $pop1 \
                 --weir-fst-pop $pop2 \
                 --fst-window-size 10000 \
                 --fst-window-step 5000 \
                 --out "$res/${prefix1}_vs_${prefix2}_fst_10kW_5kS"
    done
done
