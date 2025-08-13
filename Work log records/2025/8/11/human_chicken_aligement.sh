#!/bin/bash
human_dir="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250811_conservative_analysis/data"
output_dir="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250811_conservative_analysis/MSA"

for chr in "$human_dir"/chr*.fa;do
        chr_name=$(basename "$chr".fa)
        echo "processing $chr_name"
        lastz "$chr" /home/bzheng/project/01_evolution/convergent_evolution/07_result/20250811_conservative_analysis/data/chicken.ref.fa \
        --format=maf \
        --output="$output_dir/${chr_name}_chicken.maf"
        echo "finished $chr_name"
done
