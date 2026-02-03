nextflow run nf-core/hic 
-profile singularity \
--input ./sample.csv \
--fasta ~/HIC/ref/chicken/chicken.ref.fa \
--digestion 'hindiii' \
--bwt2_index ~/HIC/ref/chicken \
--split_fastq \
--fastq_chunks_size '10000000' \
--bin_size '20000,40000,100000,250000,1000000' \
--tads_caller 'insulation,hicexplorer' \
-c ~/HIC/script/custom.config \
--res_tads '10000,20000,40000' \
-resume

#
 nextflow run nf-core/hic \
   -profile singularity \
   --input ./sample.csv \
   --fasta ~/HIC/ref/chicken/chicken.ref.fa \
   --digestion 'dpnii' \
   --bwt2_index ~/HIC/ref/chicken \
   --split_fastq \
   --fastq_chunks_size '10000000' \
   --bin_size '20000,40000,100000,250000,1000000' \
   --tads_caller 'insulation,hicexplorer' \
   -c ~/HIC/script/custom2.config \
   --res_tads '10000,20000,40000' \
   -resume
