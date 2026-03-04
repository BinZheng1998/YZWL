cat *fastq.gz > data.fastq.gz
seqtk seq -A data.fastq.gz > data.fasta 
~/software/pbmm2 index ~/ref/huxu_v23_ref/chicken.ref.fa chicken.mmi
~/software/pbmm2 align ~/ref/huxu_v23_ref/chicken.ref.fa data.fasta res.bam --preset ISOSEQ --sort -j 24 --log-level INFO 2>pbmm2.log &
