#1
~/software/chromap/chromap \
  --preset hic \
  --low-mem \
  -t 24 \
  -x ../chromap_index/chicken_ref \
  -r ../chromap_index/chicken.ref.fa \
  -1 SRR28164016_1.clean.fastq.gz \
  -2 SRR28164016_2.clean.fastq.gz \
  -o SRR28164016.chromap.pairs
