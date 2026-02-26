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

#2
pairtools sort CRR1137296.chromap.pairs --nproc 22 --memory 128G --nproc-in 22 --nproc-out 22 -o CRR1137296.chromap_sorted.pairs

#3
pairtools dedup CRR1137296.chromap_sorted.pairs --max-mismatch 1 --method max --nproc-in 22 --nproc-out 22 -o CRR1137296.chromap_sorted_nodups.pairs
