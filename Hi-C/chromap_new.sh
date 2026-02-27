# 1. 直接用 chromap 完成比对、排序、去重 (一步到位)
~/software/chromap/chromap \
  --preset hic \
  --remove-duplicates \
  -t 24 \
  -x ../chromap_index/chicken_ref \
  -r ../chromap_index/chicken.ref.fa \
  -1 SRR28164016_1.clean.fastq.gz \
  -2 SRR28164016_2.clean.fastq.gz \
  -o SRR28164016.pairs.gz  # 直接输出 gz

# 2. 建立索引 (pairix 索引是很多下游工具的必需品)
~/software/pairix/bin/pairix SRR28164016.pairs.gz

# 3. 生成高分辨率 cool 文件 (建议 5000，即 5kb)
cooler cload pairix ../chicken.chrom.sizes:5000 \
  SRR28164016.pairs.gz \
  SRR28164016_5k.cool

# 4. 生成多分辨率 mcool 并自动归一化
# 这里如果不指定 --resolutions，它会默认从 5k 生成到几 Mb
cooler zoomify -n 25 --balance SRR28164016_5k.cool -o SRR28164016.mcool
