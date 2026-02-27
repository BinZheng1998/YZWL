# 1. Chromap 一步到位：比对 + 过滤 + 去重 + 排序 + 压缩
# 使用 --remove-pcr-duplicates 替代 pairtools dedup
# 直接输出 .gz 后缀，chromap 会自动进行 bgzip 压缩
~/software/chromap/chromap \
  --preset hic \
  --remove-pcr-duplicates \
  --low-mem \
  -t 24 \
  -x ../chromap_index/chicken_ref \
  -r ../chromap_index/chicken.ref.fa \
  -1 SRR28164016_1.clean.fastq.gz \
  -2 SRR28164016_2.clean.fastq.gz \
  --pairs \
  -o SRR28164016_nodups.pairs

# 2. 排序
pairtools sort SRR28164016_nodups.pairs --nproc 22 --memory 128G --nproc-in 22 --nproc-out 22 -o SRR28164016_sorted_nodups.pairs

# 3. 建立索引 (pairix)
# 虽然 chromap 输出了 pairs，但 cooler cload pairix 模式需要这个索引来提速
bgzip -@ 24 SRR28164016_sorted_nodups.pairs
~/software/pairix/bin/pairix SRR28164016_sorted_nodups.pairs.gz

# 4. 生成高分辨率 Cooler 文件 (以 5kb 为起始)
# 只有起始分辨率够高，后面的 zoomify 才有意义
cooler cload pairix \
  ../chicken.chrom.sizes:5000 \
  SRR28164016_sorted_nodups.pairs.gz \
  SRR28164016_5k.cool

# 5. 生成多分辨率矩阵并自动归一化 (Balance)
# 这会生成一个包含 5k, 10k, 25k, 50k, 100k... 等所有分辨率的 mcool 文件
cooler zoomify \
  -n 24 \
  --balance \
  SRR28164016_5k.cool \
  -o SRR28164016.mcool
