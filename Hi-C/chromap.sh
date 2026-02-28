#!/bin/bash

# 设置路径变量（请根据你的实际路径修改）
CHROMAP=~/software/chromap/chromap
PAIRIX=~/software/pairix/bin/pairix
INDEX=~/HIC/chicken/chromap_index/chicken_ref
FASTA=~/HIC/chicken/chromap_index/chicken.ref.fa
CHROMSIZES=~/HIC/chicken/chicken.chrom.sizes

# 遍历当前文件夹下所有 R1 文件
for r1_file in *_1.clean.fastq.gz; do
    # 提取样本前缀 (例如从 SRR123_1.clean.fastq.gz 提取 SRR123)
    sample_id=${r1_file%_1.clean.fastq.gz}
    r2_file="${sample_id}_2.clean.fastq.gz"
    
    # 检查 R2 文件是否存在
    if [[ ! -f "$r2_file" ]]; then
        echo "Warning: R2 file for $sample_id not found, skipping..."
        continue
    fi

    echo "======= Processing Sample: $sample_id ======="

    # 1. Chromap 比对、去重、输出 pairs
    echo "Step 1: Mapping with chromap..."
    $CHROMAP \
      --preset hic \
      --remove-pcr-duplicates \
      --low-mem \
      -t 24 \
      -x $INDEX \
      -r $FASTA \
      -1 "$r1_file" \
      -2 "$r2_file" \
      --pairs \
      -o "${sample_id}_nodups.pairs"

    # 2. 排序 (pairtools sort)
    echo "Step 2: Sorting pairs..."
    pairtools sort \
      "${sample_id}_nodups.pairs" \
      --nproc 22 \
      --memory 128G \
      -o "${sample_id}_sorted_nodups.pairs"

    # 3. 压缩并建立索引
    echo "Step 3: Compressing and indexing..."
    bgzip -@ 24 "${sample_id}_sorted_nodups.pairs"
    $PAIRIX "${sample_id}_sorted_nodups.pairs.gz"

    # 4. 生成 5kb 分辨率的 Cooler
    echo "Step 4: Creating 5kb cooler..."
    cooler cload pairix \
      "$CHROMSIZES:5000" \
      "${sample_id}_sorted_nodups.pairs.gz" \
      "${sample_id}_5k.cool"

    # 5. 生成多分辨率 mcool 并归一化
    echo "Step 5: Zoomifying and balancing..."
    cooler zoomify \
      -n 24 \
      --balance \
      "${sample_id}_5k.cool" \
      -o "${sample_id}.mcool"

    # 可选：清理中间大的文本文件以节省空间
    # rm "${sample_id}_nodups.pairs" 
    
    echo "======= Finished Sample: $sample_id ======="
done
