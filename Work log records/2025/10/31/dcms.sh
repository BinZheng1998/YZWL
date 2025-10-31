#!/bin/bash

# 循环遍历当前目录下所有 .windowed.weir.fst 文件
for fst in *.windowed.weir.fst; do
    # 检查文件是否存在（避免处理非文件或空匹配）
    if [[ ! -f "$fst" ]]; then
        continue
    fi
    
    # 从文件名提取前缀，去掉 .windowed.weir.fst 后缀
    prefix=$(basename "$fst" .windowed.weir.fst)
    
    # 构建对应的 pi 文件路径
    pi1="${prefix}_high.windowed.pi"
    pi2="${prefix}_low.windowed.pi"
    
    # 检查两个 pi 文件是否存在，如果缺失则跳过
    if [[ -f "$pi1" && -f "$pi2" ]]; then
        echo "Processing: $fst with $pi1 and $pi2"
        
        # 执行 R 脚本
        /usr/bin/Rscript ~/project/01_evolution/script/DCMS_202506.r \
            --input-fst "$fst" \
            --input-pi1 "$pi1" \
            --input-pi2 "$pi2" \
            --species chicken \
            --out-prefix "$prefix" \
            --output-folder ./
    else
        echo "Skipping $fst: missing $pi1 or $pi2"
    fi
done

echo "All processing completed."
