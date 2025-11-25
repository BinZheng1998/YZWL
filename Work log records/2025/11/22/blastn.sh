#!/bin/bash

# ================= 配置区域 =================
# 1. 数据目录 (包含 .fa 文件 和 已经建好的 blast 数据库)
DATA_DIR="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251120_selection_region_convergent/1.blast_db"

# 2. 结果输出目录
OUT_DIR="/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251120_selection_region_convergent/2.blast_result"
mkdir -p "$OUT_DIR"

# 3. 物种/文件名列表
# 请确保 DATA_DIR 下存在对应的 .fa 文件和数据库索引
species=(
    "sheep_BW_pos_top5"
    "pig_BW_pos_top5"
    "cattle_BW_pos_top5"
    "chicken_BW_pos_top5"
    "dog_BW_pos_top5"
    # 在此添加更多...
)

# 4. BLAST 参数
BLAST_CMD="/home/bzheng/software/ncbi-blast-2.15.0+/bin/blastn"
THREADS=64
EVALUE="1e-4"
OUTFMT="6"
MAX_TARGET_SEQS=1

# ================= 全对全相互比对 =================
echo ">>> 开始全对全比对 (跳过建库步骤)..."

# 外层循环：Query (查询序列)
for query_sp in "${species[@]}"; do
    # 内层循环：Reference/DB (参考数据库)
    for db_sp in "${species[@]}"; do
        
        # 1. 跳过自己比自己
        if [ "$query_sp" == "$db_sp" ]; then
            continue
        fi

        # 2. 定义输出文件: DB_ref_Query_query.txt
        output_file="${OUT_DIR}/${db_sp}_ref_${query_sp}_query.txt"

        # 3. 断点续跑检查
        if [ -f "$output_file" ]; then
            echo "  [已存在] 跳过: ${db_sp} (Ref) vs ${query_sp} (Query)"
            continue
        fi

        echo "  [运行中] 参考(DB)=${db_sp} vs 查询(Query)=${query_sp} ..."

        # 4. 运行 BLASTN
        # 注意：这里默认使用 -db 模式，需要 DATA_DIR 下有索引文件
        # 如果你想纯 FASTA 对比 (不依赖索引)，请看脚本下方的说明
        $BLAST_CMD \
            -query "${DATA_DIR}/${query_sp}.fa" \
            -db "${DATA_DIR}/${db_sp}" \
            -out "$output_file" \
            -evalue $EVALUE \
            -outfmt $OUTFMT \
            -num_threads $THREADS \
            -max_target_seqs $MAX_TARGET_SEQS

        if [ $? -ne 0 ]; then
            echo "  [错误] BLAST 运行失败: Ref=${db_sp}, Query=${query_sp}"
        fi

    done
done

echo ">>> 所有任务完成！"
