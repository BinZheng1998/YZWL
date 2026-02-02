#!/bin/bash

# 1. 定义样本列表文件和搜索根目录
SAMPLE_LIST="samples_list.txt"
SEARCH_DIR="../1.data/"

# 2. 获取所有唯一的器官/组织名称
TISSUES=$(awk '{print $2}' $SAMPLE_LIST | sort | uniq)

echo "开始按器官分类合并 Peak 文件..."

for TISSUE in $TISSUES; do
    echo "正在处理器官: $TISSUE"
    
    # 获取属于该器官的所有样本 ID
    SAMPLES=$(awk -v t="$TISSUE" '$2==t {print $1}' $SAMPLE_LIST)
    
    # 创建一个临时文件来存放该器官的所有 Peak 内容
    TMP_FILE="${TISSUE}_combined_raw.tmp"
    > "$TMP_FILE"
    
    for ID in $SAMPLES; do
        # 在指定文件夹下搜索该样本的文件
        # -maxdepth 可以根据你的目录深度调整，以提高速度
        FILE_PATH=$(find $SEARCH_DIR -name "${ID}_peaks.narrowPeak" | head -n 1)
        
        if [ -n "$FILE_PATH" ]; then
            echo "  找到样本 $ID: $FILE_PATH"
            cat "$FILE_PATH" >> "$TMP_FILE"
        else
            echo "  警告: 未找到样本 $ID 的文件"
        fi
    done
    
    # 3. 如果临时文件不为空，则进行合并
    if [ -s "$TMP_FILE" ]; then
        cat "$TMP_FILE" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i stdin > "${TISSUE}_peaks.bed"
        echo "  成功生成: ${TISSUE}_peaks.bed"
    fi
    
    # 删除临时文件
    rm "$TMP_FILE"
done

echo "全部任务完成！"
