import pandas as pd
import os
import sys
from tqdm import tqdm

# ================= 配置区域 =================
# 1. 包含6000个 featureCounts 结果文件路径的txt
INPUT_FILE_LIST = "/home/bzheng/project/04_rnaseq/chicken/03_featureCounts/sample_link.txt" 

# 2. 所有基因ID的列表文件 (无表头，每行一个ID)
ALL_GENES_FILE = "/home/bzheng/project/04_rnaseq/chicken/01_ref/chicken_gene_list.txt"

# 3. 输出文件名
OUTPUT_COUNTS = "Matrix_Counts.csv"
OUTPUT_TPM = "Matrix_TPM.csv"
# ===========================================

def get_clean_sample_name(bam_path_column):
    """
    featureCounts 的最后一列通常是 bam 文件的绝对路径
    例如: /home/bzheng/.../CRR087663_aligned_sorted.bam
    我们需要提取: CRR087663
    """
    base_name = os.path.basename(bam_path_column)
    # 假设样本名是第一个点或下划线之前的部分，或者去掉后缀
    # 这里处理常见的 .bam, _aligned, _sorted 等后缀
    clean_name = base_name.replace(".bam", "").replace("_aligned", "").replace("_sorted", "")
    return clean_name

def calculate_tpm_for_sample(df):
    """
    在单样本 DataFrame 内部计算 TPM
    输入 df 必须包含 'Count' 和 'Length' 两列
    """
    # 1. 计算 RPK (Reads Per Kilobase)
    # Length 是 bp，需要除以 1000 变 kb
    # 为了避免除以0，Length加一个极小值(虽然featureCounts通常不会有0长)
    rpk = df['Count'] / (df['Length'] / 1000.0)
    
    # 2. 计算缩放因子 (Per Million)
    scaling_factor = rpk.sum() / 1000000.0
    
    # 3. 计算 TPM
    if scaling_factor > 0:
        tpm = rpk / scaling_factor
    else:
        tpm = rpk * 0 # 如果样本没有任何read，TPM全为0
        
    return tpm

def main():
    # --- 1. 读取基因列表 ---
    print(f"正在读取基因列表: {ALL_GENES_FILE} ...")
    if not os.path.exists(ALL_GENES_FILE):
        print(f"错误: 找不到基因列表文件 {ALL_GENES_FILE}")
        sys.exit(1)
        
    with open(ALL_GENES_FILE, 'r') as f:
        master_gene_list = list(dict.fromkeys([line.strip() for line in f if line.strip()]))
    print(f"标准基因列表已加载: {len(master_gene_list)} 个基因")

    # --- 2. 读取文件路径 ---
    if not os.path.exists(INPUT_FILE_LIST):
        print(f"错误: 找不到路径文件 {INPUT_FILE_LIST}")
        sys.exit(1)

    with open(INPUT_FILE_LIST, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
    print(f"检测到 {len(files)} 个样本文件")

    # 存储结果的列表
    all_counts_series = []
    all_tpm_series = []

    # --- 3. 循环处理 ---
    for file_path in tqdm(files, desc="Processing featureCounts", unit="file"):
        if not os.path.exists(file_path):
            continue
        
        try:
            # featureCounts 包含注释行(#)，需要用 comment='#' 跳过
            # sep='\t' 指定制表符
            # header=0 表示第一行有效数据是表头
            df = pd.read_csv(file_path, sep='\t', comment='#')
            
            # 提取需要的列
            # Geneid 是第1列 (index 0)
            # Length 是第6列 (index 5)
            # Count 是第7列 (index 6, 也是最后一列)
            
            # 获取最后一列的列名（这就是含路径的样本名）
            raw_sample_col = df.columns[-1]
            sample_name = get_clean_sample_name(raw_sample_col)
            
            # 简化 DataFrame，只保留 Geneid, Length, Count
            sub_df = df[['Geneid', 'Length', raw_sample_col]].copy()
            sub_df.columns = ['Geneid', 'Length', 'Count']
            
            # 设为索引
            sub_df.set_index('Geneid', inplace=True)
            
            # --- 关键步骤: 对齐到 Master Gene List ---
            # reindex 会自动引入缺失的基因(填NaN)，并丢弃多余的基因
            # 我们先不填0，因为Length填0会导致TPM计算除以0错误
            aligned_df = sub_df.reindex(master_gene_list)
            
            # Count 缺失补 0
            aligned_df['Count'] = aligned_df['Count'].fillna(0)
            
            # Length 缺失问题：
            # 如果 featureCounts 文件里没有 MasterList 里的某个基因，Length 也会是 NaN。
            # 这种情况下 TPM 无法计算。
            # 通常我们假设所有 featureCounts 文件用的同一个 GTF，Length 应该是一样的。
            # 如果是 NaN，说明该基因在这个样本的 featureCounts 没出现，Length 我们暂且设为 1 (避免报错)，反正Count是0，TPM也是0。
            aligned_df['Length'] = aligned_df['Length'].fillna(1)

            # --- 收集 Counts ---
            count_series = aligned_df['Count'].astype(int)
            count_series.name = sample_name
            all_counts_series.append(count_series)
            
            # --- 计算并收集 TPM ---
            tpm_series = calculate_tpm_for_sample(aligned_df)
            tpm_series.name = sample_name
            # 转为 float32 节省内存
            tpm_series = tpm_series.astype('float32')
            all_tpm_series.append(tpm_series)

        except Exception as e:
            print(f"\n[Error] 处理文件失败 {file_path}: {e}")

    # --- 4. 合并并输出 ---
    
    # 输出 Counts 矩阵
    if all_counts_series:
        print("\n正在合并 Counts 矩阵...")
        counts_matrix = pd.concat(all_counts_series, axis=1)
        print(f"写入: {OUTPUT_COUNTS}")
        counts_matrix.to_csv(OUTPUT_COUNTS)
    
    # 输出 TPM 矩阵
    if all_tpm_series:
        print("正在合并 TPM 矩阵...")
        tpm_matrix = pd.concat(all_tpm_series, axis=1)
        print(f"写入: {OUTPUT_TPM}")
        tpm_matrix.to_csv(OUTPUT_TPM)

    print("完成！")

if __name__ == "__main__":
    main()
