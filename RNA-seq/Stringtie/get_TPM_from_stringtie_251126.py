import pandas as pd
import os
import sys
from tqdm import tqdm

# ================= 配置区域 =================
# 1. 包含6000个文件绝对路径的txt文件
INPUT_FILE_LIST = "/home/bzheng/project/04_rnaseq/chicken/03_stringtie_2.2.1/re_quant/sample_link.txt" 

# 2. 所有基因ID的列表文件 (无表头，每行一个ID)
ALL_GENES_FILE = "/home/bzheng/project/04_rnaseq/chicken/01_ref/chicken_gene_list.txt"

# 3. 输出文件名
OUTPUT_MATRIX = "TPM_matrix_v2.csv"

# StringTie -A 结果文件的列位置 (0-based)
# 第1列是 Gene ID (索引0)，第9列是 TPM (索引8)
COL_INDEX_ID = 1
COL_INDEX_TPM = 8
# ===========================================

def get_sample_name(file_path):
    """
    从文件路径提取样本名
    """
    base_name = os.path.basename(file_path)
    # 去除文件后缀 (例如 _gene_abund.tsv)
    if "_gene_abund" in base_name:
         sample_name = base_name.split("_gene_abund")[0]
    else:
        sample_name = os.path.splitext(base_name)[0]
    return sample_name

def main():
    # --- 步骤 1: 读取所有基因列表 (Master List) ---
    print(f"正在读取基因列表: {ALL_GENES_FILE} ...")
    if not os.path.exists(ALL_GENES_FILE):
        print(f"错误: 找不到基因列表文件 {ALL_GENES_FILE}")
        sys.exit(1)
    
    with open(ALL_GENES_FILE, 'r') as f:
        # 读取并去重，确保Master List本身唯一
        master_gene_list = list(dict.fromkeys([line.strip() for line in f if line.strip()]))
    
    print(f"标准基因列表已加载，共包含 {len(master_gene_list)} 个基因。")

    # --- 步骤 2: 读取样本路径 ---
    if not os.path.exists(INPUT_FILE_LIST):
        print(f"错误: 找不到路径文件 {INPUT_FILE_LIST}")
        sys.exit(1)

    with open(INPUT_FILE_LIST, 'r') as f:
        files = [line.strip() for line in f if line.strip()]
    
    print(f"检测到 {len(files)} 个样本文件，准备处理...")

    processed_columns = []
    
    # --- 步骤 3: 循环处理每个文件 ---
    for file_path in tqdm(files, desc="Processing", unit="file"):
        if not os.path.exists(file_path):
            continue

        try:
            s_name = get_sample_name(file_path)
            
            # 读取数据: 只读 ID 和 TPM 列
            df = pd.read_csv(file_path, delim_whitespace=True, comment='#', usecols=[COL_INDEX_ID, COL_INDEX_TPM])
            df.columns = ['Gene_ID', 'TPM']

            # 1. 清洗数据：去除表头行 (如果在中间出现)
            df = df[df['Gene_ID'] != 'Gene ID'] 
            df = df[df['Gene_ID'] != 'Gene'] 

            # 2. 转换数值：无法转换的变NaN
            df['TPM'] = pd.to_numeric(df['TPM'], errors='coerce')
            
            # 3. 去除 NaN (可能是原来的表头行转换失败产生的)
            df = df.dropna(subset=['TPM'])

            # ========================================================
            # 关键修复: 处理重复 Gene ID
            # ========================================================
            # 如果同一个ID出现多次，groupby().sum() 将把它们的TPM加起来
            # 这样保证了索引的唯一性
            df_unique = df.groupby('Gene_ID')['TPM'].sum()
            
            # 4. 对齐到 Master Gene List
            # reindex 现在安全了，因为 df_unique 的索引绝对唯一
            aligned_series = df_unique.reindex(master_gene_list, fill_value=0)
            
            # 5. 格式化
            aligned_series.name = s_name
            aligned_series = aligned_series.astype('float32')

            processed_columns.append(aligned_series)

        except Exception as e:
            # 打印具体哪个文件报错，但不中断整个程序
            print(f"\n[Warning] 跳过文件 {s_name}: {e}")

    # --- 步骤 4: 合并并输出 ---
    print("\n正在合并所有样本...")
    
    if processed_columns:
        final_df = pd.concat(processed_columns, axis=1)
        
        print(f"正在写入 CSV: {OUTPUT_MATRIX} ...")
        final_df.to_csv(OUTPUT_MATRIX)
        print("完成！")
    else:
        print("没有有效数据被合并。")

if __name__ == "__main__":
    main()
