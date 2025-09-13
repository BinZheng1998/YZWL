import pandas as pd
from Bio import SeqIO
import argparse

def get_protein(input_feature_txt, species, input_protein_faa, output_fa=None, output_txt=None):
    """
    提取特征表中所有CDS记录，从faa文件中提取对应蛋白序列。
    删除完全相同的转录本（重复的protein_id），保留首次出现的记录。
    同一基因的不同转录本（不同protein_id）继续保留，并以_1、_2等后缀区分。

    参数:
    input_feature_txt (str): 特征表txt文件的路径。
    species (str): 物种名称，用于输出文件名和序列ID。
    input_protein_faa (str): 蛋白质faa文件的路径。
    output_fa (str, optional): 输出fa文件的路径，默认为物种名+'.fa'。
    output_txt (str, optional): 输出txt文件的路径，默认为物种名+'.txt'。
    """
    # 设置默认输出文件名
    if output_fa is None:
        output_fa = f"{species}.fa"
    if output_txt is None:
        output_txt = f"{species}.txt"

    # 步骤1：读取并过滤特征表
    df = pd.read_csv(input_feature_txt, sep='\t', header=None)
    # 过滤第1列为 'CDS' 的行
    df_filtered = df[df[0].isin(['CDS'])].copy()

    # 步骤2：提取蛋白ID、基因ID和蛋白长度
    df_selected = df_filtered[[10, 14, 18]].copy()
    df_selected.columns = ['protein_id', 'gene_id', 'protein_length']
    # 将蛋白长度转换为数值类型，非数值转为NaN
    df_selected['protein_length'] = pd.to_numeric(df_selected['protein_length'], errors='coerce')

    # 步骤3：删除重复的protein_id，保留首次出现
    df_unique = df_selected.drop_duplicates(subset=['protein_id'], keep='first').copy()
    # 检查是否有重复的protein_id
    duplicate_proteins = df_selected[df_selected.duplicated(subset=['protein_id'], keep=False)]['protein_id'].unique()
    if duplicate_proteins.size > 0:
        print(f"警告: 以下蛋白ID（转录本ID）重复，将只保留首次出现的记录: {list(duplicate_proteins)}")

    # 步骤4：为转录本编号
    df_unique['transcript_number'] = df_unique.groupby('gene_id').cumcount() + 1
    # 创建新的序列ID列
    df_unique['new_id'] = df_unique.apply(
        lambda row: f"{row['protein_id']} {species} {row['gene_id']}_{row['transcript_number']}", axis=1
    )

    # 步骤5：保存提取的特征数据到txt文件
    # 保存 protein_id, gene_id, protein_length, transcript_number
    df_unique[['protein_id', 'gene_id', 'protein_length', 'transcript_number']].to_csv(
        output_txt, sep='\t', index=False, header=False
    )

    # 步骤6：读取蛋白质序列文件，构建字典
    protein_dict = {record.id: record.seq for record in SeqIO.parse(input_protein_faa, 'fasta')}

    # 步骤7：生成fa文件
    with open(output_fa, 'w') as out_handle:
        for _, row in df_unique.iterrows():
            protein_id = row['protein_id']
            new_id = row['new_id']
            if protein_id in protein_dict:
                # 写入序列，ID 使用 new_id
                out_handle.write(f">{new_id}\n{protein_dict[protein_id]}\n")
            else:
                print(f"警告: 蛋白ID {protein_id} 在 {input_protein_faa} 中未找到")

    print(f"输出fa文件已写入: {output_fa}")
    print(f"输出txt文件已写入: {output_txt}")

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="从特征表和蛋白序列文件中提取所有CDS数据，删除重复的protein_id")
    parser.add_argument('--input-feature-txt', required=True, help="特征表txt文件路径")
    parser.add_argument('--species', required=True, help="物种名称")
    parser.add_argument('--input-protein-faa', required=True, help="蛋白质faa文件路径")
    args = parser.parse_args()

    # 调用主函数
    get_protein(args.input_feature_txt, args.species, args.input_protein_faa)

if __name__ == "__main__":
    main()
