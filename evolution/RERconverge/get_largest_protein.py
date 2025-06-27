import pandas as pd
from Bio import SeqIO
import argparse

def get_protein(input_feature_txt, species, input_protein_faa, output_fa=None, output_txt=None):
    """
    input_feature_txt (str): 特征表txt文件的路径。
    species (str): 物种名称，用于输出文件名和序列ID。
    input_protein_faa (str): 蛋白质faa文件的路径。
    output_fa (str, optional): 输出fa文件的路径，默认为物种名+'.fa'。
    output_txt (str, optional): 输出txt文件的路径，默认为物种名+'.txt'。
    """
    if output_fa is None:
        output_fa = f"{species}.fa"
    if output_txt is None:
        output_txt = f"{species}.txt"

    df = pd.read_csv(input_feature_txt, sep='\t', header=None)
    df_filtered = df[df[0].isin(['CDS'])].copy()

    df_selected = df_filtered[[10, 14, 18]].copy()
    df_selected.columns = ['protein_id', 'gene_id', 'protein_length']

    df_selected['protein_length'] = pd.to_numeric(df_selected['protein_length'], errors='coerce')

    df_max_length = df_selected.sort_values(by=['gene_id', 'protein_length'], ascending=[True, False])
    df_max_length = df_max_length.groupby('gene_id').first().reset_index()

    df_max_length.to_csv(output_txt, sep='\t', index=False, header=False)

    protein_dict = {record.id: record.seq for record in SeqIO.parse(input_protein_faa, 'fasta')}

    with open(output_fa, 'w') as out_handle:
        for _, row in df_max_length.iterrows():
            protein_id = row['protein_id']
            gene_id = row['gene_id']
            if protein_id in protein_dict:
                # 创建新的序列ID，格式为 protein_id-gene_id-species
                new_id = f"{protein_id}-{gene_id}-{species}"
                # 写入序列
                out_handle.write(f">{new_id}\n{protein_dict[protein_id]}\n")
            else:
                print(f"警告: 蛋白ID {protein_id} 在 {input_protein_faa} 中未找到")

    print(f"输出fa文件已写入: {output_fa}")
    print(f"输出txt文件已写入: {output_txt}")

def main():
    parser = argparse.ArgumentParser(description="从特征表和蛋白序列文件中提取数据")
    parser.add_argument('--input-feature-txt', required=True, help="特征表txt文件路径")
    parser.add_argument('--species', required=True, help="物种名称")
    parser.add_argument('--input-protein-faa', required=True, help="蛋白质faa文件路径")
    args = parser.parse_args()

    get_protein(args.input_feature_txt, args.species, args.input_protein_faa)

if __name__ == "__main__":
    main()
