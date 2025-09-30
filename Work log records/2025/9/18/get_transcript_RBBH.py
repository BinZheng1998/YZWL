import pandas as pd
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Replace '-' in RBBH matrix with protein IDs from blastp results.")
    parser.add_argument('--blastp-tab', required=True, help='Path to blastp result file (tab1, mt6 format, no header, with species info in column 13)')
    parser.add_argument('--RBBH-tab', required=True, help='Path to RBBH matrix file (tab2)')
    parser.add_argument('--out', required=True, help='Path to output modified RBBH matrix')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads for parallel processing')
    return parser.parse_args()

def read_blastp_tab(blastp_file):
    # 读取tab1，无列名，mt6格式（12列）+ 第13列为species_info
    blastp_df = pd.read_csv(blastp_file, sep='\t', header=None)
    # 指定列名，索引0-11为mt6标准列，索引12为species_info
    blastp_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                         'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'species_info']
    
    # 构建比对字典
    human_to_species = defaultdict(dict)  # human_query -> species_ref -> species_protein
    species_to_human = defaultdict(dict)  # species_query -> human_ref -> human_protein
    
    for _, row in blastp_df.iterrows():
        species_info = row['species_info']
        qseqid, sseqid = row['qseqid'], row['sseqid']
        # 解析species_info，格式如 aardvark_ref_human_query 或 human_ref_ruddy_duck_query
        ref, query = species_info.split('_ref_')[0], species_info.split('_query')[0].split('_ref_')[1]
        
        if ref == 'human':
            human_to_species[qseqid][query] = sseqid
        else:
            species_to_human[query][qseqid] = sseqid
            
    return human_to_species, species_to_human

def process_species(human_id, species, tab2_row, human_to_species, species_to_human, replacements):
    # 检查该物种是否为“-”
    species_idx = tab2_row.index[tab2_row == '-'].tolist()
    if species not in species_idx:
        return None
    
    # 查找最佳相互比对
    species_protein = human_to_species.get(human_id, {}).get(species)
    if species_protein:
        # 验证反向比对
        human_protein = species_to_human.get(species, {}).get(species_protein)
        if human_protein == human_id:
            # 找到最佳相互比对，记录替换
            replacements.append((human_id, species, species_protein))
            return species, species_protein
    return None

def main():
    args = parse_args()
    
    # 读取tab2矩阵（假设tab2有列名）
    tab2 = pd.read_csv(args.RBBH_tab, sep='\t')
    species_names = tab2.columns[1:]  # 除第1列（人类蛋白ID）外的物种名
    
    # 读取tab1 blastp结果
    human_to_species, species_to_human = read_blastp_tab(args.blastp_tab)
    
    # 记录替换信息
    replacements = []
    
    # 使用多线程处理每个物种的“-”替换
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for idx, row in tab2.iterrows():
            human_id = row[0]  # 人类蛋白ID
            # 打印该人类蛋白ID对应的“-”物种
            missing_species = row[1:][row[1:] == '-'].index.tolist()
            if missing_species:
                print(f"Human protein {human_id} has '-' for species: {', '.join(missing_species)}")
            
            # 并行处理每个物种
            for species in missing_species:
                futures.append(executor.submit(process_species, human_id, species, row, 
                                             human_to_species, species_to_human, replacements))
        
        # 等待所有任务完成
        for future in futures:
            result = future.result()
            if result:
                species, protein_id = result
                tab2.loc[tab2[species] == '-', species] = protein_id
    
    # 打印替换信息
    if replacements:
        print("\nReplacements made:")
        for human_id, species, protein_id in replacements:
            print(f"Human protein {human_id}, species {species}: '-' replaced with {protein_id}")
    else:
        print("\nNo replacements were made.")
    
    # 保存修改后的tab2
    tab2.to_csv(args.out, sep='\t', index=False)
    print(f"\nModified RBBH matrix saved to {args.out}")

if __name__ == '__main__':
    main()
