#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import os
import sys

# 文件路径
fasta_file = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250828_hyphy/1.cds_data/57species_cds_ref.fa"  # 你的FASTA文件
txt_file = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250828_hyphy/0.script/rbh_31species_noMouse.txt"  # 你的TXT文件
output_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250828_hyphy/2.cds_31species"  # 输出目录
min_proteins = 30  # 最小蛋白质数量阈值

# 检查文件是否存在
if not os.path.exists(fasta_file):
    print(f"Error: FASTA file {fasta_file} does not exist")
    sys.exit(1)
if not os.path.exists(txt_file):
    print(f"Error: TXT file {txt_file} does not exist")
    sys.exit(1)

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

# 1. 读取FASTA文件，存储序列为字典
try:
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
except Exception as e:
    print(f"Error parsing FASTA file: {e}")
    sys.exit(1)

# 2. 读取TXT文件，明确指定分隔符（如制表符）
try:
    df = pd.read_csv(txt_file, sep='\t', header=0)
except Exception as e:
    print(f"Error reading TXT file: {e}")
    sys.exit(1)

# 获取列名（第一列为group_name，其余为物种名）
col_names = df.columns[1:].tolist()

# 3. 记录未提取的蛋白质
missing_proteins_log = []

# 4. 处理每一行
for index, row in df.iterrows():
    group_name = row[0]  # 第一列作为文件名
    proteins = row[1:].tolist()  # 除第一列外的蛋白质名称
    # 过滤掉“-”并计数有效蛋白质
    valid_proteins = [(col, prot) for col, prot in zip(col_names, proteins) if prot != "-"]
    
    # 检查有效蛋白质数量
    if len(valid_proteins) < min_proteins:
        print(f"Skipping {group_name}: only {len(valid_proteins)} proteins (less than {min_proteins})")
        continue

    # 检查蛋白质名称是否在FASTA文件中存在
    missing_proteins = [prot for _, prot in valid_proteins if prot not in fasta_dict]
    if missing_proteins:
        missing_proteins_log.append((group_name, missing_proteins))
        print(f"Warning: {group_name} has missing proteins in FASTA: {missing_proteins}")
        continue

    # 5. 生成子FASTA文件
    sub_fasta_file = os.path.join(output_dir, f"{group_name}.fa")
    if os.path.exists(sub_fasta_file):
        print(f"Warning: {sub_fasta_file} already exists, overwriting")
    
    try:
        with open(sub_fasta_file, "w") as f:
            for col_name, protein in valid_proteins:
                f.write(f">{col_name}\n{fasta_dict[protein]}\n")
        print(f"Created FASTA file: {sub_fasta_file}")
    except Exception as e:
        print(f"Error writing FASTA file {sub_fasta_file}: {e}")
        continue

# 6. 输出未提取的蛋白质
if missing_proteins_log:
    print("\n=== Missing Proteins Summary ===")
    for group_name, proteins in missing_proteins_log:
        print(f"Group {group_name}: {', '.join(proteins)}")
    # 保存未提取蛋白质到文件
    missing_log_file = os.path.join(output_dir, "missing_proteins.txt")
    with open(missing_log_file, "w") as f:
        for group_name, proteins in missing_proteins_log:
            f.write(f"Group {group_name}: {', '.join(proteins)}\n")
    print(f"Missing proteins saved to: {missing_log_file}")
else:
    print("\nNo missing proteins found.")

print("Processing complete.")
