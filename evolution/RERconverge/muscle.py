#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import subprocess
import os

# 文件路径
fasta_file = "/home/bzheng/project/10_RER/data/Ref_allspecies.fa"  # 你的FASTA文件
txt_file = "/home/bzheng/project/10_RER/result/RBBH_result/rbh_rename.txt"      # 你的TXT文件
output_dir = "/home/bzheng/project/10_RER/result/protein_MSA_result/alignment_muscle"       # 输出目录

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 1. 读取FASTA文件，存储序列为字典
fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# 2. 读取TXT文件，自动检测分隔符
try:
    df = pd.read_csv(txt_file, sep=None, engine='python', header=0)
except Exception as e:
    print(f"Error reading TXT file: {e}")
    exit(1)

# 获取列名（第一行）
col_names = df.columns[1:].tolist()  # 除第一列外的列名

# 3. 处理每一行
for index, row in df.iterrows():
    group_name = row[0]  # 第一列作为文件名
    proteins = row[1:].tolist()  # 除第一列外的蛋白质名称
    # 过滤掉“-”并计数有效蛋白质
    valid_proteins = [(col, prot) for col, prot in zip(col_names, proteins) if prot != "-"]
    if len(valid_proteins) < 10:  # 如果有效蛋白质数量<10，跳过
        print(f"Skipping {group_name}: only {len(valid_proteins)} proteins (less than 10)")
        continue

    # 检查所有蛋白质名称是否在FASTA文件中存在
    missing_proteins = [prot for _, prot in valid_proteins if prot not in fasta_dict]
    if missing_proteins:
        print(f"Warning: {group_name} has missing proteins in FASTA: {missing_proteins}")
        continue

    # 4. 生成子FASTA文件
    sub_fasta_file = os.path.join(output_dir, f"{group_name}.fa")
    with open(sub_fasta_file, "w") as f:
        for col_name, protein in valid_proteins:
            f.write(f">{col_name}\n{fasta_dict[protein]}\n")
    print(f"Created FASTA file: {sub_fasta_file}")

    # 5. 使用MAFFT进行多序列比对
    aligned_file = os.path.join(output_dir, f"{group_name}_aligned.fa")
    try:
        subprocess.run(["/home/bzheng/software/muscle-linux-x86.v5.3", "-threads", "24", "-align", sub_fasta_file, "-output", aligned_file], check=True)
        print(f"Generated alignment: {aligned_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE for {group_name}: {e}")
        continue

    # 6. 使用trimal修剪比对结果
    trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed.fa")
    try:
        subprocess.run(["trimal", "-in", aligned_file ,"-out", trimmed_file, "-automated1"], check=True)
        print(f"Generated trimmed alignment: {trimmed_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running trimal for {group_name}: {e}")
        continue

print("Processing complete!")
