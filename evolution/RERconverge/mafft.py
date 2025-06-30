#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import subprocess
import os

fasta_file = "/home/bzheng/project/10_RER/data/Ref_allspecies.fa"  
txt_file = "/home/bzheng/project/10_RER/result/RBBH_result/rbh_rename.txt"      
output_dir = "/home/bzheng/project/10_RER/result/protein_MSA_result"       

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

try:
    df = pd.read_csv(txt_file, sep=None, engine='python', header=0)
except Exception as e:
    print(f"Error reading TXT file: {e}")
    exit(1)

col_names = df.columns[1:].tolist() 

for index, row in df.iterrows():
    group_name = row[0] 
    proteins = row[1:].tolist()  
    valid_proteins = [(col, prot) for col, prot in zip(col_names, proteins) if prot != "-"]
    if len(valid_proteins) < 10: 
        print(f"Skipping {group_name}: only {len(valid_proteins)} proteins (less than 10)")
        continue

    missing_proteins = [prot for _, prot in valid_proteins if prot not in fasta_dict]
    if missing_proteins:
        print(f"Warning: {group_name} has missing proteins in FASTA: {missing_proteins}")
        continue

    sub_fasta_file = os.path.join(output_dir, f"{group_name}.fa")
    with open(sub_fasta_file, "w") as f:
        for col_name, protein in valid_proteins:
            f.write(f">{col_name}\n{fasta_dict[protein]}\n")
    print(f"Created FASTA file: {sub_fasta_file}")

    aligned_file = os.path.join(output_dir, f"{group_name}_aligned.fa")
    try:
        subprocess.run(["mafft", "--auto","--thread","4", sub_fasta_file], stdout=open(aligned_file, "w"), check=True)
        print(f"Generated alignment: {aligned_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running MAFFT for {group_name}: {e}")
        continue

    trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed.fa")
    try:
        subprocess.run(["trimal", "-in", aligned_file ,"-out", trimmed_file, "-automated1"], check=True)
        print(f"Generated trimmed alignment: {trimmed_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running trimal for {group_name}: {e}")
        continue

print("Processing complete!")
