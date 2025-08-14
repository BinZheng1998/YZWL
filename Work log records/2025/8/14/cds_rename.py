#!/usr/bin/env python
import argparse
from Bio import SeqIO
import re

parser = argparse.ArgumentParser(description="Rename FASTA sequence headers to protein IDs.")
parser.add_argument("--input", required=True, help="Input FASTA file")
parser.add_argument("--output", required=True, help="Output FASTA file")
args = parser.parse_args()

input_file = args.input
output_file = args.output

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        # 提取蛋白质ID（匹配 NP_、XP_ 或 YP_ 开头的ID）
        match = re.search(r'([NXY]P_\d+\.\d+)', record.id)
        if match:
            protein_id = match.group(0)
            record.id = protein_id
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")
        else:
            print(f"警告：无法从 {record.id} 中提取蛋白质ID，跳过该序列")

print(f"序列已保存到 {output_file}")
