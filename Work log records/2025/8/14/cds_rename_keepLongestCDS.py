#!/usr/bin/env python
import argparse
from Bio import SeqIO
import re
from collections import defaultdict

parser = argparse.ArgumentParser(description="Rename FASTA sequence headers to protein IDs and keep the longest sequence for duplicate IDs.")
parser.add_argument("--input", required=True, help="Input FASTA file")
parser.add_argument("--output", required=True, help="Output FASTA file")
args = parser.parse_args()

input_file = args.input
output_file = args.output

sequences = defaultdict(list)

for record in SeqIO.parse(input_file, "fasta"):
    # 提取蛋白质ID（匹配 NP_、XP_ 或 YP_ 开头的ID）
    match = re.search(r'([NXY]P_\d+\.\d+)', record.id)
    if match:
        protein_id = match.group(0)
        sequences[protein_id].append((record, len(record.seq)))
    else:
        print(f"警告：无法从 {record.id} 中提取蛋白质ID，跳过该序列")

unique_records = []
for protein_id, record_list in sequences.items():
    if len(record_list) > 1:
        longest_record = max(record_list, key=lambda x: x[1])[0]
        print(f"发现重复ID {protein_id}，保留长度为 {len(longest_record.seq)} 的序列")
    else:
        longest_record = record_list[0][0]
    longest_record.id = protein_id
    longest_record.description = ""
    unique_records.append(longest_record)

with open(output_file, "w") as output_handle:
    SeqIO.write(unique_records, output_handle, "fasta")

print(f"序列已保存到 {output_file}，共保留 {len(unique_records)} 条序列")
