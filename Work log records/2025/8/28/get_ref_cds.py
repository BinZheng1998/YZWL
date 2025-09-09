#!/usr/bin/env python
import argparse
from Bio import SeqIO
import re
import os

# 设置命令行参数解析
parser = argparse.ArgumentParser(description="Rename and combine FASTA sequence headers to protein IDs.")
parser.add_argument("--input", required=True, help="Input text file containing absolute paths to FASTA files")
parser.add_argument("--output", required=True, help="Output combined FASTA file")
args = parser.parse_args()

# 获取输入和输出文件路径
input_txt = args.input
output_file = args.output

# 确保输出文件的目录存在
output_dir = os.path.dirname(output_file)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

# 读取文本文件中的FASTA文件路径
fasta_files = []
with open(input_txt, "r") as txt_handle:
    for line in txt_handle:
        fasta_path = line.strip()
        if fasta_path and os.path.exists(fasta_path):
            fasta_files.append(fasta_path)
        else:
            print(f"警告：文件 {fasta_path} 不存在，跳过")

# 打开输出文件
with open(output_file, "w") as output_handle:
    # 遍历每个FASTA文件
    for fasta_file in fasta_files:
        try:
            # 读取FASTA文件
            for record in SeqIO.parse(fasta_file, "fasta"):
                # 提取蛋白质ID（匹配 NP_、XP_ 或 YP_ 开头的ID）
                match = re.search(r'([NXY]P_\d+\.\d+)', record.id)
                if match:
                    protein_id = match.group(0)
                    # 更新序列的ID
                    record.id = protein_id
                    record.description = ""  # 清空描述部分，只保留ID
                    # 写入新的FASTA文件
                    SeqIO.write(record, output_handle, "fasta")
                else:
                    print(f"警告：无法从 {fasta_file} 的 {record.id} 中提取蛋白质ID，跳过该序列")
        except Exception as e:
            print(f"错误：处理文件 {fasta_file} 时发生异常 - {str(e)}")

print(f"所有序列已合并并保存到 {output_file}")
