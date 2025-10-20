import pandas as pd
from Bio import SeqIO
import subprocess
import os
from multiprocessing import Pool
import logging
from itertools import islice
import re

# 设置日志
logging.basicConfig(filename="alignment.log", level=logging.INFO, 
                    format="%(asctime)s - %(levelname)s - %(message)s")

# 文件路径
fasta_file = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/1.protein_data/10species.fa"
txt_file = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/2.RBBH/rbh_combined4.txt"
output_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/3.protein_alignment"

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 定义标准氨基酸字符集合
valid_aa = set("ACDEFGHIKLMNPQRSTVWY-X?")

# 函数：检查并替换非标准字符
def replace_nonstandard_chars(sequence, group_name, stage):
    invalid_chars = set(sequence) - valid_aa
    if invalid_chars:
        logging.info(f"Non-standard characters found in {group_name} at {stage}: {invalid_chars}")
        # 替换非标准字符为 X
        for char in invalid_chars:
            sequence = sequence.replace(char, "X")
            logging.info(f"Replaced {char} with X in {group_name} at {stage}")
    return sequence

# 读取FASTA文件
try:
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    logging.info(f"Loaded {len(fasta_dict)} sequences from {fasta_file}")
except Exception as e:
    logging.error(f"Error reading FASTA file: {e}")
    exit(1)

# 读取TXT文件
try:
    df = pd.read_csv(txt_file, sep=None, engine='python', header=0)
    logging.info(f"Loaded {len(df)} gene families from {txt_file}")
except Exception as e:
    logging.error(f"Error reading TXT file: {e}")
    exit(1)

col_names = df.columns[1:].tolist()  # 物种名称

# 处理单行的函数
def process_row(row_data):
    group_name, proteins = row_data[0], row_data[1:].tolist()
    valid_proteins = [(col, prot) for col, prot in zip(col_names, proteins) if prot != "-"]
    
    if len(valid_proteins) < 3:
        logging.warning(f"Skipping {group_name}: only {len(valid_proteins)} proteins")
        return
    
    missing_proteins = [prot for _, prot in valid_proteins if prot not in fasta_dict]
    if missing_proteins:
        logging.warning(f"{group_name} has missing proteins: {missing_proteins}")
        return
    
    # 生成子FASTA文件
    sub_fasta_file = os.path.join(output_dir, f"{group_name}.fa")
    try:
        with open(sub_fasta_file, "w") as f:
            for col_name, protein in valid_proteins:
                sequence = fasta_dict[protein].replace("\n", "")
                # 替换非标准字符
                sequence = replace_nonstandard_chars(sequence, group_name, f"sub_fasta for {col_name}")
                f.write(f">{col_name}\n{sequence}\n")
        logging.info(f"Created FASTA file: {sub_fasta_file}")
    except Exception as e:
        logging.error(f"Error creating FASTA for {group_name}: {e}")
        return
    
    # MAFFT比对
    aligned_file = os.path.join(output_dir, f"{group_name}_aligned.fa")
    temp_aligned_file = os.path.join(output_dir, f"{group_name}_aligned_temp.fa")
    try:
        with open(temp_aligned_file, "w") as out:
            subprocess.run(["mafft", "--auto", "--anysymbol", "--thread", "12", sub_fasta_file], 
                          stdout=out, check=True)
        with open(aligned_file, "w") as out_handle:
            for record in SeqIO.parse(temp_aligned_file, "fasta"):
                sequence = str(record.seq).replace("\n", "")
                # 替换非标准字符
                sequence = replace_nonstandard_chars(sequence, group_name, f"aligned_fasta for {record.id}")
                out_handle.write(f">{record.id}\n{sequence}\n")
        os.remove(temp_aligned_file)
        logging.info(f"Generated alignment: {aligned_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"MAFFT error for {group_name}: {e}")
        return
    
    # trimal修剪
    trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed.fa")
    temp_trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed_temp.fa")
    try:
        subprocess.run(["trimal", "-in", aligned_file, "-out", temp_trimmed_file, 
                       "-automated1"], check=True)
        with open(trimmed_file, "w") as out_handle:
            for record in SeqIO.parse(temp_trimmed_file, "fasta"):
                sequence = str(record.seq).replace("\n", "")
                # 替换非标准字符
                sequence = replace_nonstandard_chars(sequence, group_name, f"trimmed_fasta for {record.id}")
                out_handle.write(f">{record.id}\n{sequence}\n")
        os.remove(temp_trimmed_file)
        logging.info(f"Generated trimmed alignment: {trimmed_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"trimal error for {group_name}: {e}")

# 分批处理函数
def process_batch(batch_rows, batch_num):
    logging.info(f"Starting batch {batch_num} with {len(batch_rows)} rows")
    with Pool(processes=50) as pool:  # 调整为服务器核心数/66，例如256核心用4
        pool.map(process_row, batch_rows)
    logging.info(f"Finished batch {batch_num}")

# 主程序
if __name__ == "__main__":
    batch_size = 100  # 每批100行
    processes = 20  # 并行进程数
    total_rows = len(df)
    logging.info(f"Processing {total_rows} rows in batches of {batch_size}")

    # 分批迭代
    batch_num = 0
    for start in range(0, total_rows, batch_size):
        batch_rows = [row for _, row in islice(df.iterrows(), start, start + batch_size)]
        process_batch(batch_rows, batch_num)
        batch_num += 1
    
    logging.info("All batches processed!")
