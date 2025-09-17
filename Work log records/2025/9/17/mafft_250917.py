#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import subprocess
import os
from multiprocessing import Pool

# 文件路径
fasta_file = "/home/bzheng/project/11_multic_species_homologs/protein_index_202509/57species_ref.fa"  # FASTA 文件
txt_file = "/home/bzheng/project/11_multic_species_homologs/RBBH_res_202509/31species_rbh_combined.txt"  # TXT 文件
#txt_file = "/home/bzheng/project/10_RER/result/RBBH_result/rbh_rename_20250813.txt"  # TXT 文件
output_dir = "/home/bzheng/project/11_multic_species_homologs/alignment_geneTree_2509"  # 输出目录

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 读取 FASTA 文件，存储序列为字典
fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# 读取 TXT 文件，自动检测分隔符
try:
    df = pd.read_csv(txt_file, sep=None, engine='python', header=0)
except Exception as e:
    print(f"读取 TXT 文件时出错：{e}")
    exit(1)

# 获取列名（除第一列外的列名）
col_names = df.columns[1:].tolist()

# 定义处理单个蛋白质组的函数
def process_group(row_data):
    # 使用 iloc 进行位置索引，避免 FutureWarning
    group_name, proteins = row_data.iloc[0], row_data.iloc[1:].tolist()
    print(f"处理组: {group_name}, 蛋白质: {proteins}")
    valid_proteins = [(col, prot) for col, prot in zip(col_names, proteins) if prot != "-"]
    print(f"有效蛋白质: {valid_proteins}")
    # 检查有效蛋白质数量
    if len(valid_proteins) < 20:
        print(f"跳过 {group_name}：仅有 {len(valid_proteins)} 个蛋白质（少于 10 个）")
        return

    # 检查蛋白质是否在 FASTA 文件中
    missing_proteins = [prot for _, prot in valid_proteins if prot not in fasta_dict]
    if missing_proteins:
        print(f"警告：{group_name} 在 FASTA 文件中缺少蛋白质：{missing_proteins}")
        return

    # 生成子 FASTA 文件
    sub_fasta_file = os.path.join(output_dir, f"{group_name}.fa")
    with open(sub_fasta_file, "w") as f:
        for col_name, protein in valid_proteins:
            sequence = fasta_dict[protein].replace("\n", "")  # 确保序列为单行
            f.write(f">{col_name}\n{sequence}\n")
    print(f"已创建 FASTA 文件：{sub_fasta_file}")

    # 使用 MAFFT 进行多序列比对
    aligned_file = os.path.join(output_dir, f"{group_name}_aligned.fa")
    temp_aligned_file = os.path.join(output_dir, f"{group_name}_aligned_temp.fa")
    try:
        with open(temp_aligned_file, "w") as out:
            subprocess.run(["mafft", "--auto", "--thread", "4", sub_fasta_file], stdout=out, check=True)
        with open(aligned_file, "w") as out_handle:
            for record in SeqIO.parse(temp_aligned_file, "fasta"):
                sequence = str(record.seq).replace("\n", "")  # 确保序列为单行
                out_handle.write(f">{record.id}\n{sequence}\n")
        os.remove(temp_aligned_file)  # 删除临时文件
        print(f"已生成比对文件：{aligned_file}")
    except subprocess.CalledProcessError as e:
        print(f"运行 MAFFT 时出错（{group_name}）：{e}")
        return

    # 使用 trimal 修剪比对结果
    trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed.fa")
    temp_trimmed_file = os.path.join(output_dir, f"{group_name}_trimmed_temp.fa")
    try:
        subprocess.run(["trimal", "-in", aligned_file, "-out", temp_trimmed_file, "-automated1"], check=True)
        with open(trimmed_file, "w") as out_handle:
            for record in SeqIO.parse(temp_trimmed_file, "fasta"):
                sequence = str(record.seq).replace("\n", "")  # 确保序列为单行
                out_handle.write(f">{record.id}\n{sequence}\n")
        os.remove(temp_trimmed_file)  # 删除临时文件
        print(f"已生成修剪后的比对文件：{trimmed_file}")
    except subprocess.CalledProcessError as e:
        print(f"运行 trimal 时出错（{group_name}）：{e}")
        return

    # 使用 IQ-TREE 构建系统发育树
    tree_file = os.path.join(output_dir, f"{group_name}_tree.nwk")
    try:
        subprocess.run([
            "iqtree2",
            "-s", trimmed_file,
            "-m", "MFP",  # 自动选择最佳模型
            "--ufboot", "1000",  # 1000 次超快自展法
            "--alrt", "1000",
            "--abayes",
            "--bnni",
            "-T", "5",  # 每个任务使用 4 个线程
            "-pre", os.path.join(output_dir, f"{group_name}_tree")
        ], check=True)
        print(f"已生成系统发育树：{tree_file}")
    except subprocess.CalledProcessError as e:
        print(f"运行 IQ-TREE 时出错（{group_name}）：{e}")
        return

# 使用进程池并行处理（最多 10 个任务）
if __name__ == "__main__":
    with Pool(processes=30) as pool:
        pool.map(process_group, [row for _, row in df.iterrows()])
    print("处理完成！")
