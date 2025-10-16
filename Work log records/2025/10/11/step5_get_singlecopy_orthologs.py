#!/usr/bin/env python3
import pandas as pd
import os
import shutil
import argparse

def copy_gene_files(input_txt, folder_a, folder_b):
    # 确保输出文件夹b存在
    if not os.path.exists(folder_b):
        os.makedirs(folder_b)

    # 读取txt文件，并过滤掉任何一列包含'-'的行
    try:
        df = pd.read_csv(input_txt, sep='\t')  # 读取整个txt文件
        # 过滤掉任何一列包含'-'的行
        mask = df.apply(lambda row: row.str.contains('-', na=False).any(), axis=1)
        df = df[~mask]  # 保留不含'-'的行
        gene_list = df.iloc[:, 0].tolist()  # 获取第一列（基因名）
    except Exception as e:
        print(f"Error reading {input_txt}: {e}")
        return

    # 遍历基因列表，查找并拷贝文件
    for gene in gene_list:
        source_file = os.path.join(folder_a, f"{gene}_trimmed.fa")
        destination_file = os.path.join(folder_b, f"{gene}_trimmed.fa")
        
        # 检查源文件是否存在
        if os.path.exists(source_file):
            shutil.copy(source_file, destination_file)
            print(f"Copied: {gene}_trimmed.fa")
        else:
            print(f"File not found: {gene}_trimmed.fa")

    print("Copy operation completed!")

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description="Copy gene files from folder_a to folder_b based on gene names in the first column of a txt file, excluding rows with '-' in any column.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input txt file containing gene names in the first column")
    parser.add_argument('-a', '--folder_a', required=True, help="Path to the source folder containing .trimmed.fa files")
    parser.add_argument('-b', '--folder_b', required=True, help="Path to the destination folder")

    args = parser.parse_args()

    # 调用主功能
    copy_gene_files(args.input, args.folder_a, args.folder_b)

if __name__ == "__main__":
    main()
