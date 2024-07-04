#!/usr/bin/env python3

import sys

def parse_gff(gff_file):
    """
    解析GFF文件，返回基因的起始位置、终止位置、染色体号和ID的列表
    """
    genes = []
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # 忽略注释行
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'gene':
                    start = int(parts[3])
                    end = int(parts[4])
                    chromosome = parts[0]
                    gene_id = parts[8].split(';')[0].split('=')[1]  # 解析第9列中的ID
                    genes.append((chromosome, start, end, gene_id))
    return genes

def find_overlapping_genes(genes):
    """
    查找重叠的基因并整理关系
    """
    overlapping_genes = []
    genes.sort(key=lambda x: (x[0], x[1]))  # 按照染色体号和基因的起始位置进行排序
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            gene1 = genes[i]
            gene2 = genes[j]
            if gene2[0] != gene1[0]:
                break  # 已经到达下一个染色体，退出内层循环
            if gene2[1] >= gene1[1] and gene2[2] <= gene1[2]:
                overlap_length = gene2[2] - gene2[1]
                overlapping_genes.append((gene1, gene2, "包含", overlap_length))
            elif gene1[1] >= gene2[1] and gene1[2] <= gene2[2]:
                overlap_length = gene1[2] - gene1[1]
                overlapping_genes.append((gene1, gene2, "被包含", overlap_length))
            elif gene2[1] < gene1[2] and gene2[2] > gene1[1]:
                overlap_length = min(gene1[2], gene2[2]) - max(gene1[1], gene2[1])
                overlapping_genes.append((gene1, gene2, "重叠", overlap_length))
    return overlapping_genes

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python findoverlap.py input.gff")
        sys.exit(1)

    gff_file = sys.argv[1]
    
    # 解析GFF文件并获取基因位置、染色体号和ID信息
    genes = parse_gff(gff_file)
    
    # 查找重叠的基因并整理关系
    overlapping_genes = find_overlapping_genes(genes)
    
    # 输出重叠基因的信息到标准输出
    if overlapping_genes:
        for gene_pair in overlapping_genes:
            gene1 = gene_pair[0]
            gene2 = gene_pair[1]
            relationship = gene_pair[2]
            overlap_length = gene_pair[3]
            print(f"基因1：染色体号：{gene1[0]}, 位置：{gene1[1]}-{gene1[2]}, ID：{gene1[3]}")
            print(f"基因2：染色体号：{gene2[0]}, 位置：{gene2[1]}-{gene2[2]}, ID：{gene2[3]}")
            print(f"关系：{relationship}, 重叠长度：{overlap_length} bp")
            print()
    else:
        print("没有找到重叠的基因。")
