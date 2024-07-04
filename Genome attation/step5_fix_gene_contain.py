#!/usr/bin/env python3

import sys

def parse_attributes(attributes):
    """解析属性字符串为字典"""
    return dict(attr.split('=') for attr in attributes.split(';') if attr)

def format_attributes(attributes):
    """将字典格式的属性转换回字符串"""
    return ';'.join(f'{key}={value}' for key, value in attributes.items())

def main():
    input_file = sys.argv[1]
    records = []
    genes = []

    # 读取GFF文件
    with open(input_file, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                records.append(parts)
                if parts[2] == 'gene':
                    genes.append(parts)

    # 按染色体和起始位置对基因进行排序
    genes.sort(key=lambda x: (x[0], int(x[3])))

    # 查找包含关系并记录
    covered_genes = {}
    gene_parent_map = {}
    for i in range(len(genes)):
        chrom1, start1, end1, strand1, attrs1 = genes[i][0], int(genes[i][3]), int(genes[i][4]), genes[i][6], genes[i][8]
        gene_id1 = parse_attributes(attrs1)['ID']
        
        for j in range(i + 1, len(genes)):
            chrom2, start2, end2, strand2, attrs2 = genes[j][0], int(genes[j][3]), int(genes[j][4]), genes[j][6], genes[j][8]
            gene_id2 = parse_attributes(attrs2)['ID']
            
            if chrom1 != chrom2 or strand1 != strand2:
                continue
            
            if start1 <= start2 and end1 >= end2:
                # Gene1包含Gene2
                covered_genes[gene_id2] = gene_id1
            elif start2 <= start1 and end2 >= end1:
                # Gene2包含Gene1
                covered_genes[gene_id1] = gene_id2

    # 更新Parent映射表
    for rec in records:
        if rec[2] in ['mRNA', 'lnc_RNA', 'tRNA', 'snoRNA', 'ncRNA', 'rRNA', 'snRNA', 'SRP_RNA']:
            gene_id = parse_attributes(rec[8]).get('Parent')
            if gene_id in covered_genes:
                # 更新被覆盖基因的RNA行的Parent
                original_gene_id = covered_genes[gene_id]
                # 查找最长基因的Parent ID
                while original_gene_id in covered_genes:
                    original_gene_id = covered_genes[original_gene_id]
                gene_parent_map[gene_id] = original_gene_id

    # 输出处理后的结果
    for rec in records:
        if rec[2] == 'gene' and parse_attributes(rec[8])['ID'] in covered_genes:
            continue  # 跳过被包含的基因行
        if rec[2] in ['mRNA', 'lnc_RNA', 'tRNA', 'snoRNA', 'ncRNA', 'rRNA', 'snRNA', 'SRP_RNA']:
            attr = parse_attributes(rec[8])
            parent_id = attr.get('Parent')
            if parent_id in gene_parent_map:
                attr['Parent'] = gene_parent_map[parent_id]
                rec[8] = format_attributes(attr)
        print("\t".join(rec))

if __name__ == "__main__":
    main()
