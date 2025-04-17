import sys
from collections import defaultdict

def parse_gff(gff_file):
    """
    解析GFF文件，按染色体和链存储基因信息
    参考GFF格式规范和解析方法
    """
    gene_dict = defaultdict(list)
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # 跳过注释行
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2].lower() != 'gene':
                continue

            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            gene_id = parts[8].split(';')[0].split('=')[1]  # 解析ID属性

            # 存储基因信息：染色体作为主键，同一染色体下按链分组
            gene_dict[chrom].append((start, end, strand, gene_id))
    return gene_dict

def find_chromosome_overlaps(gene_dict):
    """
    检测同一条染色体上的基因重叠关系
    优化算法：按起始位置排序后滑动窗口比较
    """
    results = []
    for chrom, genes in gene_dict.items():
        # 先按起始位置排序提升效率
        sorted_genes = sorted(genes, key=lambda x: x[0])

        # 滑动窗口两两比较
        for i in range(len(sorted_genes)):
            g1_start, g1_end, g1_strand, g1_id = sorted_genes[i]
            g1_length = g1_end - g1_start + 1

            # 提前终止条件优化
            for j in range(i+1, len(sorted_genes)):
                g2_start, g2_end, g2_strand, g2_id = sorted_genes[j]

                # 提前终止：后续基因起始位置超过当前基因结束位置
                if g2_start > g1_end:
                    break

                # 检查是否同链
                if g1_strand != g2_strand:
                    continue

                # 计算重叠区域
                overlap_start = max(g1_start, g2_start)
                overlap_end = min(g1_end, g2_end)
                overlap_length = overlap_end - overlap_start + 1

                if overlap_length <= 0:
                    continue

                # 判断包含关系
                if (g1_start <= g2_start and g1_end >= g2_end) or \
                   (g2_start <= g1_start and g2_end >= g1_end):
                    relation = '包含'
                else:
                    relation = '重叠'

                # 组装结果行
                g2_length = g2_end - g2_start + 1
                results.append([
                    chrom,
                    g1_id, g1_start, g1_end, g1_strand, g1_length,
                    g2_id, g2_start, g2_end, g2_strand, g2_length,
                    overlap_length, relation
                ])
    return results

def main():
    if len(sys.argv) != 2:
        print("Usage: python find_gene_overlaps.py input.gff")
        sys.exit(1)

    gene_dict = parse_gff(sys.argv[1])
    overlaps = find_chromosome_overlaps(gene_dict)

    # 输出表头
    header = [
        "Chromosome",
        "Gene1_ID", "Gene1_Start", "Gene1_End", "Gene1_Strand", "Gene1_Length",
        "Gene2_ID", "Gene2_Start", "Gene2_End", "Gene2_Strand", "Gene2_Length",
        "Overlap_Length", "Relationship"
    ]
    print('\t'.join(header))

    # 输出结果
    for row in overlaps:
        print('\t'.join(map(str, row)))

if __name__ == '__main__':
    main()
