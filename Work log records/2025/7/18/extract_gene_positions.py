import argparse
import os
from collections import defaultdict

def read_gene_names(gene_list_file):
    """
    从TXT文件读取基因名列表
    :param gene_list_file: 包含基因名的TXT文件路径
    :return: 基因名集合
    """
    gene_names = set()
    try:
        with open(gene_list_file, 'r') as f:
            for line in f:
                gene_name = line.strip()
                if gene_name:
                    gene_names.add(gene_name)
        if not gene_names:
            print("警告: 基因名文件为空")
        return gene_names
    except Exception as e:
        print(f"错误: 读取基因名文件 {gene_list_file} 时发生错误: {e}")
        return None

def parse_gff(gff_file, gene_names):
    """
    从GFF/GTF文件中提取指定基因的位置信息
    :param gff_file: GFF/GTF文件路径
    :param gene_names: 要查询的基因名集合（空集合表示提取所有基因）
    :return: 列表，格式为(染色体, 起始位置, 结束位置, 基因名)
    """
    gene_info = []
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                if fields[2] != 'gene':
                    continue  # 只处理gene类型的记录
                chrom = fields[0]
                start = int(fields[3]) - 1  # GFF/GTF是1-based，转换为0-based
                end = int(fields[4]) - 1    # 结束位置（0-based，闭区间）
                attributes = fields[8]
                # 解析属性字段以获取基因名
                attr_dict = dict(attr.split('=') for attr in attributes.split(';') if '=' in attr)
                gene_name = attr_dict.get('Name') or attr_dict.get('gene') or attr_dict.get('gene_id')
                if gene_name and (not gene_names or gene_name in gene_names):
                    gene_info.append((chrom, start, end, gene_name))
        return gene_info
    except Exception as e:
        print(f"错误: 读取GFF文件 {gff_file} 时发生错误: {e}")
        return None

def write_results(gene_info, output_file):
    """
    将基因位置信息写入TXT文件
    :param gene_info: 基因信息列表
    :param output_file: 输出文件路径
    """
    try:
        with open(output_file, 'w') as f:
            f.write("染色体\t基因开始位置\t基因结束位置\t基因名\n")
            for chrom, start, end, gene_name in sorted(gene_info):
                f.write(f"{chrom}\t{start}\t{end}\t{gene_name}\n")
        print(f"结果已保存至: {output_file}")
    except Exception as e:
        print(f"错误: 写入文件 {output_file} 时发生错误: {e}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="从GFF文件中提取指定基因的位置信息")
    parser.add_argument("gff_file", help="输入的GFF/GTF文件路径")
    parser.add_argument("gene_list_file", help="包含基因名的TXT文件路径")
    parser.add_argument("--output", default="gene_positions.txt", help="输出TXT文件路径（默认为gene_positions.txt）")
    
    args = parser.parse_args()
    
    # 读取基因名列表
    gene_names = read_gene_names(args.gene_list_file)
    if gene_names is None:
        print("错误: 无法继续处理，因为基因名文件解析失败")
        return
    
    # 解析GFF文件，提取基因位置信息
    gene_info = parse_gff(args.gff_file, gene_names)
    if gene_info is None:
        print("错误: 无法继续处理，因为GFF文件解析失败")
        return
    
    if not gene_info:
        print("警告: 未找到任何匹配的基因信息")
        return
    
    # 写入结果
    write_results(gene_info, args.output)

if __name__ == "__main__":
    main()
