import pysam
import argparse
from collections import defaultdict

def parse_vcf(vcf_file, window_size):
    """
    统计vcf.gz文件中每window_size长度内的SNP数量
    :param vcf_file: 输入的vcf.gz文件路径
    :param window_size: 窗口大小（bp，例如1000表示1kb，5000表示5kb）
    """
    # 初始化字典，用于存储每个窗口的SNP数量
    snp_counts = defaultdict(int)
    
    try:
        # 打开vcf.gz文件
        with pysam.VariantFile(vcf_file) as vcf:
            for record in vcf:
                # 只处理SNP（REF和ALT长度相等）
                if len(record.ref) == len(record.alts[0]):
                    chrom = record.chrom
                    pos = record.pos
                    # 计算SNP所属的窗口
                    window = pos // window_size
                    # 记录SNP计数
                    snp_counts[(chrom, window)] += 1
    
    except Exception as e:
        print(f"错误: 处理VCF文件时发生错误: {e}")
        return None
    
    return snp_counts

def print_snp_stats(snp_counts, window_size):
    """
    打印SNP数量统计结果
    """
    print(f"\n窗口大小: {window_size/1000}kb")
    print("染色体\t区间开始\t区间结束\tSNP数量")
    for (chrom, window), count in sorted(snp_counts.items()):
        window_start = window * window_size
        window_end = window_start + window_size - 1
        print(f"{chrom}\t{window_start}\t{window_end}\t{count}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="统计GATK生成的vcf.gz文件中每1kb或5kb窗口内的SNP数量")
    parser.add_argument("vcf_file", help="输入的vcf.gz文件路径")
    parser.add_argument("--window", type=int, choices=[1000, 5000, 10000, 20000, 50000], default=1000,
                        help="窗口大小（1000表示1kb，5000表示5kb，默认为1000）")
    
    args = parser.parse_args()
    
    # 统计SNP
    snp_counts = parse_vcf(args.vcf_file, args.window)
    
    if snp_counts:
        # 打印统计结果
        print_snp_stats(snp_counts, args.window)
    else:
        print("未生成统计结果，可能由于输入文件错误或无SNP记录")

if __name__ == "__main__":
    main()
