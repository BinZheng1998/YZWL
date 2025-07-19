import pysam
import argparse
from collections import defaultdict

def parse_vcf(vcf_file, window_size):
    """
    统计vcf.gz文件中每window_size长度内indel的长度分布
    :param vcf_file: 输入的vcf.gz文件路径
    :param window_size: 窗口大小（bp，例如1000表示1kb，5000表示5kb）
    """
    # 初始化字典，用于存储每个窗口的indel长度分布和总长度
    indel_lengths = defaultdict(lambda: defaultdict(int))
    total_lengths = defaultdict(lambda: {'INS': 0, 'DEL': 0})
    
    try:
        # 打开vcf.gz文件
        with pysam.VariantFile(vcf_file) as vcf:
            for record in vcf:
                # 只处理indel（通过比较REF和ALT的长度）
                if len(record.ref) != len(record.alts[0]):
                    chrom = record.chrom
                    pos = record.pos
                    # 计算indel所属的窗口
                    window = pos // window_size
                    # 计算indel长度
                    if len(record.ref) > len(record.alts[0]):
                        # 删除
                        length = len(record.ref) - len(record.alts[0])
                        indel_type = 'DEL'
                    else:
                        # 插入
                        length = len(record.alts[0]) - len(record.ref)
                        indel_type = 'INS'
                    # 记录indel长度和计数
                    indel_lengths[(chrom, window)][(indel_type, length)] += 1
                    # 更新总长度
                    total_lengths[(chrom, window)][indel_type] += length
    
    except Exception as e:
        print(f"错误: 处理VCF文件时发生错误: {e}")
        return None, None
    
    return indel_lengths, total_lengths

def print_indel_stats(indel_lengths, total_lengths, window_size):
    """
    打印indel长度分布统计结果
    """
    print(f"\n窗口大小: {window_size/1000}kb")
    print("染色体\t窗口起始\t窗口终止\t插入数量\t删除数量\t插入总长度\t删除总长度\t长度分布")
    for (chrom, window) in sorted(indel_lengths.keys()):
        window_start = window * window_size
        window_end = window_start + window_size - 1
        ins_count = sum(count for (t, l), count in indel_lengths[(chrom, window)].items() if t == 'INS')
        del_count = sum(count for (t, l), count in indel_lengths[(chrom, window)].items() if t == 'DEL')
        ins_total_length = total_lengths[(chrom, window)]['INS']
        del_total_length = total_lengths[(chrom, window)]['DEL']
        length_dist = "; ".join(f"{t}:{l}bp={c}" for (t, l), c in sorted(indel_lengths[(chrom, window)].items()))
        print(f"{chrom}\t{window_start}\t{window_end}\t{ins_count}\t{del_count}\t{ins_total_length}\t{del_total_length}\t{length_dist}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="统计GATK生成的vcf.gz文件中每1kb或5kb窗口内indel长度分布及总长度")
    parser.add_argument("vcf_file", help="输入的vcf.gz文件路径")
    parser.add_argument("--window", type=int, choices=[1000, 5000], default=1000,
                        help="窗口大小（1000表示1kb，5000表示5kb，默认为1000）")
    
    args = parser.parse_args()
    
    # 统计indel
    indel_lengths, total_lengths = parse_vcf(args.vcf_file, args.window)
    
    if indel_lengths and total_lengths:
        # 打印统计结果
        print_indel_stats(indel_lengths, total_lengths, args.window)
    else:
        print("未生成统计结果，可能由于输入文件错误或无indel记录")

if __name__ == "__main__":
    main()
