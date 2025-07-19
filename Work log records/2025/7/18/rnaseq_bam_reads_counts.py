import pysam
import argparse
import os
from collections import defaultdict
import math

def read_bam_paths(txt_file):
    """
    从TXT文件读取BAM文件的绝对路径
    :param txt_file: 包含BAM文件路径的TXT文件
    :return: BAM文件路径列表
    """
    bam_files = []
    try:
        with open(txt_file, 'r') as f:
            for line in f:
                path = line.strip()
                if path.endswith('.bam') and os.path.exists(path):
                    bam_files.append(path)
                else:
                    print(f"警告: 无效的BAM文件路径或文件不存在: {path}")
        if not bam_files:
            print("错误: TXT文件中没有有效的BAM文件路径")
        return bam_files
    except Exception as e:
        print(f"错误: 读取TXT文件 {txt_file} 时发生错误: {e}")
        return []

def get_chrom_lengths(bam_file):
    """
    获取BAM文件中所有染色体的长度
    :param bam_file: 输入的BAM文件路径
    :return: 字典，键为染色体名，值为染色体长度
    """
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            chrom_lengths = {name: length for name, length in zip(bam.references, bam.lengths)}
        return chrom_lengths
    except Exception as e:
        print(f"错误: 获取BAM文件 {bam_file} 的染色体长度时发生错误: {e}")
        return {}

def parse_bam(bam_file, window_size=5000):
    """
    统计BAM文件中每5kb区间内的reads数量，并记录所有可能的区间
    :param bam_file: 输入的BAM文件路径
    :param window_size: 窗口大小（默认5000bp，即5kb）
    :return: reads数量字典，所有可能的窗口集合
    """
    # 初始化字典，用于存储每个窗口的reads数量
    reads_counts = defaultdict(int)
    
    try:
        # 打开BAM文件
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for read in bam:
                if read.is_unmapped:
                    continue  # 跳过未比对的reads
                chrom = read.reference_name
                start = read.reference_start  # read的起始位置（0-based）
                end = read.reference_end - 1  # read的结束位置（0-based，reference_end是开区间）
                
                # 计算read覆盖的窗口范围
                start_window = start // window_size
                end_window = end // window_size
                
                # 为每个覆盖的窗口计数
                for window in range(start_window, end_window + 1):
                    reads_counts[(chrom, window)] += 1
    
    except Exception as e:
        print(f"错误: 处理BAM文件 {bam_file} 时发生错误: {e}")
        return None, None
    
    # 获取所有可能的窗口（包括reads数量为0的区间）
    chrom_lengths = get_chrom_lengths(bam_file)
    all_windows = set()
    for chrom, length in chrom_lengths.items():
        num_windows = math.ceil(length / window_size)
        for window in range(num_windows):
            all_windows.add((chrom, window))
    
    return reads_counts, all_windows

def write_stats(reads_counts, all_windows, bam_file, window_size, output_dir):
    """
    将统计结果写入TXT文件，包括reads数量为0的区间
    :param reads_counts: 每个窗口的reads数量
    :param all_windows: 所有可能的窗口集合
    :param bam_file: 输入的BAM文件路径（用于生成输出文件名）
    :param window_size: 窗口大小
    :param output_dir: 输出目录
    """
    # 生成输出文件名（基于BAM文件名）
    bam_basename = os.path.basename(bam_file).replace('.bam', '')
    output_file = os.path.join(output_dir, f"{bam_basename}_reads_counts.txt")
    
    try:
        with open(output_file, 'w') as f:
            f.write("染色体\t区间开始\t区间结束\treads数量\n")
            for chrom, window in sorted(all_windows):
                window_start = window * window_size
                window_end = window_start + window_size - 1
                count = reads_counts[(chrom, window)]
                f.write(f"{chrom}\t{window_start}\t{window_end}\t{count}\n")
        print(f"统计结果已保存至: {output_file}")
    except Exception as e:
        print(f"错误: 写入文件 {output_file} 时发生错误: {e}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="统计BAM文件中每5kb区间内的reads数量，包括reads数量为0的区间")
    parser.add_argument("bam_list_file", help="包含BAM文件绝对路径的TXT文件")
    parser.add_argument("--output-dir", default=".", help="输出TXT文件的目录（默认为当前目录）")
    
    args = parser.parse_args()
    
    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 读取BAM文件路径
    bam_files = read_bam_paths(args.bam_list_file)
    
    if not bam_files:
        print("错误: 无法继续处理，因为没有有效的BAM文件路径")
        return
    
    # 处理每个BAM文件
    for bam_file in bam_files:
        print(f"正在处理: {bam_file}")
        reads_counts, all_windows = parse_bam(bam_file)
        
        if reads_counts is not None and all_windows is not None:
            # 写入统计结果
            write_stats(reads_counts, all_windows, bam_file, 5000, args.output_dir)
        else:
            print(f"未生成统计结果，可能由于BAM文件 {bam_file} 错误或无有效reads")

if __name__ == "__main__":
    main()
