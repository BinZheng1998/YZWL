#!/usr/bin/env python3
import os
import subprocess
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

def process_sample(bam_path, gtf_file, output_dir):
    """处理单个样本的 StringTie 定量"""
    # 提取样本名（去掉路径和 'aligned.bam' 后缀）
    bam_basename = os.path.basename(bam_path)
    sample_name = bam_basename.replace('_aligned_sorted.bam', '').rstrip('_')
    
    # 创建样本输出目录
    sample_output_dir = os.path.join(output_dir, sample_name)
    os.makedirs(sample_output_dir, exist_ok=True)
    
    # 构建 StringTie 命令
    output_gtf = os.path.join(sample_output_dir, f"{sample_name}_quant.gtf")
    cmd = [
        "/home/bzheng/software/stringtie-2.2.1.Linux_x86_64/stringtie",
        bam_path,
        "-e", "-B",
        "-G", gtf_file,
        "-o", output_gtf,
        "-p", "8"  # 默认使用 8 个线程，可根据需要调整
    ]
    
    # 运行 StringTie
    print(f"Processing sample: {sample_name}")
    try:
        subprocess.run(cmd, check=True)
        print(f"Quantification completed for {sample_name}, output: {output_gtf}")
    except subprocess.CalledProcessError as e:
        print(f"Error running StringTie for {sample_name}: {e}")
    except FileNotFoundError:
        print("Error: StringTie not found. Please ensure it is installed and in your PATH.")

def run_stringtie_quant(bam_list_file, gtf_file, output_dir, max_samples):
    """并行处理多个样本"""
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 检查输入文件是否存在
    if not os.path.isfile(bam_list_file):
        raise FileNotFoundError(f"BAM list file {bam_list_file} does not exist.")
    if not os.path.isfile(gtf_file):
        raise FileNotFoundError(f"GTF file {gtf_file} does not exist.")
    
    # 读取 BAM 文件路径
    with open(bam_list_file, 'r') as f:
        bam_files = [line.strip() for line in f if line.strip()]
    
    # 使用 ProcessPoolExecutor 进行并行处理
    with ProcessPoolExecutor(max_workers=max_samples) as executor:
        futures = [
            executor.submit(process_sample, bam_path, gtf_file, output_dir)
            for bam_path in bam_files
        ]
        # 等待所有任务完成
        for future in futures:
            future.result()

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="Run StringTie for re-quantification of samples in parallel.")
    parser.add_argument("--input-bam", required=True, help="TXT file containing absolute paths to BAM files")
    parser.add_argument("--gtf", required=True, help="Absolute path to merged GTF file (stringtie_merge.gtf)")
    parser.add_argument("--output-dir", required=True, help="Output directory to store sample-specific results")
    parser.add_argument("--samples", type=int, default=8, help="Number of samples to process in parallel (default: 8)")
    
    args = parser.parse_args()
    
    # 运行主逻辑
    run_stringtie_quant(args.input_bam, args.gtf, args.output_dir, args.samples)

if __name__ == "__main__":
    main()
