import argparse
import os
import subprocess
from multiprocessing import Pool
import traceback

def run_pipeline_in_env(env_name, pipeline_script, input_files, sample_name):
    """

    :param env_name: Conda 环境名称（例如 'rnaseq' 或 'fusion'）
    :param pipeline_script: pipeline 脚本的路径（例如 'rnaseq_pipeline.sh'）
    :param input_files: 输入文件列表（单端或双端）
    :param sample_name: 样本名
    """
    try:
        cmd = f"conda run -n {env_name} {pipeline_script} {' '.join(input_files)}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"成功在 {env_name} 环境中执行 {pipeline_script} 处理 {input_files}")
    except subprocess.CalledProcessError as e:
        error_message = f"执行 {pipeline_script} 处理 {input_files} 时出错: {e}\n{traceback.format_exc()}"
        print(error_message)
        error_file = f"{sample_name}_error.txt"
        with open(error_file, 'w') as f:
            f.write(error_message)
        print(f"错误信息已记录到 {error_file}")

def extract_sample_name(input_files):
    """
    从输入文件路径中提取样本名。
    """
    basename = os.path.basename(input_files[0])
    sample_name = os.path.splitext(basename)[0].replace('_1', '').replace('_2', '')
    return sample_name

def process_sample(sample, rnaseq, fusion):
    """

    :param sample: 样本文件路径列表（单端或双端）
    :param rnaseq: 是否执行 RNA-seq 分析 ('yes' 或 'no')
    :param fusion: 是否执行 Fusion Gene 分析 ('yes' 或 'no')
    """
    sample_name = extract_sample_name(sample)
    if rnaseq == 'yes':
        run_pipeline_in_env('rnaseq', '/home/bzheng/project/06_rnaseq/02_script/rnaseq_pipeline.sh', sample, sample_name)
    if fusion == 'yes':
        run_pipeline_in_env('fusion', '/home/bzheng/project/06_rnaseq/02_script/fusion_gene_pipeline.sh', sample, sample_name)

def main():
    parser = argparse.ArgumentParser(description="自动化执行 RNA-seq 和 Fusion Gene 分析流程")
    parser.add_argument('--input-fastq-file', required=True, help="输入文件路径（如 sample.txt）")
    parser.add_argument('--sample-numbers', type=int, default=1, help="并行处理的样本数量，默认为1")
    parser.add_argument('--rnaseq', default='yes', choices=['yes', 'no'], help="是否执行 RNA-seq 分析，默认为yes")
    parser.add_argument('--fusion', default='yes', choices=['yes', 'no'], help="是否执行 Fusion Gene 分析，默认为yes")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_fastq_file):
        print(f"错误：输入文件 {args.input_fastq_file} 不存在")
        return
    
    samples = []
    with open(args.input_fastq_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 1:
                # 单端测序
                samples.append([parts[0]])
            elif len(parts) == 2:
                # 双端测序
                samples.append(parts)
            else:
                print(f"警告：输入文件中的行格式错误: {line.strip()}")

    tasks = [(sample, args.rnaseq, args.fusion) for sample in samples]
    
    with Pool(processes=args.sample_numbers) as pool:
        pool.starmap(process_sample, tasks)

if __name__ == "__main__":
    main()
