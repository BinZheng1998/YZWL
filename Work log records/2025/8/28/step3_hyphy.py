import os
import subprocess
import json
import glob
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq  # Explicit import for Seq

# 解析命令行参数
parser = argparse.ArgumentParser(description="Run CDS to protein alignment, codon alignment, BUSTED, and RELAX analysis")
parser.add_argument("--input-cds", required=True, help="Directory with CDS files (*.fa)")
parser.add_argument("--tree-file", required=True, help="Path to tree file")
parser.add_argument("--batch", type=int, default=10, help="Number of genes to process concurrently")
parser.add_argument("--output", required=True, help="Output directory for all results")
args = parser.parse_args()

# 定义路径
cds_dir = args.input_cds
output_dir = args.output
protein_output_dir = os.path.join(output_dir, "protein_sequences")
mafft_output_dir = os.path.join(output_dir, "mafft_alignments")
codon_output_dir = os.path.join(output_dir, "codon_alignments")
busted_output_dir = os.path.join(output_dir, "busted_results")
relax_output_dir = os.path.join(output_dir, "relax_results")
tree_file = args.tree_file
busted_summary = os.path.join(output_dir, "busted_summary.txt")
relax_summary = os.path.join(output_dir, "relax_summary.txt")
error_log = os.path.join(output_dir, "error_log.txt")

# 确保输出目录存在
for directory in [output_dir, protein_output_dir, mafft_output_dir, codon_output_dir, busted_output_dir, relax_output_dir]:
    os.makedirs(directory, exist_ok=True)

# 初始化 BUSTED 汇总结果文件
with open(busted_summary, "w") as summary_file:
    summary_file.write("Gene\tp-value\tLRT\tTest_Omega1\tTest_Prop1\tTest_Omega2\tTest_Prop2\tTest_Omega3\tTest_Prop3\tLogL\tAIC-c\n")

# 初始化 RELAX 汇总结果文件
with open(relax_summary, "w") as summary_file:
    summary_file.write("Gene\tp-value\tLRT\tK\tLogL_Null\tLogL_Alternative\tAIC-c_Null\tAIC-c_Alternative\n")

# 初始化错误日志文件
with open(error_log, "w") as error_file:
    error_file.write("Gene\tError_Message\n")

# 安全的 JSON 键访问函数
def safe_get(data, *keys, default="NA"):
    for key in keys:
        data = data.get(key, default) if isinstance(data, dict) else default
        if data == default:
            return default
    return data

# 处理单个基因的函数
def process_gene(cds_file):
    gene_name = os.path.basename(cds_file).replace(".fa", "")
    
    # 1. 使用 transeq 转换 CDS 为蛋白质序列
    protein_output = os.path.join(protein_output_dir, f"{gene_name}_protein.fa")
    transeq_cmd = ["transeq", "-sequence", cds_file, "-outseq", protein_output, "-frame", "1", "-table", "1"]
    try:
        subprocess.run(transeq_cmd, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated protein sequence: {protein_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"transeq error: {e.stderr}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 2. 去除蛋白质序列中的 _1 后缀并将末尾 * 替换为 X
    cleaned_protein_output = os.path.join(protein_output_dir, f"{gene_name}_cleaned_protein.fa")
    try:
        records = []
        with open(protein_output, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # 去除序列 ID 中的 _1 后缀
                record.id = record.id.replace("_1", "")
                # 将序列中所有 * 替换为 X
                seq = str(record.seq).replace("*", "X")
                record.seq = Seq(seq)  # 使用 Bio.Seq.Seq
                records.append(record)
        with open(cleaned_protein_output, "w") as handle:
            SeqIO.write(records, handle, "fasta")
        print(f"Cleaned protein sequence: {cleaned_protein_output}")
    except Exception as e:
        error_msg = f"Protein sequence cleaning error: {str(e)}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 3. 使用 MAFFT 进行蛋白质序列比对
    mafft_output = os.path.join(mafft_output_dir, f"{gene_name}_aligned.fa")
    mafft_cmd = ["mafft", "--localpair", "--maxiterate", "1000", cleaned_protein_output]
    try:
        with open(mafft_output, "w") as outfile:
            subprocess.run(mafft_cmd, stdout=outfile, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated MAFFT alignment: {mafft_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"MAFFT error: {e.stderr}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 4. 使用 pal2nal 生成密码子比对
    codon_output = os.path.join(codon_output_dir, f"{gene_name}_codon.fas")
    pal2nal_cmd = ["perl", "/home/bzheng/software/pal2nal.v14/pal2nal.pl", mafft_output, cds_file, "-nogap", "-output", "fasta"]
    try:
        with open(codon_output, "w") as outfile:
            subprocess.run(pal2nal_cmd, stdout=outfile, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated codon alignment: {codon_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"pal2nal error: {e.stderr}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 5. 运行 HyPhy BUSTED 分析
    busted_hyphy_output = os.path.join(busted_output_dir, f"{gene_name}.json")
    busted_hyphy_cmd = ["hyphy", "/home/bzheng/software/hyphy-analyses/BUSTED-PH/BUSTED-PH.bf", "--srv", "Yes","CPU=1", "--alignment", codon_output, "--tree", tree_file, "--output", busted_hyphy_output, "--branches", "FG"]
    try:
        subprocess.run(busted_hyphy_cmd, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated HyPhy BUSTED output: {busted_hyphy_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"HyPhy BUSTED error: {e.stderr}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 6. 提取 HyPhy BUSTED 结果
    try:
        with open(busted_hyphy_output, "r") as json_file:
            data = json.load(json_file)
        p_value = safe_get(data, "test results", "p-value")
        lrt = safe_get(data, "test results", "LRT")
        log_likelihood = safe_get(data, "fits", "Unconstrained model", "logL")
        aic_c = safe_get(data, "fits", "Unconstrained model", "AIC-c")
        test_rates = safe_get(data, "fits", "Unconstrained model", "rate distributions", "FG", default=[])
        test_omegas = ["NA", "NA", "NA"]
        test_proportions = ["NA", "NA", "NA"]
        for i, rate in enumerate(test_rates[:3]):
            test_omegas[i] = rate.get("omega", "NA")
            test_proportions[i] = rate.get("proportion", "NA")
        with open(busted_summary, "a") as summary_file:
            summary_file.write(f"{gene_name}\t{p_value}\t{lrt}\t{test_omegas[0]}\t{test_proportions[0]}\t"
                              f"{test_omegas[1]}\t{test_proportions[1]}\t{test_omegas[2]}\t{test_proportions[2]}\t"
                              f"{log_likelihood}\t{aic_c}\n")
        print(f"Extracted BUSTED results for {gene_name}")
    except Exception as e:
        error_msg = f"BUSTED JSON processing error: {str(e)}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
    
    # 7. 运行 HyPhy RELAX 分析
    relax_hyphy_output = os.path.join(relax_output_dir, f"{gene_name}.json")
    relax_hyphy_cmd = ["hyphy", "relax", "--alignment", codon_output, "--tree", tree_file, "--output", relax_hyphy_output, "--srv","Yes", "CPU=1" ,"--test", "FG", "--reference", "BG"]
    try:
        subprocess.run(relax_hyphy_cmd, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated HyPhy RELAX output: {relax_hyphy_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"HyPhy RELAX error: {e.stderr}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        return
    
    # 8. 提取 HyPhy RELAX 结果
    try:
        with open(relax_hyphy_output, "r") as json_file:
            data = json.load(json_file)
        p_value = safe_get(data, "test results", "p-value")
        lrt = safe_get(data, "test results", "LRT")
        k_value = safe_get(data, "test results", "relaxation or intensification parameter")
        logl_null = safe_get(data, "fits", "Null", "logL")
        logl_alternative = safe_get(data, "fits", "Alternative", "logL")
        aicc_null = safe_get(data, "fits", "Null", "AIC-c")
        aicc_alternative = safe_get(data, "fits", "Alternative", "AIC-c")
        with open(relax_summary, "a") as summary_file:
            summary_file.write(f"{gene_name}\t{p_value}\t{lrt}\t{k_value}\t{logl_null}\t{logl_alternative}\t"
                              f"{aicc_null}\t{aicc_alternative}\n")
        print(f"Extracted RELAX results for {gene_name}")
    except Exception as e:
        error_msg = f"RELAX JSON processing error: {str(e)}"
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")

# 主程序：并发处理基因，限制同时运行的任务数
cds_files = glob.glob(os.path.join(cds_dir, "*.fa"))
max_concurrent = args.batch

print(f"Found {len(cds_files)} CDS files. Processing up to {max_concurrent} genes concurrently.")
with ProcessPoolExecutor(max_workers=max_concurrent) as executor:
    executor.map(process_gene, cds_files)

# 打印错误日志
if os.path.exists(error_log):
    with open(error_log, "r") as error_file:
        errors = error_file.readlines()[1:]  # 跳过表头
        if errors:
            print("\nGenes with errors:")
            for error in errors:
                gene, msg = error.strip().split("\t", 1)
                print(f"  {gene}: {msg}")
        else:
            print("\nNo errors encountered.")
else:
    print("\nNo errors encountered.")

print(f"BUSTED summary results written to {busted_summary}")
print(f"RELAX summary results written to {relax_summary}")
print(f"Error log written to {error_log}")
