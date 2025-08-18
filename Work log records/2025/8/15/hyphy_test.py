import os
import subprocess
import json
import glob
from pathlib import Path

# 定义路径
protein_msa_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250816_RER/mafft/alignment"
cds_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250814_hyphy/3.cds_12species"
codon_output_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250814_hyphy/4.codon_alignment"
hyphy_output_dir = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250814_hyphy/5.busted_res"
tree_file = "/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250814_hyphy/1.tree/12species_tree.nwk"
pal2nal_script = "/home/bzheng/software/pal2nal.v14/pal2nal.pl"
output_summary = "busted_summary.txt"
error_log = "error_log.txt"

# 确保输出目录存在
os.makedirs(codon_output_dir, exist_ok=True)
os.makedirs(hyphy_output_dir, exist_ok=True)

# 初始化汇总结果文件
with open(output_summary, "w") as summary_file:
    summary_file.write("Gene\tp-value\tLRT\tTest_Omega1\tTest_Prop1\tTest_Omega2\tTest_Prop2\tTest_Omega3\tTest_Prop3\tLogL\tAIC-c\n")

# 初始化错误日志文件
with open(error_log, "w") as error_file:
    error_file.write("Gene\tError_Message\n")

# 查找所有以 _aligned.fa 结尾的蛋白质比对文件
protein_files = glob.glob(os.path.join(protein_msa_dir, "*_aligned.fa"))

for protein_file in protein_files:
    # 提取基因名（去掉路径和 _aligned.fa 后缀）
    gene_name = os.path.basename(protein_file).replace("_aligned.fa", "")
    
    # 对应的 CDS 文件路径
    cds_file = os.path.join(cds_dir, f"{gene_name}.fa")
    
    # 检查 CDS 文件是否存在
    if not os.path.exists(cds_file):
        error_msg = f"CDS file {cds_file} not found"
        print(f"Error for {gene_name}: {error_msg}")
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        continue
    
    # 生成密码子比对输出文件名
    codon_output = os.path.join(codon_output_dir, f"{gene_name}_codon.fas")
    
    # 运行 pal2nal.pl
    pal2nal_cmd = [
        "perl", pal2nal_script,
        protein_file, cds_file,
        "-nogap", "-output", "fasta"
    ]
    try:
        with open(codon_output, "w") as outfile:
            subprocess.run(pal2nal_cmd, stdout=outfile, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated codon alignment: {codon_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"pal2nal error: {e.stderr}"
        print(f"Error for {gene_name}: {error_msg}")
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        continue
    
    # 运行 HyPhy BUSTED 分析
    hyphy_output = os.path.join(hyphy_output_dir, f"{gene_name}.json")
    hyphy_cmd = [
        "hyphy", "busted",
        "--alignment", codon_output,
        "--tree", tree_file,
        "--output", hyphy_output,
        "--branches", "FG"
    ]
    try:
        subprocess.run(hyphy_cmd, check=True, stderr=subprocess.PIPE, text=True)
        print(f"Generated HyPhy BUSTED output: {hyphy_output}")
    except subprocess.CalledProcessError as e:
        error_msg = f"HyPhy BUSTED error: {e.stderr}"
        print(f"Error for {gene_name}: {error_msg}")
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        continue
    
    # 从 JSON 文件中提取关键结果
    try:
        with open(hyphy_output, "r") as json_file:
            data = json.load(json_file)
        
        # 提取关键字段
        p_value = data.get("test results", {}).get("p-value", "NA")
        lrt = data.get("test results", {}).get("LRT", "NA")
        log_likelihood = data.get("fits", {}).get("Unconstrained model", {}).get("logL", "NA")
        aic_c = data.get("fits", {}).get("Unconstrained model", {}).get("AIC-c", "NA")
        
        # 提取测试分支的 ω 值和占比
        test_rates = data.get("fits", {}).get("Unconstrained model", {}).get("rate distributions", {}).get("FG", [])
        test_omegas = ["NA", "NA", "NA"]
        test_proportions = ["NA", "NA", "NA"]
        for i, rate in enumerate(test_rates[:3]):  # 最多取3个类别
            test_omegas[i] = rate.get("omega", "NA")
            test_proportions[i] = rate.get("proportion", "NA")
        
        # 写入汇总结果
        with open(output_summary, "a") as summary_file:
            summary_file.write(f"{gene_name}\t{p_value}\t{lrt}\t{test_omegas[0]}\t{test_proportions[0]}\t"
                              f"{test_omegas[1]}\t{test_proportions[1]}\t{test_omegas[2]}\t{test_proportions[2]}\t"
                              f"{log_likelihood}\t{aic_c}\n")
        print(f"Extracted results for {gene_name}")
    except Exception as e:
        error_msg = f"JSON processing error: {str(e)}"
        print(f"Error for {gene_name}: {error_msg}")
        with open(error_log, "a") as error_file:
            error_file.write(f"{gene_name}\t{error_msg}\n")
        continue

# 打印报错基因列表
if os.path.exists(error_log):
    with open(error_log, "r") as error_file:
        errors = error_file.readlines()[1:]  # 跳过标题行
        if errors:
            print("\nGenes with errors:")
            for error in errors:
                gene, msg = error.strip().split("\t", 1)
                print(f"  {gene}: {msg}")
        else:
            print("\nNo errors encountered.")
else:
    print("\nNo errors encountered.")

print(f"Summary results written to {output_summary}")
print(f"Error log written to {error_log}")
