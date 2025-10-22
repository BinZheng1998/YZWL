#!/usr/bin/env python3
import json
import os
import glob
import argparse
from typing import Dict, Any
from pathlib import Path
from collections import Counter
from statsmodels.stats.multitest import multipletests

def extract_busted_results(data: Dict[str, Any], gene_name: str) -> Dict[str, Any]:
    """
    从 BUSTED JSON 文件中提取关键结果并科学分类选择类型。
    """
    results = {"gene_name": gene_name}
    
    if "test results" not in data:
        return results  # 不是 BUSTED 结果文件
    
    # 提取基本测试结果
    results["p_value"] = data["test results"].get("p-value", None)
    results["lrt"] = data["test results"].get("LRT", None)
    
    # 背景分支测试结果
    bg_test = data.get("test results background", {})
    results["bg_p_value"] = bg_test.get("p-value", None)
    results["bg_lrt"] = bg_test.get("LRT", None)
    
    # 共享分布测试结果
    shared_test = data.get("test results shared distributions", {})
    results["shared_p_value"] = shared_test.get("p-value", None)
    results["shared_lrt"] = shared_test.get("LRT", None)
    
    # 从 Unconstrained 模型提取 ω 分布
    unconstrained = data.get("fits", {}).get("Unconstrained model", {})
    rate_dist = unconstrained.get("Rate Distributions", {})
    
    # 前景分支 ω
    test_omega = rate_dist.get("Test", {})
    results["test_omega_2"] = test_omega.get("2", {}).get("omega", None)
    results["test_prop_2"] = test_omega.get("2", {}).get("proportion", None)
    results["test_omega_1"] = test_omega.get("1", {}).get("omega", None)
    results["test_prop_1"] = test_omega.get("1", {}).get("proportion", None)
    results["test_omega_0"] = test_omega.get("0", {}).get("omega", None)
    results["test_prop_0"] = test_omega.get("0", {}).get("proportion", None)
    
    # 背景分支 ω
    bg_omega = rate_dist.get("Background", {})
    results["bg_omega_2"] = bg_omega.get("2", {}).get("omega", None)
    results["bg_prop_2"] = bg_omega.get("2", {}).get("proportion", None)
    results["bg_omega_1"] = bg_omega.get("1", {}).get("omega", None)
    results["bg_prop_1"] = bg_omega.get("1", {}).get("proportion", None)
    results["bg_omega_0"] = bg_omega.get("0", {}).get("omega", None)
    results["bg_prop_0"] = bg_omega.get("0", {}).get("proportion", None)
    
    # BUSTED 选择类型分类
    p_value = results.get("p_value")
    test_omega_2 = results.get("test_omega_2")
    test_prop_2 = results.get("test_prop_2")
    bg_omega_1 = results.get("bg_omega_1")
    bg_prop_1 = results.get("bg_prop_1")
    
    selection_type = "Unknown"
    selection_strength = "Unknown"
    selection_notes = "No_data"
    
    if p_value is None or test_omega_2 is None:
        selection_notes = "Incomplete_data"
    else:
        # 正选择
        if p_value < 0.05 and test_omega_2 > 1:
            selection_type = "Positive_selection"
            if test_omega_2 > 10:
                selection_strength = "Very_strong"
            elif test_omega_2 > 5:
                selection_strength = "Strong"
            elif test_omega_2 > 2:
                selection_strength = "Moderate"
            else:
                selection_strength = "Weak"
            if test_prop_2 is not None:
                if test_prop_2 > 0.3:
                    selection_notes = f"Widespread_selection({test_prop_2:.3f})"
                elif test_prop_2 > 0.1:
                    selection_notes = f"Moderate_proportion({test_prop_2:.3f})"
                else:
                    selection_notes = f"Limited_sites({test_prop_2:.3f})"
            else:
                selection_notes = "Positive_detected_no_proportion"
        else:
            # 无显著正选择，检查背景 ω
            if bg_omega_1 is not None and bg_prop_1 is not None and bg_prop_1 > 0.5:
                if bg_omega_1 < 0.3:
                    selection_type = "Strong_purifying"
                    selection_strength = "Very_constrained"
                    selection_notes = f"Background_omega={bg_omega_1:.3f}"
                elif bg_omega_1 < 0.5:
                    selection_type = "Moderate_purifying"
                    selection_strength = "Constrained"
                    selection_notes = f"Background_omega={bg_omega_1:.3f}"
                elif bg_omega_1 < 0.7:
                    selection_type = "Weak_purifying"
                    selection_strength = "Mildly_constrained"
                    selection_notes = f"Background_omega={bg_omega_1:.3f}"
                elif 0.7 <= bg_omega_1 <= 1.3:
                    selection_type = "Neutral_evolution"
                    selection_strength = "Neutral"
                    selection_notes = f"Background_omega={bg_omega_1:.3f}"
                else:
                    selection_type = "Elevated_evolution"
                    selection_strength = "Accelerated_no_significance"
                    selection_notes = f"Background_omega={bg_omega_1:.3f}_no_sig_positive"
            elif test_omega_2 > 1 and p_value >= 0.05:
                selection_type = "Suggestive_positive"
                selection_strength = "Marginal"
                selection_notes = f"omega={test_omega_2:.3f}_p={p_value:.4f}(not_sig)"
            else:
                selection_type = "No_positive_detected"
                selection_notes = f"p_value={p_value:.4f}_test_omega={test_omega_2:.3f}"
    
    results["busted_selection_type"] = selection_type
    results["busted_selection_strength"] = selection_strength
    results["busted_selection_notes"] = selection_notes
    
    return results

def extract_relax_results(data: Dict[str, Any], gene_name: str) -> Dict[str, Any]:
    """
    从 RELAX JSON 文件中提取关键结果并确定选择类型。
    """
    results = {
        "gene_name": gene_name,
        "relax_selection_type": "Unknown",
        "relax_selection_notes": "No_data"
    }
    
    # 提取测试结果
    test_results = data.get('test results', {})
    results["relax_lrt"] = test_results.get('LRT', None)
    results["relax_p_value"] = test_results.get('p-value', None)
    results["K"] = test_results.get('relaxation or intensification parameter', None)
    
    # 提取测试分支 ω 分布
    fits = data.get('fits', {})
    partitioned_descriptive = fits.get('RELAX partitioned descriptive', {})
    test_dist = partitioned_descriptive.get('Rate Distributions', {}).get('Test', {})
    
    omega_values = []
    proportions = []
    for i in range(3):
        category = test_dist.get(str(i), {})
        omega = category.get('omega', None)
        proportion = category.get('proportion', None)
        omega_values.append(omega)
        proportions.append(proportion)
    
    results["relax_omega_0"] = omega_values[0]
    results["relax_prop_0"] = proportions[0]
    results["relax_omega_1"] = omega_values[1]
    results["relax_prop_1"] = proportions[1]
    results["relax_omega_2"] = omega_values[2]
    results["relax_prop_2"] = proportions[2]
    
    # RELAX 选择类型分类
    k_value = results.get("K")
    p_value = results.get("relax_p_value")
    
    if k_value is not None and p_value is not None:
        if p_value < 0.05:
            if k_value < 1:
                results["relax_selection_type"] = "Relaxed_selection"
                results["relax_selection_notes"] = f"K={k_value:.3f}_p={p_value:.4f}"
            elif k_value > 1:
                results["relax_selection_type"] = "Intensified_selection"
                results["relax_selection_notes"] = f"K={k_value:.3f}_p={p_value:.4f}"
        else:
            results["relax_selection_type"] = "No_significant_change"
            results["relax_selection_notes"] = f"K={k_value:.3f}_p={p_value:.4f}"
    else:
        results["relax_selection_notes"] = "Incomplete_data"
    
    return results

def combine_selection_types(busted_result: Dict[str, Any], relax_result: Dict[str, Any]) -> Dict[str, Any]:
    """
    合并 BUSTED 和 RELAX 的选择类型，生成最终分类。
    """
    combined_type = "Unknown"
    combined_notes = ""
    
    busted_type = busted_result.get("busted_selection_type", "Unknown")
    relax_type = relax_result.get("relax_selection_type", "Unknown")
    busted_notes = busted_result.get("busted_selection_notes", "No_data")
    relax_notes = relax_result.get("relax_selection_notes", "No_data")
    
    if busted_type == "Positive_selection":
        combined_type = "Positive_selection"
        combined_notes = f"BUSTED_{busted_notes}"
    elif relax_type == "Relaxed_selection":
        combined_type = "Relaxed_selection"
        combined_notes = f"RELAX_{relax_notes}"
    elif relax_type == "Intensified_selection":
        combined_type = "Intensified_selection"
        combined_notes = f"RELAX_{relax_notes}"
    elif busted_type in ["Strong_purifying", "Moderate_purifying", "Weak_purifying"]:
        combined_type = busted_type
        combined_notes = f"BUSTED_{busted_notes}"
    elif busted_type == "Neutral_evolution" or relax_type == "No_significant_change":
        combined_type = "Neutral_evolution"
        combined_notes = f"BUSTED_{busted_notes}_RELAX_{relax_notes}"
    elif busted_type == "Suggestive_positive":
        combined_type = "Suggestive_positive"
        combined_notes = f"BUSTED_{busted_notes}"
    else:
        combined_type = "Ambiguous"
        combined_notes = f"BUSTED_{busted_notes}_RELAX_{relax_notes}"
    
    return {"combined_selection_type": combined_type, "combined_selection_notes": combined_notes}

def process_single_file(file_path: str, analysis_type: str) -> tuple[Dict[str, Any], bool]:
    """
    处理单个 JSON 文件（BUSTED 或 RELAX）。
    返回 (结果字典, 是否有效)。
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        gene_name = Path(file_path).stem
        if analysis_type == "busted":
            results = extract_busted_results(data, gene_name)
        else:  # relax
            results = extract_relax_results(data, gene_name)
        return results, bool(results)
    except json.JSONDecodeError as e:
        print(f"错误：{file_path} 的 JSON 格式无效: {e}")
        return {"gene_name": Path(file_path).stem, "relax_selection_type": "Invalid_JSON", "relax_selection_notes": f"JSON_error_{str(e)}"} if analysis_type == "relax" else {"gene_name": Path(file_path).stem}, False
    except Exception as e:
        print(f"错误：处理 {file_path} 时发生错误: {e}")
        return {"gene_name": Path(file_path).stem, "relax_selection_type": "Processing_error", "relax_selection_notes": f"Error_{str(e)}"} if analysis_type == "relax" else {"gene_name": Path(file_path).stem}, False

def main():
    parser = argparse.ArgumentParser(description="Extract and combine BUSTED and RELAX results from JSON files.")
    parser.add_argument("--busted-folder", required=True, help="Folder containing BUSTED JSON files")
    parser.add_argument("--relax-folder", required=True, help="Folder containing RELAX JSON files")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    parser.add_argument("--verbose", action="store_true", help="Print detailed processing information")
    parser.add_argument("--neutral-range", default="0.7-1.3", help="Neutral omega range (e.g., 0.7-1.3)")
    
    args = parser.parse_args()
    
    # 解析中性进化范围
    try:
        neutral_min, neutral_max = map(float, args.neutral_range.split('-'))
    except ValueError:
        print(f"Error: Invalid neutral range format '{args.neutral_range}'. Using default 0.7-1.3")
        neutral_min, neutral_max = 0.7, 1.3
    
    # 检查输入文件夹
    if not os.path.isdir(args.busted_folder):
        print(f"Error: BUSTED folder {args.busted_folder} does not exist")
        return
    if not os.path.isdir(args.relax_folder):
        print(f"Error: RELAX folder {args.relax_folder} does not exist")
        return
    
    # 初始化计数器
    processed_count = 0
    error_count = 0
    busted_results = {}
    relax_results = {}
    
    # 处理 BUSTED 文件
    busted_files = glob.glob(os.path.join(args.busted_folder, "*.json"))
    if args.verbose:
        print(f"Found {len(busted_files)} BUSTED JSON files in {args.busted_folder}")
    
    for file_path in busted_files:
        result, is_valid = process_single_file(file_path, "busted")
        if is_valid:
            busted_results[result["gene_name"]] = result
            processed_count += 1
            if args.verbose:
                print(f"Processed BUSTED {result['gene_name']}: {result.get('busted_selection_type', 'Unknown')}")
        else:
            error_count += 1
    
    # 处理 RELAX 文件
    relax_files = glob.glob(os.path.join(args.relax_folder, "*.json"))
    if args.verbose:
        print(f"Found {len(relax_files)} RELAX JSON files in {args.relax_folder}")
    
    for file_path in relax_files:
        result, is_valid = process_single_file(file_path, "relax")
        if is_valid:
            relax_results[result["gene_name"]] = result
            processed_count += 1
            if args.verbose:
                print(f"Processed RELAX {result['gene_name']}: {result.get('relax_selection_type', 'Unknown')}")
        else:
            relax_results[result["gene_name"]] = result  # 保留无效 JSON 的基本信息
            error_count += 1
    
    # 合并结果
    all_genes = set(busted_results.keys()) | set(relax_results.keys())
    summary = []
    busted_p_values = []
    relax_p_values = []
    
    for gene in all_genes:
        busted_result = busted_results.get(gene, {"gene_name": gene})
        relax_result = relax_results.get(gene, {"gene_name": gene, "relax_selection_type": "Missing", "relax_selection_notes": "No_RELAX_data"})
        
        # 合并选择类型
        combined = combine_selection_types(busted_result, relax_result)
        
        # 收集 p 值用于 FDR 校正
        if busted_result.get("p_value") is not None:
            busted_p_values.append((gene, float(busted_result["p_value"])))
        if relax_result.get("relax_p_value") is not None:
            relax_p_values.append((gene, float(relax_result["relax_p_value"])))
        
        # 创建输出行
        row = [
            gene,
            str(busted_result.get("p_value", "NA")),
            str(busted_result.get("lrt", "NA")),
            str(busted_result.get("bg_p_value", "NA")),
            str(busted_result.get("bg_lrt", "NA")),
            str(busted_result.get("shared_p_value", "NA")),
            str(busted_result.get("shared_lrt", "NA")),
            str(busted_result.get("test_omega_2", "NA")),
            str(busted_result.get("test_prop_2", "NA")),
            str(busted_result.get("test_omega_1", "NA")),
            str(busted_result.get("test_prop_1", "NA")),
            str(busted_result.get("test_omega_0", "NA")),
            str(busted_result.get("test_prop_0", "NA")),
            str(busted_result.get("bg_omega_2", "NA")),
            str(busted_result.get("bg_prop_2", "NA")),
            str(busted_result.get("bg_omega_1", "NA")),
            str(busted_result.get("bg_prop_1", "NA")),
            str(busted_result.get("bg_omega_0", "NA")),
            str(busted_result.get("bg_prop_0", "NA")),
            str(busted_result.get("busted_selection_type", "NA")),
            str(busted_result.get("busted_selection_strength", "NA")),
            str(busted_result.get("busted_selection_notes", "NA")),
            str(relax_result.get("relax_lrt", "NA")),
            str(relax_result.get("relax_p_value", "NA")),
            str(relax_result.get("K", "NA")),
            str(relax_result.get("relax_omega_0", "NA")),
            str(relax_result.get("relax_prop_0", "NA")),
            str(relax_result.get("relax_omega_1", "NA")),
            str(relax_result.get("relax_prop_1", "NA")),
            str(relax_result.get("relax_omega_2", "NA")),
            str(relax_result.get("relax_prop_2", "NA")),
            str(relax_result.get("relax_selection_type", "NA")),
            str(relax_result.get("relax_selection_notes", "NA")),
            str(combined["combined_selection_type"]),
            str(combined["combined_selection_notes"])
        ]
        summary.append(row)
    
    # FDR 校正
    header = [
        "Gene_Name", "BUSTED_p_value", "BUSTED_LRT", "BUSTED_bg_p_value", "BUSTED_bg_LRT",
        "BUSTED_shared_p_value", "BUSTED_shared_LRT",
        "BUSTED_test_omega_2", "BUSTED_test_prop_2", "BUSTED_test_omega_1", "BUSTED_test_prop_1",
        "BUSTED_test_omega_0", "BUSTED_test_prop_0",
        "BUSTED_bg_omega_2", "BUSTED_bg_prop_2", "BUSTED_bg_omega_1", "BUSTED_bg_prop_1",
        "BUSTED_bg_omega_0", "BUSTED_bg_prop_0",
        "BUSTED_Selection_Type", "BUSTED_Selection_Strength", "BUSTED_Selection_Notes",
        "RELAX_LRT", "RELAX_p_value", "RELAX_K",
        "RELAX_omega_0", "RELAX_prop_0", "RELAX_omega_1", "RELAX_prop_1", "RELAX_omega_2", "RELAX_prop_2",
        "RELAX_Selection_Type", "RELAX_Selection_Notes",
        "Combined_Selection_Type", "Combined_Selection_Notes"
    ]
    
    if busted_p_values:
        genes, p_vals = zip(*busted_p_values)
        _, busted_fdr, _, _ = multipletests(p_vals, method='fdr_bh')
        busted_fdr_dict = dict(zip(genes, busted_fdr))
        header.append("BUSTED_FDR_p_value")
    
    if relax_p_values:
        genes, p_vals = zip(*relax_p_values)
        _, relax_fdr, _, _ = multipletests(p_vals, method='fdr_bh')
        relax_fdr_dict = dict(zip(genes, relax_fdr))
        header.append("RELAX_FDR_p_value")
    
    # 写入输出文件
    with open(args.output, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in summary:
            if "BUSTED_FDR_p_value" in header:
                row.append(str(busted_fdr_dict.get(row[0], "NA")))
            if "RELAX_FDR_p_value" in header:
                row.append(str(relax_fdr_dict.get(row[0], "NA")))
            f.write("\t".join(row) + "\n")
    
    # 打印总结
    print(f"\n=== BUSTED and RELAX Results Extraction Complete ===")
    print(f"BUSTED folder: {args.busted_folder}")
    print(f"RELAX folder: {args.relax_folder}")
    print(f"Output file: {args.output}")
    print(f"Total files processed: {len(busted_files) + len(relax_files)}")
    print(f"Successful analyses: {processed_count}")
    print(f"Errors/skipped files: {error_count}")
    
    if summary:
        selection_types = Counter([row[header.index("Combined_Selection_Type")] for row in summary])
        print(f"\nCombined selection type distribution:")
        for sel_type, count in selection_types.most_common():
            print(f"  {sel_type}: {count} genes")

if __name__ == "__main__":
    main()
