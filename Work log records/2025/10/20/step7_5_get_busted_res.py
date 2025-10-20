import os
import json
import argparse
from pathlib import Path

def extract_relax_results(json_file):
    """
    从单个RELAX JSON文件中提取关键结果。
    返回字典，包含基因名、LRT、p值、K值等。
    """
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # 提取基因名（直接使用文件名，不包括.json后缀）
        gene_name = Path(json_file).stem
        
        # 提取test results
        test_results = data.get('test results', {})
        lrt = test_results.get('LRT', 'N/A')
        p_value = test_results.get('p-value', 'N/A')
        k_value = test_results.get('relaxation or intensification parameter', 'N/A')
        
        # 提取Test分支的ω分布（RELAX partitioned descriptive模型）
        fits = data.get('fits', {})
        partitioned_descriptive = fits.get('RELAX partitioned descriptive', {})
        test_dist = partitioned_descriptive.get('Rate Distributions', {}).get('Test', {})
        
        omega_values = []
        proportions = []
        for i in range(3):  # 假设有0,1,2三个类别
            category = test_dist.get(str(i), {})
            omega = category.get('omega', 'N/A')
            proportion = category.get('proportion', 'N/A')
            omega_values.append(omega)
            proportions.append(proportion)
        
        return {
            'gene_name': gene_name,
            'LRT': lrt,
            'p_value': p_value,
            'K': k_value,
            'omega_0': omega_values[0],
            'prop_0': proportions[0],
            'omega_1': omega_values[1],
            'prop_1': proportions[1],
            'omega_2': omega_values[2],
            'prop_2': proportions[2]
        }
    except Exception as e:
        print(f"Error processing {json_file}: {str(e)}")
        return None

def process_folder(input_folder, output_file):
    """
    处理指定文件夹中的所有.json文件，提取RELAX结果并写入输出文件。
    """
    # 确保输入文件夹存在
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder {input_folder} does not exist.")
        return
    
    # 准备输出文件
    with open(output_file, 'w') as out_f:
        # 写入表头
        header = ("Gene_Name\tLRT\tp_value\tK\t"
                  "Omega_0\tProp_0\tOmega_1\tProp_1\tOmega_2\tProp_2\n")
        out_f.write(header)
        
        # 遍历文件夹中的.json文件
        for file_name in os.listdir(input_folder):
            if file_name.endswith('.json'):
                json_file = os.path.join(input_folder, file_name)
                result = extract_relax_results(json_file)
                if result:
                    # 写入结果
                    line = (f"{result['gene_name']}\t{result['LRT']}\t{result['p_value']}\t{result['K']}\t"
                            f"{result['omega_0']}\t{result['prop_0']}\t"
                            f"{result['omega_1']}\t{result['prop_1']}\t"
                            f"{result['omega_2']}\t{result['prop_2']}\n")
                    out_f.write(line)
                    print(f"Processed: {file_name}")
    
    print(f"Results written to {output_file}")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="Extract RELAX results from JSON files.")
    parser.add_argument('--input-folder', type=str, required=True, 
                        help="Folder containing RELAX JSON files")
    parser.add_argument('--output', type=str, required=True, 
                        help="Output file to save results")
    args = parser.parse_args()
    
    # 处理文件夹
    process_folder(args.input_folder, args.output)

if __name__ == '__main__':
    main()
