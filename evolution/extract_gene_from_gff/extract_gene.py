import pandas as pd
import pybedtools
import sys
import os

def detect_separator(file_path):
    # 尝试读取文件的前几行以检测分隔符
    with open(file_path, 'r') as f:
        first_line = f.readline()
        if '^I' in first_line:
            return '^I'
        elif '\t' in first_line:
            return '\t'
        elif ',' in first_line:
            return ','
        elif ' ' in first_line:
            return ' '
        else:
            print("Warning: Unrecognized separator. Using tab as default.")
            return '\t'  # 默认使用制表符

def main(region_file, gff_file, extension, output_file):
    # 检查输出目录
    if not os.path.isdir(os.path.dirname(output_file)):
        print(f"Output directory does not exist: {os.path.dirname(output_file)}")
        sys.exit(1)

    # 自动检测分隔符
    sep = detect_separator(region_file)

    # 读取没有表头的区间文件，并手动指定列名
    try:
        # 读取时去掉行尾的 $ 符号
        regions_df = pd.read_csv(region_file, sep=sep, header=None, names=['chromosome', 'start', 'end'], 
                                 converters={'start': lambda x: int(x.strip('$')), 'end': lambda x: int(x.strip('$'))})
    except Exception as e:
        print(f"Error reading region file: {e}")
        sys.exit(1)

    # 解析扩展区间的大小
    if extension.endswith('kb'):
        extension_distance = int(extension[:-2]) * 1000
    elif extension.endswith('Mb'):
        extension_distance = int(extension[:-2]) * 1000000
    else:
        raise ValueError("Extension should be in the format of '200kb' or '1Mb'")

    # 存储结果
    results = []

    # 处理每个区间
    for idx, row in regions_df.iterrows():
        chromosome, start, end = row['chromosome'], row['start'], row['end']
        
        # 扩展上下游区间
        region_start = max(0, start - extension_distance)
        region_end = end + extension_distance

        # 创建bedtools对象
        region = pybedtools.BedTool(f'{chromosome}\t{region_start}\t{region_end}', from_string=True)

        try:
            gff = pybedtools.BedTool(gff_file)
            overlapping_genes = gff.intersect(region, wa=True)
        except Exception as e:
            print(f"Error intersecting with GFF file: {e}")
            continue

        # 提取基因名
        gene_names = set()
        for feature in overlapping_genes:
            attributes = feature.fields[8]
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=')
                    if key.strip() == 'gene':
                        gene_names.add(value.strip())
                        break
        
        # 记录结果
        results.append({
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'extension': extension,
            'genes': ', '.join(gene_names)
        })

    # 将结果写入输出文件
    output_df = pd.DataFrame(results)
    output_df.to_csv(output_file, sep='\t', index=False, header=True)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extrct_gene.py <input.txt> <input.gff> <extension> <output_gene.txt>")
        sys.exit(1)

    region_file = sys.argv[1]
    gff_file = sys.argv[2]
    extension = sys.argv[3]
    output_file = sys.argv[4]

    main(region_file, gff_file, extension, output_file)
