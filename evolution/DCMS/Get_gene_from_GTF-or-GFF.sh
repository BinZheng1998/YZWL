import pandas as pd
import pybedtools
import sys
import os
import re

def detect_separator(file_path):
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
            return '\t'

def extract_gene_name(attribute_str):
    # 尝试从GFF格式中提取gene名 (gene=XXX 或 Name=XXX)
    match = re.search(r'(gene=|Name=)([^;]+)', attribute_str)
    if match:
        return match.group(2).strip()

    # 尝试从GTF格式中提取 gene_name 或 gene_id
    match = re.search(r'gene_name\s+"([^"]+)"', attribute_str)
    if match:
        return match.group(1).strip()
    match = re.search(r'gene_id\s+"([^"]+)"', attribute_str)
    if match:
        return match.group(1).strip()
    
    return None  # 未找到基因名

def main(region_file, gff_file, extension, output_file):
    if not os.path.isdir(os.path.dirname(output_file)):
        print(f"Output directory does not exist: {os.path.dirname(output_file)}")
        sys.exit(1)

    sep = detect_separator(region_file)

    try:
        regions_df = pd.read_csv(region_file, sep=sep, header=None, names=['chromosome', 'start', 'end'],
                                 converters={'start': lambda x: int(x.strip('$')), 'end': lambda x: int(x.strip('$'))})
    except Exception as e:
        print(f"Error reading region file: {e}")
        sys.exit(1)

    if extension.endswith('kb'):
        extension_distance = int(extension[:-2]) * 1000
    elif extension.endswith('Mb'):
        extension_distance = int(extension[:-2]) * 1000000
    else:
        raise ValueError("Extension should be in the format of '200kb' or '1Mb'")

    results = []

    for idx, row in regions_df.iterrows():
        chromosome, start, end = row['chromosome'], row['start'], row['end']
        region_start = max(0, start - extension_distance)
        region_end = end + extension_distance

        region = pybedtools.BedTool(f'{chromosome}\t{region_start}\t{region_end}', from_string=True)

        try:
            gff = pybedtools.BedTool(gff_file)
            overlapping_genes = gff.intersect(region, wa=True)
        except Exception as e:
            print(f"Error intersecting with annotation file: {e}")
            continue

        gene_names = set()
        for feature in overlapping_genes:
            if len(feature.fields) >= 9:
                gene_name = extract_gene_name(feature.fields[8])
                if gene_name:
                    gene_names.add(gene_name)

        results.append({
            'chromosome': chromosome,
            'start': start,
            'end': end,
            'extension': extension,
            'genes': ', '.join(sorted(gene_names)) if gene_names else "None"
        })

    output_df = pd.DataFrame(results)
    output_df.to_csv(output_file, sep='\t', index=False, header=True)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_gene.py <input.txt> <input.gff/gtf> <extension> <output_gene.txt>")
        sys.exit(1)

    region_file = sys.argv[1]
    gff_file = sys.argv[2]
    extension = sys.argv[3]
    output_file = sys.argv[4]

    main(region_file, gff_file, extension, output_file)
