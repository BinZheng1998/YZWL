import sys
from Bio import SeqIO
import re

def parse_attributes(attributes, file_type):
    """
    Parse GFF or GTF attribute column to extract gene ID and other attributes.

    Args:
        attributes (str): 9th column attribute string
        file_type (str): 'gff' or 'gtf'

    Returns:
        dict: Parsed attribute key-value pairs
    """
    attr_dict = {}
    try:
        if file_type == 'gff':
            for item in attributes.split(';'):
                item = item.strip()
                if not item:
                    continue
                if '=' in item:
                    key, value = item.split('=', 1)
                    attr_dict[key.strip()] = value.strip()
                else:
                    attr_dict[item] = ""
        elif file_type == 'gtf':
            parts = attributes.strip().split(';')
            for part in parts:
                part = part.strip()
                if not part:
                    continue
                match = re.match(r'(\S+)\s+("?[^"]+"?|\S+)', part)
                if match:
                    key, value = match.groups()
                    value = value.strip().strip('"')
                    attr_dict[key] = value
                else:
                    attr_dict[part] = ""
    except Exception as e:
        print(f"警告: 属性解析错误: {attributes[:50]}... ({e})")
    return attr_dict

def main(gene_ids_file, annotation_file, genome_file, output_file, species):
    # Read gene IDs, stripping quotes and ignoring comments
    gene_ids = set()
    try:
        with open(gene_ids_file, 'r') as f:
            for line in f:
                stripped = line.strip()
                if stripped and not stripped.startswith('#'):
                    gene_id = stripped.split()[0].strip('"').strip("'")
                    gene_ids.add(gene_id)
        print(f"读取到 {len(gene_ids)} 个唯一基因 ID")
        print(f"前5个基因 ID: {list(gene_ids)[:5]}")
    except Exception as e:
        print(f"错误: 无法读取基因 ID 文件 {gene_ids_file}: {e}")
        sys.exit(1)

    # Determine file type
    file_type = 'gff' if annotation_file.lower().endswith(('.gff', '.gff3')) else 'gtf'

    # Parse annotation file
    gene_locations = {}
    genes_found = set()
    attribute_samples = []
    try:
        with open(annotation_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                seqid = parts[0]
                feature_type = parts[2]
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]
                except (IndexError, ValueError):
                    print(f"警告: 无效的 GTF 行: {line[:50]}...")
                    continue

                if feature_type.lower() not in ['gene', 'pseudogene']:
                    continue

                if len(attribute_samples) < 5 and file_type == 'gtf':
                    attribute_samples.append(attributes)

                attr_dict = parse_attributes(attributes, file_type)

                gene_id = None
                if file_type == 'gff':
                    gene_id = attr_dict.get('gene')
                elif file_type == 'gtf':
                    gene_id = (attr_dict.get('gene_id') or
                               attr_dict.get('gene_name') or
                               attr_dict.get('Name') or
                               attr_dict.get('ID'))

                if gene_id and gene_id in gene_ids:
                    if start > end:
                        print(f"警告: 基因 {gene_id} 坐标无效 (start={start} > end={end}), 跳过")
                        continue
                    gene_name = attr_dict.get('gene_name', gene_id)
                    genes_found.add(gene_id)
                    gene_locations[gene_id] = (seqid, start, end, strand, gene_name)
    except Exception as e:
        print(f"错误: 无法解析注释文件 {annotation_file}: {e}")
        sys.exit(1)

    if file_type == 'gtf' and attribute_samples:
        print("\n调试: GTF 文件中前5个基因特征的属性列示例:")
        for i, attr in enumerate(attribute_samples, 1):
            print(f"样本 {i}: {attr}")

    not_found = gene_ids - genes_found
    if not_found:
        print(f"\n警告: {len(not_found)} 个基因在注释文件中未找到")
        print("前10个未找到的基因:", list(not_found)[:10])

    try:
        genome_index = SeqIO.index(genome_file, "fasta")
        print(f"索引了 {len(genome_index)} 条染色体序列")
    except Exception as e:
        print(f"错误: 无法索引基因组文件 {genome_file}: {e}")
        sys.exit(1)

    extracted_count = 0
    try:
        with open(output_file, 'w') as out_f:
            for gene_id in sorted(gene_ids):
                if gene_id in gene_locations:
                    seqid, start, end, strand, gene_name = gene_locations[gene_id]
                    if seqid in genome_index:
                        record = genome_index[seqid]
                        if start < 1:
                            print(f"警告: 基因 {gene_id} 起始位置调整为 1")
                            start = 1
                        if end > len(record):
                            print(f"警告: 基因 {gene_id} 终止位置调整为染色体末端 ({len(record)})")
                            end = len(record)
                        sequence = record.seq[start-1:end]
                        if strand == '-':
                            sequence = sequence.reverse_complement()
                        # Modified: Updated header format to gene_id-species
                        header = f">{gene_id}-{species} | {seqid}:{start}-{end}({strand})"
                        out_f.write(f"{header}\n")
                        for i in range(0, len(sequence), 80):
                            out_f.write(f"{sequence[i:i+80]}\n")
                        extracted_count += 1
                    else:
                        print(f"错误: 基因 {gene_id} 的染色体 {seqid} 在基因组文件中未找到")
    except Exception as e:
        print(f"错误: 无法写入输出文件 {output_file}: {e}")
        sys.exit(1)
    finally:
        genome_index.close()

    print(f"成功提取 {extracted_count}/{len(gene_ids)} 个基因序列到 {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("用法: python gene_extractor.py <gene_ids.txt> <annotation.gff/gtf> <genome.fa> <output.fa> <species>")
        print("示例: python gene_extractor.py genes.txt annotation.gtf genome.fasta output_genes.fa Gallus_gallus")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
