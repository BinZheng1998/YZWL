import argparse
import subprocess
import sys
from collections import defaultdict
import re

from Bio import SeqIO

def parse_gtf(gtf_file):
    """
    解析GTF文件，提取基因、转录本和CDS信息。
    返回每个基因的转录本列表，每个转录本包含transcript_id, protein_id (if available), cds_length。
    """
    genes = defaultdict(list)  # gene_id -> list of (transcript_id, protein_id, cds_length)
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, phase, attributes = fields
            
            # 解析attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if ' ' in attr:
                    key, value = attr.split(' ', 1)
                    attr_dict[key] = value.strip('"')
            
            if feature in ['transcript', 'mRNA']:
                transcript_id = attr_dict.get('transcript_id', '')
                gene_id = attr_dict.get('gene_id', '')
                if gene_id and transcript_id:
                    # 初始化CDS长度为0
                    genes[gene_id].append((transcript_id, None, 0))  # protein_id later
            elif feature == 'CDS':
                transcript_id = attr_dict.get('transcript_id', '')
                gene_id = attr_dict.get('gene_id', '')
                protein_id = attr_dict.get('protein_id', None)
                if transcript_id and gene_id:
                    # 找到对应的转录本，更新CDS长度和protein_id
                    for i in range(len(genes[gene_id])):
                        tid, pid, length = genes[gene_id][i]
                        if tid == transcript_id:
                            # CDS长度：累加每个CDS的长度（假设无重叠）
                            new_length = length + (int(end) - int(start) + 1)
                            new_pid = protein_id if protein_id else pid
                            genes[gene_id][i] = (tid, new_pid, new_length)
                            break
    
    # 清理：只保留有CDS长度的转录本
    cleaned_genes = {}
    for gene_id, trans in genes.items():
        cleaned_trans = [(tid, pid, length) for tid, pid, length in trans if length > 0]
        if cleaned_trans:
            cleaned_genes[gene_id] = cleaned_trans
    
    return cleaned_genes

def select_longest_transcript(genes):
    """
    对于每个基因，选择最长转录本（基于CDS长度）。如果长度相同，保留第一个。
    返回 dict: gene_id -> (transcript_id, protein_id2: real protein_id or modified transcript_id)
    """
    selected = {}
    for gene_id, trans_list in genes.items():
        if not trans_list:
            continue
        # max with key returns the first occurrence of the max value
        longest = max(trans_list, key=lambda x: x[2])
        transcript_id, protein_id, _ = longest
        if protein_id:
            prot_id2 = protein_id
        else:
            # 用transcript_id代替，并去掉rna-前缀
            prot_id2 = transcript_id.replace("rna-", "")
        selected[gene_id] = (transcript_id, prot_id2)
    
    return selected

def extract_protein_sequences(selected, protein_fa_file, output_fa_file, species):
    """
    从蛋白质FASTA中提取选定转录本对应的序列。
    假设FASTA header 以 prot_id2 或 transcript_id 开头（取第一个部分作为ID）。
    如果prot_id2重复，选择序列最长的那个（长度相同时取第一个）。
    生成新FASTA和TXT文件。FASTA sequence ID 为 species_protein_{i+1}，后跟空格分隔的描述。
    """
    # 解析FASTA，记录 id -> seq
    fa_dict = {}
    for record in SeqIO.parse(protein_fa_file, 'fasta'):
        header = str(record.id).split()[0]  # 取第一个部分作为ID
        fa_dict[header] = str(record.seq)
    
    # 收集 per prot_id2: list of (seq_len, transcript_id, gene_id, seq)
    prot_to_candidates = defaultdict(list)
    
    for gene_id, (transcript_id, prot_id2) in selected.items():
        # 尝试匹配 prot_id2 或 transcript_id
        seq = None
        if prot_id2 in fa_dict:
            seq = fa_dict[prot_id2]
        elif transcript_id in fa_dict:
            seq = fa_dict[transcript_id]
        
        if seq:
            seq_len = len(seq)
            prot_to_candidates[prot_id2].append((seq_len, transcript_id, gene_id, seq))
    
    # 对于每个prot_id2，选择最长的
    final_entries = []
    for prot_id2, candidates in prot_to_candidates.items():
        if candidates:
            # 选择seq_len最大的，tie时第一个
            best = max(candidates, key=lambda x: x[0])
            _, transcript_id, gene_id, seq = best
            final_entries.append((prot_id2, transcript_id, gene_id, seq))
    
    # 收集数据用于TXT，并分配unique_id
    txt_data = []
    index = 1
    with open(output_fa_file, 'w') as out_fa:
        for prot_id2, transcript_id, gene_id, seq in final_entries:
            unique_id = f"{species}_protein_{index}"
            index += 1
            # 写FASTA: >unique_id prot_id2 transcript_id gene_id species (空格分隔)
            fa_header = f">{unique_id} {prot_id2} {transcript_id} {gene_id} {species}"
            out_fa.write(f"{fa_header}\n{seq}\n")
            
            # TXT行
            txt_data.append(f"{unique_id}\t{prot_id2}\t{transcript_id}\t{species}\t{gene_id}")
    
    # 写TXT
    txt_file = output_fa_file.replace('.fa', '_info.txt')
    with open(txt_file, 'w') as out_txt:
        out_txt.write("protein_id1\tprotein_id2\ttranscript_id\tspecies\tgene_id\n")
        for line in txt_data:
            out_txt.write(line + "\n")
    
    print(f"Generated {output_fa_file} and {txt_file}")

def build_blastdb(output_fa_file, out_dir):
    """
    使用makeblastdb构建数据库。
    假设makeblastdb在PATH中。
    """
    cmd = [
        'makeblastdb',
        '-in', output_fa_file,
        '-dbtype', 'prot',
        '-out', out_dir
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"Built BLAST database in {out_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error building BLAST db: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Extract longest transcript proteins and build BLAST db')
    parser.add_argument('--gtf', required=True, help='Input GTF file')
    parser.add_argument('--species', required=True, help='Species name')
    parser.add_argument('--protein', dest='protein_fa', required=True, help='Input protein FASTA file')
    parser.add_argument('-out', required=True, help='Output FASTA file (e.g., human.fa)')
    
    args = parser.parse_args()
    
    # 解析GTF
    genes = parse_gtf(args.gtf)
    
    # 选择最长转录本
    selected = select_longest_transcript(genes)
    
    if not selected:
        print("No transcripts found!", file=sys.stderr)
        sys.exit(1)
    
    # 提取序列
    extract_protein_sequences(selected, args.protein_fa, args.out, args.species)
    
    # 构建BLAST db
    out_dir = args.out.replace('.fa', '')
    build_blastdb(args.out, out_dir)

if __name__ == '__main__':
    main()
