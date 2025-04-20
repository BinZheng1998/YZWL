import pandas as pd

def parse_attributes(attr_str):
    if pd.isna(attr_str):
        return {}
    pairs = attr_str.split(';')
    attr_dict = {}
    for pair in pairs:
        if '=' in pair:
            key, value = pair.split('=', 1)
            attr_dict[key.strip()] = value.strip()
    return attr_dict

def find_containing_exon(tid, S, E, exon_mapping):
    if tid in exon_mapping:
        for es, ee, eid in exon_mapping[tid]:
            if es <= S and E <= ee:
                return eid
    return None

# 读取 GFF 文件
gff_file = 'chicken_final.gff'  # 替换为您的 GFF 文件路径
gff_df = pd.read_csv(gff_file, sep='\t', header=None, names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
gff_df['attributes_dict'] = gff_df['attributes'].apply(parse_attributes)
gff_df['ID'] = gff_df['attributes_dict'].apply(lambda x: x.get('ID', ''))
gff_df['Parent'] = gff_df['attributes_dict'].apply(lambda x: x.get('Parent', ''))

# 创建基因和转录本映射
genes_df = gff_df[gff_df['type'] == 'gene'][['ID', 'attributes_dict']]
genes_df['gene_biotype'] = genes_df['attributes_dict'].apply(lambda x: x.get('gene_biotype', ''))
genes_df['symbol'] = genes_df['attributes_dict'].apply(lambda x: x.get('Name', ''))
genes_df.set_index('ID', inplace=True)

transcript_types = ['mRNA', 'lnc_RNA']
transcripts_df = gff_df[gff_df['type'].isin(transcript_types)][['ID', 'Parent', 'type']]
transcripts_df.rename(columns={'ID': 'transcript_id', 'Parent': 'gene_id'}, inplace=True)
transcripts_df['is_coding'] = transcripts_df['type'] == 'mRNA'
transcripts_df.set_index('transcript_id', inplace=True)

# 创建外显子映射
exons_df = gff_df[gff_df['type'] == 'exon'][['Parent', 'start', 'end', 'ID']]
exons_df.rename(columns={'Parent': 'transcript_id'}, inplace=True)
exon_mapping = {}
for _, row in exons_df.iterrows():
    tid = row['transcript_id']
    if tid not in exon_mapping:
        exon_mapping[tid] = []
    exon_mapping[tid].append((row['start'], row['end'], row['ID']))

# 提取子特征（CDS, five_prime_UTR, three_prime_UTR）
subfeature_types = ['CDS', 'five_prime_UTR', 'three_prime_UTR']
subfeatures_df = gff_df[gff_df['type'].isin(subfeature_types)]
subfeatures_df['feature'] = subfeatures_df['type'].map({
    'CDS': 'protein_coding',
    'five_prime_UTR': 'utr5',
    'three_prime_UTR': 'utr3'
})
subfeatures_df['transcript_id'] = subfeatures_df['Parent']
subfeatures_df = subfeatures_df.merge(transcripts_df[['gene_id']], left_on='transcript_id', right_index=True)
subfeatures_df = subfeatures_df.merge(genes_df[['symbol']], left_on='gene_id', right_index=True)
subfeatures_df['exon'] = subfeatures_df.apply(lambda row: find_containing_exon(row['transcript_id'], row['start'], row['end'], exon_mapping), axis=1)

# 提取非编码外显子
all_exons_df = gff_df[gff_df['type'] == 'exon']
all_exons_df['transcript_id'] = all_exons_df['Parent']
all_exons_df = all_exons_df.merge(transcripts_df[['is_coding', 'gene_id']], left_on='transcript_id', right_index=True)
non_coding_exons_df = all_exons_df[all_exons_df['is_coding'] == False]
non_coding_exons_df = non_coding_exons_df.merge(genes_df[['gene_biotype', 'symbol']], left_on='gene_id', right_index=True)
non_coding_exons_df['feature'] = non_coding_exons_df['gene_biotype']
non_coding_exons_df['exon'] = non_coding_exons_df['ID']

# 合并结果
subfeatures_result = subfeatures_df[['seqid', 'start', 'end', 'strand', 'feature', 'gene_id', 'exon', 'transcript_id', 'symbol']]
subfeatures_result.rename(columns={'seqid': 'chromosome', 'gene_id': 'gene', 'transcript_id': 'transcript'}, inplace=True)
subfeatures_result['width'] = subfeatures_result['end'] - subfeatures_result['start'] + 1

non_coding_exons_result = non_coding_exons_df[['seqid', 'start', 'end', 'strand', 'feature', 'gene_id', 'exon', 'transcript_id', 'symbol']]
non_coding_exons_result.rename(columns={'seqid': 'chromosome', 'gene_id': 'gene', 'transcript_id': 'transcript'}, inplace=True)
non_coding_exons_result['width'] = non_coding_exons_result['end'] - non_coding_exons_result['start'] + 1

result_df = pd.concat([subfeatures_result, non_coding_exons_result], ignore_index=True)
result_df = result_df.sort_values(by=['chromosome', 'start'])

# 保存结果
output_file = 'output_for_gviz.txt'
result_df.to_csv(output_file, sep='\t', index=False, header=True)

print(f"转换完成。输出已保存至 {output_file}")
