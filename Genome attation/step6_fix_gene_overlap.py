import pandas as pd
import sys

def parse_gff(file):
    return pd.read_csv(file, sep='\t', header=None, comment='#', names=[
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'
    ])

def process_genes(df):
    df.sort_values(['seqname', 'strand', 'start'], inplace=True)
    result_df = df.copy()

    # Initialize variables to track the current cluster of overlapping genes
    current_cluster = []
    current_cluster_ids = []
    parent_updates = {}

    for index, gene in df[df['feature'] == 'gene'].iterrows():
        gene_id = get_id(gene['attribute'])
        if not current_cluster or (gene['seqname'] == current_cluster[-1]['seqname'] and 
                                   gene['strand'] == current_cluster[-1]['strand'] and 
                                   gene['start'] <= current_cluster[-1]['end']):
            # Overlap detected or starting a new cluster
            current_cluster.append(gene)
            current_cluster_ids.append(gene_id)
            if len(current_cluster) > 1:
                # Extend the end position for the first gene in the cluster
                current_cluster[0]['end'] = max(current_cluster[0]['end'], gene['end'])
                # Update the gene ID
                merged_id = "_".join(current_cluster_ids)
                result_df.at[current_cluster[0].name, 'end'] = current_cluster[0]['end']
                result_df.at[current_cluster[0].name, 'attribute'] = update_id(current_cluster[0]['attribute'], merged_id)
                for gid in current_cluster_ids:
                    parent_updates[gid] = merged_id
                result_df.drop(index, inplace=True)
        else:
            # Process the previous cluster and start a new one
            current_cluster = [gene]
            current_cluster_ids = [gene_id]

    # Ensure the last cluster is processed
    if len(current_cluster) > 1:
        merged_id = "_".join(current_cluster_ids)
        result_df.at[current_cluster[0].name, 'end'] = current_cluster[0]['end']
        result_df.at[current_cluster[0].name, 'attribute'] = update_id(current_cluster[0]['attribute'], merged_id)
        for gid in current_cluster_ids:
            parent_updates[gid] = merged_id

    # Update Parent attributes in RNA features
    def update_parent_attributes(row):
        if row['feature'] != 'gene' and 'Parent=' in row['attribute']:
            attrs = parse_attributes(row['attribute'])
            if 'Parent' in attrs and attrs['Parent'] in parent_updates:
                attrs['Parent'] = parent_updates[attrs['Parent']]
            return format_attributes(attrs)
        return row['attribute']

    result_df['attribute'] = result_df.apply(update_parent_attributes, axis=1)
    return result_df

def get_id(attribute):
    return attribute.split(';')[0].split('=')[1]

def update_id(attribute, new_id):
    parts = attribute.split(';')
    for i, part in enumerate(parts):
        if part.startswith('ID='):
            parts[i] = f"ID={new_id}"
    return ';'.join(parts)

def parse_attributes(attr_string):
    return {key: value for part in attr_string.split(';') if '=' in part for key, value in [part.split('=')]}

def format_attributes(attr_dict):
    return ';'.join([f'{key}={value}' for key, value in attr_dict.items()])

def main(input_file, output_file):
    df = parse_gff(input_file)
    updated_df = process_genes(df)
    updated_df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test.py input.gff output.gff")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
