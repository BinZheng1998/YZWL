import os
import glob
import pandas as pd
from collections import Counter

# For pos
files = glob.glob("*_pos_DCMS_regions.txt")
total_count = Counter()
sig_count = Counter()

for file_path in files:
    df = pd.read_csv(file_path, sep='\t', header=0, usecols=['CHROM', 'BIN_START', 'BIN_END', 'p_values', 'threshold_005'])
    df['is_significant'] = df['p_values'] < df['threshold_005']
    
    grouped = df.groupby(['CHROM', 'BIN_START', 'BIN_END'])['is_significant'].sum()
    
    # Update total_count: +1 for each unique group
    total_count.update(grouped.index)
    
    # Update sig_count: add the sum for each group
    sig_count.update(grouped.to_dict())
    
    print(f"已读取文件: {os.path.basename(file_path)}")

# Build the dataframe
pos = pd.DataFrame(list(total_count.keys()), columns=['CHROM', 'BIN_START', 'BIN_END'])
pos['total_files'] = [total_count[k] for k in total_count]
pos['sig_files'] = [sig_count[k] for k in total_count]
pos['prop'] = pos['sig_files'] / pos['total_files']
pos = pos.sort_values(['CHROM', 'BIN_START']).reset_index(drop=True)
print(pos.head(10))
pos.to_csv("dog_pairwise_pos_summary.tsv", sep='\t', index=False, na_rep='')

# For neg
files = glob.glob("*_neg_DCMS_regions.txt")
total_count = Counter()
sig_count = Counter()

for file_path in files:
    df = pd.read_csv(file_path, sep='\t', header=0, usecols=['CHROM', 'BIN_START', 'BIN_END', 'p_values', 'threshold_005'])
    df['is_significant'] = df['p_values'] < df['threshold_005']
    
    grouped = df.groupby(['CHROM', 'BIN_START', 'BIN_END'])['is_significant'].sum()
    
    # Update total_count: +1 for each unique group
    total_count.update(grouped.index)
    
    # Update sig_count: add the sum for each group
    sig_count.update(grouped.to_dict())
    
    print(f"已读取文件: {os.path.basename(file_path)}")

# Build the dataframe
neg = pd.DataFrame(list(total_count.keys()), columns=['CHROM', 'BIN_START', 'BIN_END'])
neg['total_files'] = [total_count[k] for k in total_count]
neg['sig_files'] = [sig_count[k] for k in total_count]
neg['prop'] = neg['sig_files'] / neg['total_files']
neg = neg.sort_values(['CHROM', 'BIN_START']).reset_index(drop=True)
print(neg.head(10))
neg.to_csv("dog_pairwise_neg_summary.tsv", sep='\t', index=False, na_rep='')
