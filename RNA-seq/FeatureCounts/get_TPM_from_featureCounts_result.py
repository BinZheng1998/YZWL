import pandas as pd
import os

with open('featureCounts_tab.txt', 'r') as f:
    file_paths = f.read().splitlines()

output_file_prefix = 'merged_gene_abund_results_part_' 
final_output_file = 'merged_gene_abund_results.csv' 
batch_size = 100  

file_count = 0 
part_count = 1  
current_batch = [] 

for file_path in file_paths:
    current_batch.append(file_path)
    file_count += 1

    if file_count == batch_size or file_path == file_paths[-1]:
        output_file = f"{output_file_prefix}{part_count}.csv"
        first_file = True 
        
        for batch_file in current_batch:
            df = pd.read_csv(batch_file, sep='\t', usecols=[0, 4, 5, 8], header=None, dtype={0: str, 4: str, 5: str})
            
            df[8] = pd.to_numeric(df[8], errors='coerce')

            sample_name = os.path.basename(batch_file).replace('_gene_abund.tab', '')
            df.columns = ['Col1', 'Col5', 'Col6', sample_name]

            if first_file:
                df.to_csv(output_file, index=False)
                first_file = False
            else:
                merged_df = pd.read_csv(output_file, dtype=str)

                merged_df = pd.merge(merged_df, df, on=['Col1', 'Col5', 'Col6'], how='outer')

                for col in merged_df.columns[3:]:
                    merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce')

                merged_df.to_csv(output_file, index=False)

        current_batch = []
        file_count = 0
        part_count += 1

first_file = True
for i in range(1, part_count):
    part_file = f"{output_file_prefix}{i}.csv"
    df = pd.read_csv(part_file, dtype=str)
    
    if first_file:
        df.to_csv(final_output_file, index=False)
        first_file = False
    else:
        merged_df = pd.read_csv(final_output_file, dtype=str)
        merged_df = pd.merge(merged_df, df, on=['Col1', 'Col5', 'Col6'], how='outer')

        for col in merged_df.columns[3:]:
            merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce')

        merged_df.to_csv(final_output_file, index=False)

print(f"All files have been successfully merged and saved to '{final_output_file}'.")
