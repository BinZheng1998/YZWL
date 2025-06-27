# RERconverge
## Step1 
下载蛋白质序列文件和注释结果，比如pig：https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/ 
GCF_000003025.6_Sscrofa11.1_protein.faa.gz  
GCF_000003025.6_Sscrofa11.1_feature_table.txt.gz  
使用get_largest_protein.py从faa蛋白质序列文件中提取每个基因的最长蛋白序列作为代表。
```
Usage
get_protein.py [-h] --input-feature-txt INPUT_FEATURE_TXT --species SPECIES --input-protein-faa INPUT_PROTEIN_FAA
```
