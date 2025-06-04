# RNA-seq and fusion
## Step1
Create two conda environments, rnaseq and fusion. you can use the appropriate yaml files.
## Step2
  Use find_fastq.py find fastq.gz/fq.gz files in target folder  
  command : 
  ```
  python find_fastq.py /path/to/target_folder
  ```
will create a file(fastq_files_list.txt) like :
|fastq1 | fastq2 |
|---|---|
|/path/to/input_1.fq.gz|/path/to/input_2.fq.gz|
|/path/to/input.fq.gz||
|/path/to/input_1.fastq.gz|/path/to/input_2.fastq.gz|
|/path/to/input.fastq.gz||

## Step3
use run.py to analysis RNA-seq fastq data  
command:
```
python run.py --input-fastq-file fastq_files_list.txt --sample-numbers 8 --rnaseq yes --fusion yes
```
--input-fastq-file : pleasr see step1  
--sample-numbers : Number of samples in parallel  
--rnaseq Perform : rnaseq analysis  
--fusion Perform : fusion analysis([arriba](https://github.com/suhrig/arriba))  
