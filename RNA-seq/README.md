# RNA-seq
## Step1 
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
## Step2
