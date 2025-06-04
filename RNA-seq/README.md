# RNA-seq
# Step1 
  use find_fastq.py find fastq.gz/fq.gz files in target folder
  command : 
  ```
  python find_fastq.py /path/to/target_folder
  ```
will create a file like :
|-----:|-----------|
| /path/to/input_1.fq.gz |/path/to/input_2.fq.gz|
|/path/to/input.fq.gz||
|/path/to/input_1.fastq.gz| /path/to/input_2.fastq.gz|
|/path/to/input.fastq.gz||
# Step2
