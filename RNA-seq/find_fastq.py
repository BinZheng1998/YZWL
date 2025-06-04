import os

def find_fastq_files():
    current_dir = os.getcwd()
    
    single_end_files = []
    paired_end_files = {}
    
    for root, dirs, files in os.walk(current_dir):
        for file in files:
            if file.endswith(".fastq.gz"):
                file_path = os.path.join(root, file)
                if file.endswith("_1.fastq.gz"):
                    base_name = file.rsplit("_1.fastq.gz", 1)[0]
                    if base_name not in paired_end_files:
                        paired_end_files[base_name] = ["", ""]
                    paired_end_files[base_name][0] = file_path
                elif file.endswith("_2.fastq.gz"):
                    base_name = file.rsplit("_2.fastq.gz", 1)[0]
                    if base_name not in paired_end_files:
                        paired_end_files[base_name] = ["", ""]
                    paired_end_files[base_name][1] = file_path
                elif file.endswith(".fastq.gz"):
                    base_name = file.rsplit(".fastq.gz", 1)[0]
                    single_end_files.append(file_path)
    
    output_file = "fastq_files_list.txt"
    with open(output_file, 'w') as f:
        for file in single_end_files:
            f.write(f"{file}\n")
        
        for base_name, files in paired_end_files.items():
            if files[0] and files[1]:  
                f.write(f"{files[0]}\t{files[1]}\n")
            elif files[0]: 
                f.write(f"{files[0]}\n")
            elif files[1]:  
                f.write(f"{files[1]}\n")
    
    print(f"文件列表已保存到 {output_file}")

if __name__ == "__main__":
    find_fastq_files()
