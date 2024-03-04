#python3
def process_fasta(input_file, output_file):
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # 处理序列名称，去掉第一个_及其后面所有内容
                seq_name = line.split('_')[0]
                f_out.write(seq_name + '\n')
            else:
                f_out.write(line)

input_file = 'all.aln.fa'
output_file = 'output.fa'
process_fasta(input_file, output_file)
print("处理完成！输出文件：" + output_file)
