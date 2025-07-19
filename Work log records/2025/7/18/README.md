# 1.rnaseq_bam_reads_counts.py  
脚本是用来从bam文件中提取每个区间的reads数量的。输入文件是多个bam文件绝对路径。
# 2.rnaseq_bam_reads_count_merge.sh  
用来合并多个bam文件reads的结果。  
# 3.snp_stats.py
用来统计vcf文件固定窗口内的snp数量
# 4.indel_stats.py
用来统计vcf文件固定窗口内的indel数量
# 5.repeat_GC_stats.py
用来统计基因组fasta文件中窗口内GC含量和重复序列含量，注意统计重复序列需要fasta文件的重复序列是小写字母也就是repeatmasker掩盖过的结果。
