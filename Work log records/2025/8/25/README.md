# 20250825

#修改下载的cds文件名 
find . -name "*cds_from_genomic.fna*" -exec rename 's/cds_from_genomic\.fna/cds.fa/' {} \; 
#提取features表、蛋白质序列文件绝对路径和物种名的对应关系 
#用于提取每个基因的最长蛋白质序列并构建索引 
find "$(pwd)" -type f -name '*_feature_table.txt' | while read -r ft; do 
    dir=$(dirname "$ft") 
    protein=$(find "$dir" -type f -name '*protein.faa' | head -n 1) 
    if [ -n "$protein" ]; then 
        echo -e "$ft\t$(basename "$dir")\t$protein" 
    fi 
done 

************************************************************* 
#blastp脚本生成的结果表头可能不是\t分割而是空白分割，所有要进行替换 
sed -i 's/ /\t/g' rbh_combined2.txt  
#blastp.sh脚本生成的rbh_combined.txt文件中物种名不是\t分割的，用上述代码进行分割，不然报错 
