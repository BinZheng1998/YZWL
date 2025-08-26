#首先使用AMAS程序合并trimal修剪后的氨基酸比对结果，如下：
python /home/bzheng/.conda/envs/evolution/lib/python3.12/site-packages/amas/AMAS.py concat -f fasta -d aa -i *_trimmed.fa --part-format nexus
#生成concatenated.out和partitions.txt两个结果，其中concatenated.out是我们需要的合并后的结果
/home/bzheng/project/11_multic_species_homologs/tree/concatenate

#iqtree2构建物种树
iqtree2 -s concatenated.out -p partitions.txt -m MFP+MERGE -mset LG,JTT,WAG -mrate G,C20,C60 --ufboot 1000 --alrt 1000 --bnni -T 12 -pre concatenate
