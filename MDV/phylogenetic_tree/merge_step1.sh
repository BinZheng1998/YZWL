#!/bin/bash
#将单拷贝同源基因序列合并
#https://www.jianshu.com/p/92c9e0c70d5b
ls *fa | while read file; do 
    mafft --auto "$file" > "${file}.aln"
done &&

ls *aln | while read file; do 
    seqkit sort "$file" | seqkit seq -w 0 > "${file}.format"
done &&

paste -d " " *format > all.aln.fa
