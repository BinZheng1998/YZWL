awk -F "\t" '$3 == "gene" {
    if ($9 ~ /gene ".*HOX/) {
        
        # 提取 gene " 之后的内容
        split($9, a, "gene \"");
        split(a[2], b, "\"");
        name = b[1];
        
        if (name ~ /HOX/) {
            print $1, $4, $5, name
        }
    }
}' OFS="\t" GCF_000001405.40_GRCh38.p14_genomic.gtf > hox_genes.txt

bedtools getfasta -fi GCF_000001405.40_GRCh38.p14_genomic.fna -bed hox_genes.txt -name -fo hox_genes.fa

~/software/ncbi-blast-2.15.0+/bin/blastn -query hox_genes.fa -db /home/Ref/ref/v23/chicken -outfmt 6 -evalue 1e-4 -out hox_genes_blastn.txt -num_threads 4

awk 'BEGIN{OFS="\t"} {
    if ($9 < $10) {
        print $2, $9, $10, $1
    } else {
        print $2, $10, $9, $1
    }
}' hox_genes_blastn.txt > hox_genes_blastn.bed

bedtools intersect \
  -a hox_genes_corrected.bed \
  -b <(awk '$3 == "gene"' /home/Ref/ref/huxu/20250529_chicken.gtf) \
  -wa -wb > final_gene_overlap.txt
