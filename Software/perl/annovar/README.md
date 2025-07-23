# https://github.com/WGLab/doc-ANNOVAR

# Step1
gffread pig.gff -T -o pig.gtf
~/software/gtfToGenePred -genePredExt pig.gtf pig_refGene.txt
perl retrieve_seq_from_fasta.pl --format refGene --seqfile pig.fa pig_refGene.txt --outfile pig_refGeneMrna.fa  
# Step2
perl convert2annovar.pl -format vcf4old pig_DP4_2allele_0.05miss_0.01maf_snps_rmdup_rmsingleSNPchr.vcf.gz -out pig
# Step3
