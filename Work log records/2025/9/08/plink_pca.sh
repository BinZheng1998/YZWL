plink --vcf ~/Pop/dna/snp/05_clean_vcf_phased/2036samples/merged_beagle.vcf.gz \
  --recode \
  --out chicken_2036samples \
  --const-fid \
  --allow-extra-chr \
  --chr-set 41

plink \
  --allow-extra-chr \
  --chr-set 41 \
  --file chicken_2036samples \
  --noweb \
  --make-bed \
  --out chicken_2036samples

plink \
  --allow-extra-chr \
  --chr-set 41 \
  --threads 48 \
  -bfile chicken_2036samples \
  --pca 30 \
  --out chicken_2036samples
