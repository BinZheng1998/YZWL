# 请看完整个README
##  1
```bash
awk '{a[$2] = a[$2] ? a[$2] "," $1 : $1} END {for (i in a) print i, a[i]}' cattle_metadata.txt > cattle_metadata2.txt ```

#这个代码可以将下方的格式改为beagle.sh需要的input_sample.txt格式

#原始格式：

ACKR000001      Cocker-Spaniel
ACKR000002      Cocker-Spaniel
ACKR000003      Cocker-Spaniel
AESK000001      American-Eskimo-Dog
AESK000002      American-Eskimo-Dog
AESK000003      American-Eskimo-Dog
AESK000004      American-Eskimo-Dog
AESK000005      American-Eskimo-Dog
AESK000006      American-Eskimo-Dog
AFFN000001      Affenpinscher
AFFN000002      Affenpinscher
AFFN000003      Affenpinscher
AFFN000004      Affenpinscher
AFFN000005      Affenpinscher
AFFN000006      Affenpinscher
AFGH000001      Afghan-Hound
AFOX000001      American-Foxhound
AFOX000002      American-Foxhound
AFOX000003      American-Foxhound
AFOX000004      American-Foxhound
AFOX000005      American-Foxhound
AFOX000006      American-Foxhound
AFOX000007      American-Foxhound

#修改后的：

Wire-Fox-Terrier FXTE000001,FXTE000002,WIFX000001,WIFX000002
Polish-Lowland-Sheepdog POLS000001,POLS000002
Boxer BOXR000001,BOXR000002,BOXR000003,BOXR000004,BOXR000005
American-Water-Spaniel AMWS000001,AMWS000002,AMWS000003,AMWS000004,AMWS000005,AMWS000006
Estrela-Mountain-Dog ESMD000001,ESMD000002,ESMD000003,ESMD000004,ESMD000005,ESMD000006
American-Hairless-Terrier AHRT000001,AHRT000002,AHRT000003,AHRT000004,AHRT000005
Slovak-Cuvac SLCU000001,SLCU000002,SLCU000003,SLCU000004,SLCU000005,SLCU000006,SLCU000007
Giant-Schnauzer GSNZ000001,GSNZ000002,GSNZ000003,GSNZ000004,GSNZ000005,GSNZ000006
Thailand VILLTH000001,VILLTH000002,VILLTH000003,VILLTH000004,VILLTH000005,VILLTH000006,VILLTH000007,VILLTH000008,VILLTH000009,VILLTH000010
English-Foxhound EFXH000002,EFXH000003

#  2
#如果vcf文件中出现单倍体使用下述代码进行修改，修改前请对vcf文件进行简化

```zcat chr36_DP4_SNP_2allele_50miss_maf001_HWE1e6_simply.vcf.gz | awk 'BEGIN {FS=OFS="\t"} {for (i=10; i<=NF; i++) if ($i == ".") $i = "./."; print}' | bgzip -c > test.useimputation.vcf.gz```

#简化vcf
```bcftools annotate --remove QUAL,FILTER,INFO,^FORMAT/GT snp.pass.vcf.gz -Oz -o snp.simply.vcf.gz```
