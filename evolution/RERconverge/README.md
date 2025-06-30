# RERconverge
## Step1 
下载蛋白质序列文件和注释结果，比如pig：https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/ 
GCF_000003025.6_Sscrofa11.1_protein.faa.gz  
GCF_000003025.6_Sscrofa11.1_feature_table.txt.gz  

```
使用get_largest_protein.py从faa蛋白质序列文件中提取每个基因的最长蛋白序列作为代表。
Usage
python get_protein.py \
  --input-feature-txt GCF_000003025.6_Sscrofa11.1_feature_table.txt \
  --species pig \
  --input-protein-faa GCF_000003025.6_Sscrofa11.1_protein.faa
```
## Step2
```
blastp比对获得与人类同源基因的相互最佳比对
blastp.sh
```
## Step3
这一步是使用mafft进行多序列比对，并进行trimal修剪序列。  
注意：  
1.首先将所有物种的最长蛋白质序列的fa文件合并为一个总的fa文件。  
2.step2生成的rbh_combined.txt，需要将其物种名改为跟后续进化树中的物种名保持一致。如下：第3列就是进化树中的物种名，来源于ucsc 100-way alignment（https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/）
3.具体物种名如下：
|name1|name2|name3|
|---|---|---|
|Chlorocebus_sabaeus	|green_monkey	|chlSab2|
|Tupaia chinensis 	|chinese_tree_shrew	|tupChi1|
|Jaculus jaculus	|lesser_Egyptian_jerboa	|jacJac1|
|Sus scrofa	|pig	|susScr3|
|Bos taurus	|cattle	|bosTau8|
|Ovis aries|	sheep	|oviAri3|
|Canis lupus familiaris|	dog	|canFam3|
|Pteropus alecto 	|black_flying_fox	|pteAle1|
|Pteropus vampyrus	|large_flying_fox	|pteVam1|
|Eptesicus fuscus |	big_brown_bat	|eptFus1|
|Sorex araneus	|european_shrew	|sorAra2|
|Condylura cristata	|star_nosed_mole	|conCri1|
|Elephantulus edwardii	|cape_elephant_shrew	|eleEdw1|
|Chrysochloris asiatica	|cape_golden_mole	|chrAsi1|
|Orycteropus afer afer	|aardvark	|oryAfe1|
|Dasypus novemcinctus	|armadillo	|dasNov3|
|Monodelphis domestica	|opossum	|monDom5|
|Sarcophilus harrisii	|tasmanian_devil	|sarHar1|
|Ficedula albicollis 	|collared_flycatcher	|ficAlb2|
|Zonotrichia albicollis	|white_throated_sparrow	|zonAlb1|
|Geospiza fortis	|medium_ground_finch	|geoFor1|
|Pseudopodoces humilis	|tibetan_ground_jay	|pseHum1|
|Anolis carolinensis	|green_anole	|anoCar2|
|Xenopus tropicalis	|tropical_clawed_frog	|xenTro7|
|Latimeria chalumnae	|coelacanth	|latCha1|
|Otolemur garnettii	|small_eared_galago	|otoGar3|
|Pantholops hodgsonii	|tibetan_antelope	|panHod1|
|Gallus gallus	|chicken	|galGal4|
```
mafft比对
Usage:
python mafft.py
```
