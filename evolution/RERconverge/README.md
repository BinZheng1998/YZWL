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
## Step4
从ucsc100 way alignment中提取所需的进化树拓扑结果，或者从Timetree网站上直接构建进化树。
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/  
http://www.timetree.org/  
#代码如下：  
```
library(ape)
#从ucsc-100way主树中提取28个物种的子树
ucsctree <- read.tree('hg38.100way.nh')
species <- read.table('species_name.txt',sep = '\t',header = T)
subtree <- keep.tip(ucsctree,species$name3)
plot(subtree)
write.tree(subtree, file = "ucsc100way_species_tree.nwk")
```
## Step5
1.将mafft多序列比对、trimal进行修剪后的结果放置一个新的文件夹，比如RERconverge_input  
2.接下来使用R包RERconverge进行分析，具体见RERconverge.R内容。
#### 注意：如果某个基因的多序列比对结果存在不同物种长度不一致，会导致RERconverge报错，可以将其里面的特殊氨基酸（比如J）转化为X。使用MSA-length-different.r代码进行检测。或者直接跳过那些基因。  
标准氨基酸符号如下：  
https://zh.wikipedia.org/wiki/%E6%A8%99%E6%BA%96%E8%9B%8B%E7%99%BD%E8%83%BA%E5%9F%BA%E9%85%B8%E5%88%97%E8%A1%A8
