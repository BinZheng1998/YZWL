# get gene from gff/gtf
```
python extract.py input.txt ref.gff/gtf 0kb ./output.gene 
```
#input.bed  
|chr|start|end|
|---|---|---|
|chr1|300001|315000|
|chr1|355001|410000|
|chr1|705001|740000|
|chr1|1990001|2000000|
|chr1|2080001|2090000|
|chr1|2245001|2280000|
|chr1|2375001|2445000|
|chr1|3270001|3285000|
|chr1|3615001|3630000|
|chr1|3735001|3745000|

# get genes from 5species gtf/gff
Usage
```
用法: get_genes_202506.sh --input <输入文件> [选项]
选项:
  --input    输入文件路径（必需）
  --species  物种（cattle/dog/pig/sheep/chicken）[默认: cattle]
  --output   输出目录路径 [默认: 当前目录]
```
