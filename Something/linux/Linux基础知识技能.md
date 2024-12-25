 

# Linux基础知识技能

###### --来源：B站：生信技能树

## 一、一些重要的基础命令

| 命令                    | 功能                                                         | 举例                                                         | 含义           |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | -------------- |
| echo                    | 内容打印到屏幕                                               |                                                              |                |
| w                       | 查看用户数量                                                 |                                                              |                |
| top                     | 查看当前运行命令                                             |                                                              |                |
| sleep                   | 系统暂停一顿时间                                             | sleep 10                                                     | 暂停10s        |
| &                       | 放入后台运行                                                 | sleep 10 &                                                   | 后台暂停10s    |
| which                   | 应用在那个地点                                               | which ls                                                     | ls在那个文件夹 |
| \|                      | 管道符，连接命令                                             |                                                              |                |
| pwd                     | 当前完整目录                                                 |                                                              |                |
| ls                      | 列出当前目录与文件                                           | ls boy                                                       |                |
| cd                      | 进入当前目录                                                 |                                                              |                |
| \> $log_file 2>&1       | 可以在程序末尾加上，运行的流程就会加到log文件中              |                                                              |                |
| >                       | 重定向                                                       |                                                              |                |
| 2>&1                    | 错误信息也会标准重定向输出                                   |                                                              |                |
| du -sh                  | 显示当前总目录的大小                                         |                                                              |                |
| ll -h                   | 列出当前目录下所有文件和文件夹的大小                         |                                                              |                |
| mkdir -p  ./data/raw    | 创建嵌套的文件目录，上级目录不存在也会 一起创建              |                                                              |                |
| cp -r  /O_file  /N_file | 递归复制目录及内容                                           | cp -r source_directory destination_directory                 |                |
| rm  -i                  | 删除前进行提示删除                                           |                                                              |                |
| rm -r                   | 递归删除目录                                                 |                                                              |                |
| rm -f                   | 强制删除                                                     |                                                              |                |
| ln     -f               | 创建软连接，-f:目标文件夹已经存在，强制删除 后再创建软连接   |                                                              |                |
| tree                    | 树状结构显示；-d：只显示目录；-f:显示文件和目录的完整结构    |                                                              |                |
| paste                   | 多个文件进行行连接； -d:使用分隔符;-s:所有内容合并成一行；   | paste -d  "," sample1.exp.bed sample2.exp.bed >combine.list2 |                |
| find                    | 路径表示要搜索的目录路径。如果未指定路径，则默认为当前目录   | find ./ -name "leaf1.html" -type f # 搜索文件 find ./ -name "leaf1" -type d # 搜索目录 |                |
| cat                     | 查看合并文件                                                 | cat file1.txt file2.txt >file3.txt                           |                |
| less -i file.txt        | 搜索关键词并高亮显示匹配的结果                               |                                                              |                |
| cut                     | -c：按字符位置列提字段；-f:按照字段位置；-d：指定分隔符      | 提取1-3列：cut -f 1,3 test.gff3  ; 第二个字符到第五个字符：cut -c 2-5 test.fa；逗号作为字段分隔符提取文件的第二列：cut -d ',' -f 2 exp.csv；提除第一列以外的所有列：cut --complement -f1 test.gff3 |                |
| sort                    | -g:字典顺序排序；-n:按数值顺序进行排序；-u：去重复行；-k:指定键排序； | sort -k4,4n test.gff3：按第4列，n表示按数字进行排序。        |                |
| uniq                    | -c:显示重复行出现次数；-d:仅显示重复行；-u：仅显示不重复的行； | cut -f 1 peak.xls\|sort\|uniq：去重复 ；cut -f 1 peak.xls\|sort\|uniq -c：显示文件中的重复行和出现次数；cut -f 1 peak.xls\|sort\|uniq -u：新鲜事不重复的行 |                |
| wc                      | -c:输出字节数；-w:词数；-l:输出行数                          | wc file.txt：统计文件的行数、词数和字符数；wc file1.txt file2.txt：输出多个文件的总行数、词数和字符数 |                |
| grep                    | -i:忽略大小写；-v:仅显示不匹配行；-r:递归目录及子目录搜索；-n:显示匹配结果所在行号。 -w:匹配完整的单词；-c:匹配解雇哦统计行数； | grep -i -n "gene4" exp.csv：忽略大小搜索模式并显示行号； grep -v "chrM" peak.xls：显示不匹配模式的行 |                |
| sed                     | -i:直接修改文件内容；-r:启用正则表达式语法；                 | sed 's/pattern/replacement/' file.txt：替换文件中的文本并输出；sed -i 's/_gene//g' exp.csv：使用正则表达式进行文件替换；sed "s/\r//" test.txt：替换文中的换行符；sed '2,3d' file.txt：删除文件中第2-3行 |                |
| awk                     | '$0':表示整行内容；'$1'，'$2':表示第1，2字段内容；-F:指定字段分隔符； | awk '/pattern/' filename：打印匹配整行；awk '{if($2 > 50) print $0 }' sample1.exp.bed：通过’if‘条件来控制行的输出；awk '{ total += $2 } END { print total }' stat.txt：将第二列结果进行累加；awk '{if($3=="gene"&&$5-$4>5000)print }' test.gff3：对gene长度大于5000的基因行进行输出；awk '{if(FNR>=2&&FNR <= 4)print }' file.txt：打印当前文件前三行； |                |
| gzip                    | -d：解压文件；-r:递归处理目录下所有文件；-k:保留原始文件；   | gzip -d R1_1.fq.gz：解压文件；gzip file1 file2 file3：压缩多个文件；gzip -c Peacomb_sample.recode.vcf > Peacomb_sample.recode.vcf.gz ：生成另一个压缩文件 |                |
| tar                     |                                                              | tar -zxvf  archive.tar.gz 解压缩.tar.gz 文件                 |                |
| unzip                   | 解压目录                                                     |                                                              |                |
| bzip2                   | -d 解压                                                      | 解压缩文件：bzip2 -d test.sam.bz2                            |                |
| nohup                   | 放入后台，关闭终端 后仍能使用                                |                                                              |                |
| jobs                    | 将后台任务显示：和fg %任务序号，可以将后台任务放入前台       |                                                              |                |
| ps                      | -e :显示所有进程；-u:指定用户名进程；-l:长格式显示进程信息   |                                                              |                |
| top                     | -d;指定刷新时间；-i:忽略僵尸进程；-u:指定用户形式；          |                                                              |                |
| kill                    |                                                              | kill 123 :请求终止ID为123的进程；kill -19 12345：进程挂起；kill -18 12345：请求进程回复 |                |
| df                      | -h:显示磁盘空间使用情况;                                     |                                                              |                |
| du                      | -a:每个目录和文件使用情况；-sh:指定目录使用情况；            |                                                              |                |
| passwd                  | 更改用户密码（自己）                                         |                                                              |                |
| chmod                   | +:添加权限；                                                 | chmod  u+x file.sh:给用户加上可执行权限                      |                |
| alias                   | 可以添加一些快捷方式操作                                     | alias gp='grep'；alias cp='cp -r'                            |                |
|                         |                                                              |                                                              |                |
|                         |                                                              |                                                              |                |

## 二、一些重要的流程命令组合

|      | 命令组合                                  | 功能                           |
| ---- | ----------------------------------------- | ------------------------------ |
| 1    | ps -ef \|grep                             | 用于linux中查找正在运行的进程  |
|      | ps -ef \|grep dds                         | 查找dds开头的进程              |
|      | kill PID                                  | 杀掉该进程                     |
| 2    | 电脑"TAB"键                               | 自动补全命令                   |
| 3    | sed -i '2,3d' file                        | 删除文件的第2行和第3行         |
| 4    | sed "s:^:`pwd`/:"                         | 显示所有文件的绝对路径         |
|      | sed 's/.//' temp.txt                      | 删除文件第一个字符             |
|      |                                           |                                |
|      | cat -T pop.list.txt                       | 判断文件是否使用制表符进行连接 |
|      | rename 's/.R1.fastq.gz/_S1' *.R1.fastq.gz | 批量进行文件后缀修改           |
|      | rename 's/^旧前缀_/新前缀_/' 旧前缀_*     |                                |

**#后台运行与远程传输**

```shell
#讲运行命令放入后台

nohup sh aspera_ATAC.sh &

#命令重新放回前台
jobs
#实例
(base) fguo@falcon:~/gf_workstation/02.ATAC-seq_analysis/Raw_data$ jobs
[1]-  Running                 nohup sh sra_explorer_fastq_aspera_download.sh &  (wd: ~/gf_workstation/03.Single_cell_analysis/Raw_data)
[2]+  Running                 nohup sh aspera_ATAC.sh &
#将2号命令放入前台
fg %2
#按Ctrl+Z即可以停止命令


##scp远程文件传输
scp -r directory_path remote_user@remote_host:remote_directory_path
```

**#硬盘挂载（把U盘文件拷贝到服务器中）：**

```shell


######硬盘挂载与传输####
# 1.插入存储盘，根据大小查看哪个插入盘
lsblk
# 2.创建挂载点,例如，创建‘/mmt/usb’
sudo mkdir /mmt/usb
# 3.挂载U盘
## 使用'mount'命令将u盘挂载到目录中。假如设备名为'/dev/sdb1'
sudo mount /dev/sdb1 /mnt/usb
# 4.复制文件 
sudo cp -r /mnt/usb/* /hoem/user
## 传输完之后  使用'umount'卸载硬盘
sudo umount /mnt/usb
# 5.清理挂载点 （不需要挂载点，删除创建目录）
sudo rmdir /mnt/usb
 

```

**#替换染色体名字和gff文件保持一致**

```shell

#!/bin/bash

########将gff文件名替换成fa文件名####
# 检查输入参数数量：fa文件、gff文件，New_fasta文件名
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <genome.fa> <annotations.gff> <output.fa>"
    exit 1
fi

# 文件名
fa="$1"
gff="$2"
fasta="$3"

# 存储fa文件的染色体名字
less -S $fa |grep ">"|sed 's/^>//' >temp.fa.name

# 存储gff染色体名字
  #按数字排序
grep '^chr[0-9]' $gff | sort -n -t 'r' -k 2 > sorted_numbers.txt
# 按ZW排序
grep -E '^chr[WZ]' $gff | sort > sorted_wz.txt
#最后用mito
grep '^mito' $gff > mito.txt
# 合并结果到一个文件
cat sorted_numbers.txt sorted_wz.txt mito.txt > temp.gff.name
# 删除中间文件
rm sorted_numbers.txt sorted_wz.txt mito.txt

##3.进行替换
# 创建新的文件，避免覆盖原始文件
cp $fa $fasta
、 
# 使用 paste 命令合并两个文件，使用竖线 | 作为分隔符
paste -d "|" temp.fa.name temp.gff.name | while IFS='|' read -r fa_name gff_name; do
    # 输出正在处理的信息
    echo "替换 $fa_name 为 $gff_name"
    
    # 使用 sed 命令替换内容，并将结果写入 huxu_new_genome.fa
    sed -i "s/$fa_name/$gff_name/g" $fasta
done
# 删除中间文件
rm temp.fa.name temp.gff.name

  echo "Fasta文件染色体名字替换完成!"
 
```

**批量修改的前缀名**

```shell
for file in ref*; do
  mv "$file" "huxu${file#ref}"
done
```

**并行执行命令**

```shell

#!/bin/bash

### 可以通过nohup和&实现任务的并行处理,使用’&‘将命令放入后台。可以启动多个进程，加快处理速度

# 会直接把进程放入后台进行并行运算
for i in $(cat /home/fguo/gf_workrstation/02.mapping/sample.txt)
do
   nohup fastqc /home/fguo/gf_workstation/01.Hi-C/data/"${i}_R1.clean.fq.gz" -o ./ &
   nohup fastqc /home/fguo/gf_workstation/01.Hi-C/data/"${i}_R2.clean.fq.gz" -o ./ &
done

# 等待所有后台任务的完成，避免尚未处理完成就提前推出脚本
wait

```

**#作封装成一个函数，并在脚本中调用，这种方式不仅使得代码更清晰易读，还能让你更容易地修改和扩展功能。（可以用来编写workflow流程）**

```shell
#!/bin/bash

# 定义一个函数来运行 fastqc
run_fastqc() {
    local sample="$1"
    nohup fastqc "/home/fguo/gf_workstation/01.Hi-C/data/${sample}_R1.clean.fq.gz" -o ./ &
    nohup fastqc "/home/fguo/gf_workstation/01.Hi-C/data/${sample}_R2.clean.fq.gz" -o ./ &
}

# 主程序
main() {
    # 读取 sample.txt 中的每一行作为参数调用函数，设置sample的输入路径和输出路径即可
    while IFS= read -r sample; do
        # 指示当前正在处理哪个样本。
        echo “run_fastqc $sample Start !” 
        run_fastqc "$sample"
    done < "/home/fguo/gf_workstation/02.mapping/sample.txt"

    # 等待所有后台任务完成
    wait
}

# 执行主程序
main


```

**写批量基本框架**

```shell

#!/bin/bash

# 定义一个函数来运行 fastqc
run_hic_buildMatrix() {

#定义局部变量样本名、参考基因组、模拟酶切位点
    local sample="$1"
    local ref="/home/fguo/gf_workstation/01.Hi-C/ref/huxu_genome.fa"
    local rest_site="/home/fguo/gf_workstation/01.Hi-C/ref/huxu_rest_site.bed" 
           
    hicBuildMatrix --samFiles "../01.mapping/02.HiC/aligned/${sample}_R1.bam" "../01.mapping/02.HiC/aligned/${sample}_R2.bam"  --genomeAssembly ${ref}  --outBam "./Raw/${sample}.bam" --QCfolder "./Raw/$sample" --restrictionCutFile ${rest_site} --binSize 1000  --threads 140 --minDistance 500  --maxLibraryInsertSize 1500 --inputBufferSize 400000 --restrictionSequence AGCT --danglingSequence AGCT  1>logs/${sample}_hicBuildMatrix.log 2>&1 &
}

# 主程序
main() {
    # 读取 sample.txt 中的每一行作为参数调用函数，设置sample的输入路径和输出路径即可
    while IFS= read -r sample; do
        # 指示当前正在处理哪个样本。
        echo “run_hic_buildMatrix $sample Start !” 
       run_hic_buildMatrix "$sample"
    done < "/home/fguo/gf_workstation/02.mapping/sample.txt"

    # 等待所有后台任务完成
    wait
}

# 执行主程序
main
```

**运行进程的控制**

```shell
## 终止进程
ctrl + C

## 暂停进程
ctrl + z

## 将后台程序放入前台
# 1.查看进程
jobs

# 结果
[6]-  Running                 samtools sort -n -@ 32 -o sorted_E55_G_Z_2_R1.bam E55_G_Z_2_R1.bam > sort_logs/E55_G_Z_2_sort_R1_bam.log 2>&1 &
[7]+  Running                 samtools sort -n -@ 32 -o sorted_E55_G_Z_2_R1.bam E55_G_Z_2_R1.bam > sort_logs/E55_G_Z_2_sort_R2_bam.log 2>&1 &

# 2.将进程放入前台，：使用 fg 命令，后面跟上作业的编号
fg %6

```

**终止后台进程**

```shell
# 方法1：使用’kill‘命令终止特定进程

# 1.列出所有后台
job -l

# 2.终止特定的进程
kill -9 <PID> #替换成需要终止的进程ID

# 方法2：终止所有后台进程
kill $(jiobs -p) # 获取后台所有进程PID，并终止它们

# 方法3：使用'pkill'命令
## 根据进程名称或者其他属性终止进程。
pkill hicBuildMatrix  # 终止hicBuildMatrix进程

# 方法4：终止所有用户进程 （谨慎使用！！！）
## 终止当前用户进程
pkill -u2（whoami） #用户名

# 清理后台进程脚本
#!/bin/bash

#### 自动清理后台脚本，脚本函数 ###
cleanup() {
      # 函数开始运行
      echo "Terminating all background processes..."
      kill $(jobs -p)
      # 进程运行完成后退出进程
      wait
}

# 在脚本中调用 cleanup 函数
main 
cleanup 



```

**juicer配置与分析**

```shell
**1.安装依赖软件**

java

perl 

python

GNU utils

bwa

samtools

**2.建立目录结构**

juicer软件要求一个固定目录，新建”juicer“目录，为软件安装目录，下面有4个子目录：

mkdir  reference \  ## 存放参考基因组

 work \                      ## 存放样本序列文件和分析结果

scripts \                      ## 存放需要的脚本文件

restriction_sites \      ##  存放参考基因组酶切图谱

3.下载juicer源代码和cuda源代码，放到**scripts**目录下。目录下juicer可以进行单机或者集群系统上运行，


```

**conda日常使用**

```shell
基础操作
查询 conda 版本
conda --version

更新 conda
conda update conda

查看conda环境详细信息
conda info

虚拟环境管理
查看当前有哪些虚拟环境
conda env list

或者使用如下命令：

conda info --envs

创建一个新的虚拟环境
conda create --name jupyter_venv python=3.8

其中，通过 -n或--name 来自定义的环境名称，如：jupyter_venv；同时，指定Python的版本。

激活虚拟环境
conda activate jupyter_venv

退出当前虚拟环境
conda deactivate

删除某个虚拟环境
conda remove -n your_env_name --all 其中，-n与--name等价，表示虚拟环境名

复制某个虚拟环境
conda create --name new_env_name --clone old_env_name

分享/备份一个虚拟环境
一个分享环境的快速方法就是给他一个你的环境的.yml文件。

首先激活要分享的环境，在当前工作目录下生成一个environment.yml文件。

conda env export > environment.yml

对方拿到environment.yml文件后，将该文件放在工作目录下，可以通过以下命令从该文件创建环境即可。

conda env create -f environment.yml

# 卸载软件
# 1.首先可以通过conda list查看环境中已经安装的软件包
conda list
# 2.进行卸载
conda uninstall your_package
# 3.使用pip卸载
pip uninstall your_package
# 安装指定版本的软件包
# 1.先进行搜索
conda search --full-name protobuf
# 2.再安装对应的软件
conda install protobuf=3.19.6
# 或着指定安装channel
conda install -c conda-forge protobu;f=3.19.6

```

#### SCp文件远程传输

```shell
# 使用scp命令从本地复制到远程服务器(文件) 
scp /path/to/local/file username@remote_host:/path/to/remote/directory
# 使用scp命令从本地复制到远程服务器(目录)
scp -r /path/to/local/directory username@remote_host:/path/to/remote/directory
# 使用scp命令从远程服务器复制到本地(文件)     
scp username@remote_host:/path/to/remote/file /path/to/local/directory 
# 使用scp命令从远程服务器复制到本地(目录)
scp -r username@remote_host:/path/to/remote/directory /path/to/local/directory

# 使用scp命令时跳过确认
scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no /path/to/local/file username@remote_host:/path/to/remote/directory 
# 传输文件时显示传输进度
scp -r -v /path/to/local/directory username@remote_host:/path/to/remote/directory
# 使用scp命令使用特定端口     
scp -P port /path/to/local/file username@remote_host:/path/to/remote/directory 
# 使用scp命令使用特定SSH密钥进行身份验证     
scp -i /path/to/key_file /path/to/local/file username@remote_host:/path/to/remote/directory 
1
```

### 原始数据下载

```
# 1.进入NCBI,选择左边BioSample
# 2.输入关键字。例如“gonad chicken”
# 3.进入点击BioProject
# 4.点击SRA Experiments右边
# 5.选择自己想要的数据点击进入“All runs”
# 6.选择自己需要的数据，点击AccessionList
# 7.获得表格后
# 进入SRA_Explorer输入BioProject号选择需要的数据,点击右上角购物车，然后选择Aspera下载，再对文件名字进行修改，复制到服务器下载即可
```



### 输入型流程程序使用

```shell
#!/bin/bash

# 检查输入参数的数量
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <VCF file> <Region>"
    exit 1
fi

# 获取输入参数
VCF_FILE=$1
REGION=$2

# 设置输出目录和文件名
OUTPUT_DIR="LD_block"
OUTPUT_FILE="${OUTPUT_DIR}/$(basename ${VCF_FILE} .vcf.gz)_$(echo ${REGION} | sed 's/:/-/g').png"

# 创建输出目录（如果不存在的话）
mkdir -p ${OUTPUT_DIR}

# 运行 LDBlockShow 工具
/home/bzheng/software/LDBlockShow-1.40/bin/LDBlockShow \
    -InVCF ${VCF_FILE} \
    -OutPut ${OUTPUT_FILE} \
    -Region ${REGION} \
    -OutPng \
    -SeleVar 1

# 提示用户操作完成
echo "LDBlockShow analysis complete. Output saved to ${OUTPUT_FILE}"
```

#### 文件查看与修改

```
## 1.对文件某一行进行查看
sed -n '337p' your_file  ## 对第337行内容进行查看

```

#### 在Linux环境下，安装R,和安装大量的R包

```
#### 1、R包的安装与下载 ####
# 1.常规下载
CRAN来源的一般可以使用常规下载（并且可以自动下载依赖包）：
install.packages('ggplot2')

# 2.下载源文件本地安装
如果不能顺利进行安装，可以下载R包的源码进行本地安装，下载的方法就是找到包的相关网站进行源文件的下载源文件。
```

![image-20241213102227619](C:\Users\11952\AppData\Roaming\Typora\typora-user-images\image-20241213102227619.png)

```
# 将源文件下载到本地之后，再使用如下的方法安装：
install.packages('ggplot2_3.5.1.tar.gz',type = 'source',repo=NULL)
```

### 查看文件数量

#### 1.包含子目录下的相关文件

```shell
ls -lR |grep "^_" |wc -l
```

#### 2.不包含子目录的目录

```shell
ls -l |grep "^d" |
```

#### 3.查看当前目录下的文件个数

#### 3.1 不包含目录中的目录

```shell
ls -l | grep "^d" |wc -l
```

#### 3.2 包含目录中的目录

```1111
ls -lR |grep "^d" |wc -l
```

#### 4.查看当前文件夹下较某某的文件的数量

```shell
find . -name filename |wc -l

# 统计当前目录下面所有的.png图片使用
find -name "*.png" |wc -l

## 查看某个文件夹中各个文件夹的大小
du -sh *
```



#### github的基本操作