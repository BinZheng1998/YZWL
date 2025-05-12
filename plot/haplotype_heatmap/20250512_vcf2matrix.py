#!/usr/bin/env python
import sys
import pandas as pd

# 检查命令行参数
if len(sys.argv) != 3:
    print("用法: python vcf2matrix.py <input.vcf> <out.txt>")
    sys.exit(1)

vcf_file = sys.argv[1]  # 输入 VCF 文件
txt_file = sys.argv[2]  # 输出 TXT 文件

# 读取 VCF 文件
try:
    with open(vcf_file, 'r') as f:
        # 找到头部，提取样本名
        for line in f:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.strip().split('\t')
                samples = header[9:]  # 样本名从第 10 列开始
                break
        else:
            raise ValueError("VCF 文件中未找到头部")

        # 初始化数据结构
        data = {sample: [] for sample in samples}
        pos_list = []

        # 处理每个变异行
        for line in f:
            if line.strip() == '':
                continue
            fields = line.strip().split('\t')
            
            # 提取相关字段
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            
            # 过滤为双等位基因 SNP
            if len(ref) != 1 or len(alt.split(',')) != 1 or len(alt) != 1:
                continue
            
            # 获取 FORMAT 和 GT 索引
            format_field = fields[8]
            if 'GT' not in format_field.split(':'):
                continue
            gt_index = format_field.split(':').index('GT')
            
            # 为每个样本提取并转换基因型
            for i, sample in enumerate(samples):
                sample_field = fields[9 + i]
                sample_values = sample_field.split(':')
                gt = sample_values[gt_index]
                
                # 处理分相（|）或未分相（/）基因型
                if '/' in gt:
                    alle1, alle2 = map(int, gt.split('/'))
                elif '|' in gt:
                    alle1, alle2 = map(int, gt.split('|'))
                else:
                    continue
                
                # 转换基因型：0|0 -> 0, 0|1/1|0 -> 1, 1|1 -> 2
                genotype = alle1 + alle2
                data[sample].append(genotype)
            
            pos_list.append(pos)

    # 检查是否有有效数据
    if not pos_list:
        raise ValueError("未找到符合条件的 SNP 位点")

    # 创建 DataFrame
    df = pd.DataFrame(data).T
    df.index = samples  # 设置样本为行
    df.columns = pos_list  # 设置 SNP 位置为列

    # 写入 TXT 文件
    df.to_csv(txt_file, sep='\t', index_label='Sample')

except FileNotFoundError:
    print(f"错误: 找不到输入文件 {vcf_file}")
    sys.exit(1)
except Exception as e:
    print(f"错误: {str(e)}")
    sys.exit(1)
