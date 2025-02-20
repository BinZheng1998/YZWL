# 2025-01-06
## r code for manhattan plot

### fst
#Using vcftools to generate fst with window and step  
input fst file  
|CHROM	|BIN_START	|BIN_END	|N_VARIANTS	|WEIGHTED_FST	|MEAN_FST|
| --- | --- | --- | --- | --- | ---| 
|1	|1	|50000	|293	|0.0162935	|0.0135958|
|1	|10001	|60000	|277	|0.0158246	|0.0131814|
|1	|20001	|70000	|178	|0.019789	|0.0163616|
|1	|30001	|80000	|92	|0.023359	|0.016428|

### pi
#Using vcftools to generate pi with window and step
input pi files<br>

|CHROM	|BIN_START	|BIN_END	|N_VARIANTS	|PI|  
| --- | --- | --- | --- | --- | ---|  
|1	|1	|10000	|31	|0.000272786|  
|1	|5001	|15000	|99	|0.000776788|  
|1	|10001	|20000	|100	|0.000820545|  
|1	|15001	|25000	|103	|0.00112051|  


### xpehh xpnsl
selscan2.3.0<br>
https://github.com/szpiech/selscan<br>
norm

input xpehh file<br>
chr10	1	5001	0	-1	-1	-1	-1	NA	NA<br>
chr10	5001	10001	48	0	0	100	100	-0.695601	-1.9367<br>
chr10	10001	15001	10	0	0	5	5	-1.8235	-1.93116<br>
chr10	15001	20001	0	-1	-1	-1	-1	NA	NA<br>
chr10	20001	25001	3	-1	-1	-1	-1	NA	NA

input xpnsl file<br>
chr10	1	5001	0	-1	-1	-1	-1	NA	NA<br>
chr10	5001	10001	48	0	0	100	100	-0.753844	-1.55403<br>
chr10	10001	15001	10	0	0	1	5	-1.53138	-1.66723<br>
chr10	15001	20001	0	-1	-1	-1	-1	NA	NA<br>
chr10	20001	25001	3	-1	-1	-1	-1	NA	NA<br>
chr10	25001	30001	12	0	0	1	5	-1.331	-1.53136

#plot

![image](https://github.com/binzhengbin/YZWL/blob/main/evolution/R/manhattan/png/FST.png)
