参考教程：https://htmlpreview.github.io/?https://github.com/inzilico/imputeqc/blob/master/vignettes/k_selection.html
plink --bfile cattle_BW --chr chr1 --recode fastphase-1chr --out cattle_BW_chr1 --chr-set 31
# bed bim fam文件
Rscript make_test_files.R cattle_BW_chr1.recode.phase.inp
#这一步非常慢 因为snp数量多
