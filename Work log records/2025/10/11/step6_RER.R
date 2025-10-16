setwd("/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250926_new5species/5.RER")
library(RERconverge)
library(dplyr)
library(ape)
library(purrr)
estimatePhangornTreeAll(alndir = '/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/3.protein_singlecopy_orthologs',treefile = '/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/4.tree/10species.tree',output.file = '/home/bzheng/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/5.RER/geneTree/10species_geneTree')

tree <- readTrees('../5.RER/geneTree/10species_geneTree')
res <- getAllResiduals(tree,transform = 'sqrt',weighted = T,scale = T)
marineb <- read.tree('../4.tree/10species_noCGM.tree')
phenvMarine=tree2Paths(marineb, tree)
cor_res=correlateWithBinaryPhenotype(res, phenvMarine, min.sp=10, min.pos=5,weighted="auto")
write.table(cor_res,file = '../5.RER/RER_res_251017.txt',sep = '\t')

foreground = c("cattle","pig","dog","chicken","sheep")
root_sp<-"rock_ptarmigan"
masterTree <- tree$masterTree
perms = getPermsBinary(numperms=10000,
                       fg_vec=foreground,
                       sisters_list=NULL,
                       root_sp=root_sp,
                       RERmat=res,
                       trees=tree,
                       mastertree=masterTree,
                       permmode="cc",
                       min.pos=5,
                       calculateenrich=F,
                       annotlist=NULL)
corpermpvals=permpvalcor(cor_res,perms)
result=cor_res
result$permpval = corpermpvals[match(rownames(result), rownames(corpermpvals)), "permpval"]
result$permpvaladj=p.adjust(result$permpval, method="BH")
write.table(result,file = '../5.RER/RER_perms_res_251017.txt',sep = '\t')
