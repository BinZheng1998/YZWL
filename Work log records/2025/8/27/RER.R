setwd("/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250827_RER/res_geneTree2speciesTree")
library(RERconverge)
library(dplyr)
library(ape)
library(purrr)
#estimatePhangornTreeAll(alndir = '../mafft_trimal_0905/', treefile = '../tree/31species_noMouse.tree',output.file = '31species_geneTree',type="AA", format="fasta")
tree <- readTrees('31species_geneTree')
res <- getAllResiduals(tree,transform = 'sqrt',weighted = T,scale = T)
marineb <- read.tree('../tree/31species_noMouse_noCGM.tree')
phenvMarine=tree2Paths(marineb, tree)
cor_res=correlateWithBinaryPhenotype(res, phenvMarine, min.sp=20, min.pos=10,weighted="auto")
write.table(cor_res,file = 'RERconverge_cor.txt',sep = '\t')

foreground = c("cattle","pig","dog","chicken","duck","goose","sheep","goat","turkey","rabbit","cat","donkey","horse","water_buffalo","japanese_quail")
root_sp<-"japanese_quail"
masterTree <- tree$masterTree
perms = getPermsBinary(numperms=10000,
                       fg_vec=foreground,
                       sisters_list=NULL,
                       root_sp=root_sp,
                       RERmat=res,
                       trees=tree,
                       mastertree=masterTree,
                       permmode="cc",
                       min.pos=10,
                       calculateenrich=F,
                       annotlist=NULL)
corpermpvals=permpvalcor(cor_res,perms)
result=cor_res
result$permpval = corpermpvals[match(rownames(result), rownames(corpermpvals)), "permpval"]
result$permpvaladj=p.adjust(result$permpval, method="BH")
write.table(cor_res,file = 'RERconverge_cor_perms.txt',sep = '\t')
