setwd("~/project/01_evolution/convergent_evolution/07_result/20250923_RER_5species/0.script/")
library(RERconverge)
library(dplyr)
library(ape)
library(purrr)

estimatePhangornTreeAll(alndir = '/home/bzheng/project/01_evolution/convergent_evolution/07_result/20250923_RER_5species/2.input.fa/trimmed', 
                        treefile = '../1.tree/14_species_concatenate_tree.nwk',
                        output.file = '../3.res/geneTree/14species_geneTree')
tree <- readTrees('../3.res/geneTree/14species_geneTree')
res <- getAllResiduals(tree,transform = 'sqrt',weighted = T,scale = T)
marineb <- read.tree('../1.tree/14_species_concatenate_tree_noCGM.nwk')
phenvMarine=tree2Paths(marineb, tree)
cor_res=correlateWithBinaryPhenotype(res, phenvMarine, min.sp=10, min.pos=5,weighted="auto")
write.table(cor_res,file = '../3.res/RER/RER_transcript_250924.txt',sep = '\t')
