setwd("~/project/10_RER/result/protein_MSA_result/")
library(RERconverge)
library(ape)
#从ucsc-100way主树中提取28个物种的子树
#ucsctree <- read.tree('hg38.100way.nh')
#species <- read.table('species_name.txt',sep = '\t',header = T)
#subtree <- keep.tip(ucsctree,species$name3)
#plot(subtree)
#write.tree(subtree, file = "ucsc100way_species_tree.nwk")

#统计序列长度
#raw_seqs <- read.FASTA("test/NP_005036.2_trimmed.fa", type = "AA")
#lengths <- nchar(raw_seqs)
#print(table(lengths))
estimatePhangornTreeAll(alndir = 'RERconverge_input/',treefile = 'ucsc100way_species_tree.nwk',output.file = 'RERconverge_output/geneTree')
tree <- readTrees('RERconverge_output/geneTree')
res <- getAllResiduals(tree,transform = 'sqrt',weighted = T,scale = T)
#这里的树其它物种分支长度为0，目标物种为1，需要对tree进行修改。但是拓扑结果不能变。ucsc100way_species_tree_noCGM.nwk
marineb <- read.tree('ucsc100way_species_tree_noCGM.nwk')

phenvMarine=tree2Paths(marineb, tree)
cor_res=correlateWithBinaryPhenotype(res, phenvMarine, min.sp=10, min.pos=4,weighted="auto")
cor_res
write.table(cor_res,file = 'RERconverge_res.txt',sep = '\t')

#置换检验
foreground = c("bosTau8","susScr3","canFam3","galGal4","oviAri3")
root_sp<-"latCha1"
masterTree <- tree$masterTree
perms = getPermsBinary(numperms=10, 
                       fg_vec=foreground, 
                       sisters_list=NULL, 
                       root_sp=root_sp, 
                       RERmat=res, 
                       trees=tree,
                       mastertree=masterTree,
                       permmode="cc",
                       min.pos=4,
                       calculateenrich=F,
                       annotlist=NULL)
corpermpvals=permpvalcor(cor_res,perms)
result=cor_res
result$permpval = corpermpvals[match(rownames(result), rownames(corpermpvals)), "permpval"]
result$permpvaladj=p.adjust(result$permpval, method="BH")
 
#可视化
par(mfrow=c(1,1))
phenvExample<-foreground2Paths(c("susScr3","bosTau8","oviAri3","canFam3","galGal4"),tree,clade="terminal")
plotRers(res,"XP_047281837.1_trimmed",phenv=phenvExample) #plotRERs

#添加人类基因名称
humanid <- read.table('~/project/10_RER/data/human.txt',sep = '\t')
colnames(humanid) <- c('Gene',"NCBIid",'Length')
result$Gene <- rownames(result)
result$Gene <- gsub('_trimmed','',result$Gene)
RER_res <- merge(result,humanid,by.x='Gene',by.y='NCBIid')
colnames(RER_res)[8] <- 'human'

#添加动物基因名称
animalid <- read.table('~/project/10_RER/result/RBBH_result/rbh_rename.txt',sep = '\t',header = T)
animalid <- animalid[,c(1:6)]
colnames(animalid)<- c('human','pig','dog','sheep','cattle','chicken')
geneid <- read.table('~/project/10_RER/data/human_chicken_pig_cattle_sheep_dog.txt',sep = '\t')
map_vec <- setNames(geneid[, 1], geneid[, 2])
cols_to_replace <- 2:6
animalid[, cols_to_replace] <- lapply(animalid[, cols_to_replace], function(col) {
  ifelse(col %in% names(map_vec), map_vec[col], col)
})

#最终表格
RER_res <- merge(RER_res,animalid,by.x='Gene',by.y='human')
RER_res <- RER_res[,-9]
write.table(RER_res,file = 'RERconverge_res_20250708.txt',sep = '\t')

#GO
rer_gene <- RER_res[RER_res$P<0.05,]
rer_gene<-na.omit(rer_gene)
library(clusterProfiler)
library(org.Gg.eg.db)
library(DOSE)
library(GO.db)
library(aPEAR)
chicken_genes<- rer_gene$chicken
eg <- bitr(chicken_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
id <- as.character(eg[, 2])
ego <- enrichGO(gene = id, OrgDb = "org.Gg.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
ego <- as.data.frame(ego)
ego <- ego[ego$pvalue < 0.05,]

p<-enrichmentNetwork(ego,colorBy = 'pvalue',
                     colorType = 'pval',
                     nodeSize = 'zScore',
                     fontSize = 2,
                     drawEllipses = T,
                     verbose = F)
p
