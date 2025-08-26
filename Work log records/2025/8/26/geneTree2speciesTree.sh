cat *treefile > ../tree/geneTree2speciesTree/geneTree.nwk
#参照ASTER软件的建议，先用wastral构建树，然后用astral4进行分支计算
wastral -S -i geneTree.nwk -t 48 -o wrastral.nwk
astral4 -i geneTree.nwk -g wrastral.nwk -t 48 -o astral4.nwk
