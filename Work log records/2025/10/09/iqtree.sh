iqtree2 -s concatenated.out -m MFP --ufboot 1000 --bnni -T 128 -pre concatenate
iqtree2 -s concatenated.out -m MFP+MERGE --ufboot 1000 --bnni -T 128 -pre concatenate_2
