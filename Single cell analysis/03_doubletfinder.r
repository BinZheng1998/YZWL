#!/usr/bin/env Rscript
options(future.globals.maxSize = 60 * 1024^3)

args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat("
Usage: Rscript test.R input.rds --pc-num <PC_NUM> --resolution <RES> --DoubletRate <RATE> --pN <PN> --sct <SCT> --out <OUTPUT>

Arguments:
  input.rds       The input Seurat .rds file
  --pc-num        Number of PCs to use for dimensionality reduction (e.g., 15)
  --resolution    Resolution for clustering (e.g., 1.0)
  --DoubletRate   Expected doublet rate (e.g., 0.075 for 7.5%)
  --pN            pN parameter for DoubletFinder (e.g., 0.25)
  --sct           Use SCTransform for normalization (TRUE or FALSE)
  --out           Output RDS file path

Example:
  Rscript test.R input.rds --pc-num 15 --resolution 1.0 --DoubletRate 0.075 --pN 0.25 --sct TRUE --out output.rds
  \n")
  quit(save = "no")
}

library(Seurat)
library(DoubletFinder)

if (length(args) < 12) {
  stop("Usage: Rscript test.R input.rds --pc-num <PC_NUM> --resolution <RES> --DoubletRate <RATE> --pN <PN> --sct <SCT> --out <OUTPUT>")
}

input_file <- args[1]

pc_num_idx <- which(args == "--pc-num") + 1
pc_num <- as.numeric(args[pc_num_idx])

resolution_idx <- which(args == "--resolution") + 1
resolution <- as.numeric(args[resolution_idx])

doublet_rate_idx <- which(args == "--DoubletRate") + 1
doublet_rate <- as.numeric(args[doublet_rate_idx])

pn_idx <- which(args == "--pN") + 1
pn <- as.numeric(args[pn_idx])

sct_idx <- which(args == "--sct") + 1
sct <- as.logical(args[sct_idx])

output_idx <- which(args == "--out") + 1
output_file <- args[output_idx]

cat("Using the following parameters:\n")
cat("PC num:", pc_num, "\n")
cat("Resolution:", resolution, "\n")
cat("DoubletRate:", doublet_rate, "\n")
cat("pN:", pn, "\n")
cat("SCT:", sct, "\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")

scRNA <- readRDS(input_file)
scRNA <- RunPCA(scRNA, verbose = FALSE)
scRNA <- RunUMAP(scRNA, dims = 1:pc_num, verbose = FALSE)
scRNA <- FindNeighbors(scRNA, dims = 1:pc_num)
scRNA <- FindClusters(scRNA, resolution = resolution)
sweep.res <- paramSweep(scRNA, PCs = 1:pc_num, sct = sct)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
pK <- as.numeric(as.character(pK))

homotypic.prop <- modelHomotypic(scRNA$seurat_clusters)

nExp_poi <- round(doublet_rate * nrow(scRNA@meta.data))

nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

scRNA <- doubletFinder(scRNA, PCs = 1:pc_num, pN = pn, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = sct)

classification_col <- grep("^DF.classifications_", colnames(scRNA@meta.data), value = TRUE)
pann_col <- grep("^pANN_", colnames(scRNA@meta.data), value = TRUE)

if (length(classification_col) > 0) {
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == classification_col] <- "Is_Double"
}
if (length(pann_col) > 0) {
  colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == pann_col] <- "Double_score"
}

doublet_rate_calc <- sum(scRNA@meta.data$Is_Double == "Doublet") / nrow(scRNA@meta.data)
cat("Doublet detection complete. Doublet rate:", round(doublet_rate_calc * 100, 2), "%\n")

scRNA_singlets <- subset(scRNA, Is_Double == "Singlet")

saveRDS(scRNA_singlets, file = output_file)
cat("Seurat object with singlets saved to", output_file, "\n")
