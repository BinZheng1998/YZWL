#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if ("-h" %in% args || "--help" %in% args) {
  cat("
Usage: Rscript 02_normalize.r --rds <input.rds> --method <sct/NormalizeData> --vars-to-regress <vars> --out <output.rds>

Options:
  --rds              Input RDS file containing the Seurat object
  --method           Normalization method to use: 'sct' for SCTransform or 'NormalizeData' for LogNormalize
  --vars-to-regress  (Optional) Variables to regress out, e.g. percent.mt nCount_RNA nFeature_RNA
  --out              Output RDS file to save the processed Seurat object
  -h, --help         Show this help message and exit

Example:
  Rscript 02_normalize.r --rds input.rds --method sct --vars-to-regress percent.mt nCount_RNA nFeature_RNA --out output.rds
")
  quit(save = "no")
}

if (length(args) < 6) {
  stop("Error: You need to provide at least 6 arguments: --rds <input.rds> --method <sct/NormalizeData> --out <output.rds>, with optional vars-to-regress.")
}

rds_idx <- which(args == "--rds") + 1
method_idx <- which(args == "--method") + 1
out_idx <- which(args == "--out") + 1
vars_to_regress_idx <- which(args == "--vars-to-regress") + 1

if (rds_idx > length(args) | method_idx > length(args) | out_idx > length(args)) {
  stop("Error: Missing required parameters: --rds, --method, or --out.")
}

input_rds <- args[rds_idx]
method <- args[method_idx]
output_rds <- args[out_idx]

vars_to_regress <- NULL
if (vars_to_regress_idx < out_idx) {
  vars_to_regress <- args[vars_to_regress_idx:(out_idx - 2)]
  cat("Variables to regress: ", paste(vars_to_regress, collapse=", "), "\n")
}

library(Seurat)
options(future.globals.maxSize = 15 * 1024^2 * 1024)
seurat_obj <- readRDS(input_rds)

if (method == "NormalizeData") {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst",nfeatures = 3000)
  if (!is.null(vars_to_regress)) {
    seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), vars.to.regress = vars_to_regress)
  } else {
    seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
  }
  cat("LogNormalize with ScaleData\n")
} else if (method == "sct") {
  if (!is.null(vars_to_regress)) {
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars_to_regress, verbose = FALSE)
  } else {
    seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  }
  cat("SCTransform with specified variables\n")
} else {
  stop("Unknown method: Use 'NormalizeData' or 'sct'.")
}

# 保存结果
saveRDS(seurat_obj, file = output_rds)
cat("Seurat object saved to:", output_rds, "\n")
