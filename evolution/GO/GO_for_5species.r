#!/usr/bin/env Rscript

library(dplyr)
library(clusterProfiler)
library(DOSE)
library(GO.db)
library(org.Bt.eg.db)
library(org.Cf.eg.db)
library(org.Gg.eg.db)
library(org.Ss.eg.db)
library(AnnotationHub)
library(simplifyEnrichment)
library(rrvgo)
library(optparse)
library(ggplot2)

# 定义命令行参数解析
option_list <- list(
  make_option(c("-i", "--input-path"), type="character", default=NULL, 
              help="输入文件路径", metavar="character"),
  make_option(c("-s", "--species"), type="character", default=NULL,
              help="物种选择: chicken/dog/sheep/pig/cattle", metavar="character"),
  make_option(c("-o", "--output-folder"), type="character", default="./",
              help="输出文件夹路径", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 验证必要参数
if (is.null(opt$`input-path`) || is.null(opt$species)) {
  print_help(opt_parser)
  stop("必须提供 --input-path 和 --species 参数", call.=FALSE)
}

# 主分析函数
run_go_analysis <- function(input_path, species, output_folder) {
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  org_db <- switch(species,
    "chicken" = org.Gg.eg.db,
    "dog" = org.Cf.eg.db,
    "sheep" = {
      hub <- AnnotationHub::AnnotationHub()
      hub[['AH114632']]
    },
    "pig" = org.Ss.eg.db,
    "cattle" = org.Bt.eg.db,
    stop("不支持的物种，请选择: chicken/dog/sheep/pig/cattle")
  )

  gene_files <- list.files(path = input_path, pattern = "*genes.txt", full.names = TRUE)

  for (file in gene_files) {
    message("\n", rep("=", 50))
    message("开始执行文件: ", basename(file))
    message(rep("=", 50), "\n")

    gene <- read.table(file, header = TRUE, sep = "\t")
    genes_list <- strsplit(gene[[5]], ",")
    all_genes <- unlist(genes_list)
    unique_genes <- unique(all_genes)
    symbol <- gsub(" ", "", as.character(unique_genes))
    eg <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
    id <- as.character(eg[, 2])

    ego <- enrichGO(gene = id, OrgDb = org_db, ont = "all",
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                    readable = TRUE)

    ego <- as.data.frame(ego)
    ego <- ego[ego$pvalue < 0.05,]
    base_name <- tools::file_path_sans_ext(basename(file))
    write.table(as.data.frame(ego), file = file.path(output_folder, paste0(base_name, "_GO.txt")),
                sep = "\t", row.names = FALSE)

    go_id_bp <- ego[ego$ONTOLOGY == "BP", "ID"]
    go_id_cc <- ego[ego$ONTOLOGY == "CC", "ID"]
    go_id_mf <- ego[ego$ONTOLOGY == "MF", "ID"]

    if (length(go_id_bp) > 1) {
      tryCatch({
        mat1 <- GO_similarity(go_id_bp, ont = "BP", db = org_db, measure = "Sim_Resnik_1999")
        output_file2 <- file.path(output_folder, paste0(base_name, "_GO_BP_simplifyEnrichment.pdf"))
        pdf(output_file2, width = 12, height = 8)
        simplifyGO(mat1, plot = TRUE)
        dev.off()
        message("已生成BP simplifyEnrichment结果: ", output_file2)
      }, error = function(e) {
        message("BP GO简化失败: ", e$message)
      })
    }

    if (length(go_id_cc) > 1) {
      tryCatch({
        mat2 <- GO_similarity(go_id_cc, ont = "CC", db = org_db, measure = "Sim_Resnik_1999")
        output_file3 <- file.path(output_folder, paste0(base_name, "_GO_CC_simplifyEnrichment.pdf"))
        pdf(output_file3, width = 12, height = 8)
        simplifyGO(mat2, plot = TRUE)
        dev.off()
        message("已生成CC simplifyEnrichment结果: ", output_file3)
      }, error = function(e) {
        message("CC GO简化失败: ", e$message)
      })
    }

    if (length(go_id_mf) > 1) {
      tryCatch({
        mat3 <- GO_similarity(go_id_mf, ont = "MF", db = org_db, measure = "Sim_Resnik_1999")
        output_file4 <- file.path(output_folder, paste0(base_name, "_GO_MF_simplifyEnrichment.pdf"))
        pdf(output_file4, width = 12, height = 8)
        simplifyGO(mat3, plot = TRUE)
        dev.off()
        message("已生成MF simplifyEnrichment结果: ", output_file4)
      }, error = function(e) {
        message("MF GO简化失败: ", e$message)
      })
    }

    ego_df <- ego
    if (nrow(ego_df) > 0) {
      scores <- setNames(-log10(ego_df$qvalue), ego_df$ID)

      if ("BP" %in% ego_df$ONTOLOGY) {
        bp_ids <- ego_df$ID[ego_df$ONTOLOGY == "BP"]
        if (length(bp_ids) > 1) {
          tryCatch({
            simMatrix1 <- calculateSimMatrix(bp_ids, orgdb = org_db, ont = "BP", method = "Resnik")
            reducedTerms1 <- reduceSimMatrix(simMatrix1, scores, threshold = 0.8, orgdb = org_db)
            write.table(reducedTerms1, file = file.path(output_folder, paste0(base_name, "_GO_BP_rrvgo.txt")),
                        sep = "\t", row.names = FALSE)

            if (length(bp_ids) > 2) {
              plot1 <- scatterPlot(simMatrix1, reducedTerms1, algorithm = "pca")
              ggsave(file.path(output_folder, paste0(base_name, "_GO_BP_rrvgo.pdf")),
                     plot = plot1, width = 12, height = 12)
              message("已生成BP rrvgo散点图")
            }
          }, error = function(e) {
            message("BP rrvgo失败: ", e$message)
          })
        }
      }

      if ("CC" %in% ego_df$ONTOLOGY) {
        cc_ids <- ego_df$ID[ego_df$ONTOLOGY == "CC"]
        if (length(cc_ids) > 1) {
          tryCatch({
            simMatrix2 <- calculateSimMatrix(cc_ids, orgdb = org_db, ont = "CC", method = "Resnik")
            reducedTerms2 <- reduceSimMatrix(simMatrix2, scores, threshold = 0.8, orgdb = org_db)
            write.table(reducedTerms2, file = file.path(output_folder, paste0(base_name, "_GO_CC_rrvgo.txt")),
                        sep = "\t", row.names = FALSE)

            if (length(cc_ids) > 2) {
              plot2 <- scatterPlot(simMatrix2, reducedTerms2, algorithm = "pca")
              ggsave(file.path(output_folder, paste0(base_name, "_GO_CC_rrvgo.pdf")),
                     plot = plot2, width = 12, height = 12)
              message("已生成CC rrvgo散点图")
            }
          }, error = function(e) {
            message("CC rrvgo失败: ", e$message)
          })
        }
      }

      if ("MF" %in% ego_df$ONTOLOGY) {
        mf_ids <- ego_df$ID[ego_df$ONTOLOGY == "MF"]
        if (length(mf_ids) > 1) {
          tryCatch({
            simMatrix3 <- calculateSimMatrix(mf_ids, orgdb = org_db, ont = "MF", method = "Resnik")
            reducedTerms3 <- reduceSimMatrix(simMatrix3, scores, threshold = 0.8, orgdb = org_db)
            write.table(reducedTerms3, file = file.path(output_folder, paste0(base_name, "_GO_MF_rrvgo.txt")),
                        sep = "\t", row.names = FALSE)

            if (length(mf_ids) > 2) {
              plot3 <- scatterPlot(simMatrix3, reducedTerms3, algorithm = "pca")
              ggsave(file.path(output_folder, paste0(base_name, "_GO_MF_rrvgo.pdf")),
                     plot = plot3, width = 12, height = 12)
              message("已生成MF rrvgo散点图")
            }
          }, error = function(e) {
            message("MF rrvgo失败: ", e$message)
          })
        }
      }
    }

    message("\n", rep("-", 50))
    message("完成执行文件: ", basename(file))
    message(rep("-", 50), "\n")
  }
}

run_go_analysis(opt$`input-path`, opt$species, opt$`output-folder`)
