setwd('~/project/01_evolution/convergent_evolution/07_result/20251029_BWhigh_vs_low_pairwise_species/01result/')
library(dplyr)
library(readr)
library(GenomicRanges)

#pos test
#chicken_pos_1 <- read_tsv('chicken/chicken_pairwise_pos_top5_regions.tsv')
#chicken_pos_1 <- chicken_pos_1[,c(1,2,3)]
#colnames(chicken_pos_1)
#chicken_pos_2 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS/chicken/chicken_BW_pos_DCMS_top_005_raw_regions.txt',sep = '\t')
#colnames(chicken_pos_2) <- colnames(chicken_pos_1)
#identical_intervals <- inner_join(chicken_pos_1, chicken_pos_2, by = colnames(chicken_pos_1))
#write_tsv(identical_intervals,'chicken_intersect_BW_pos_DCMS_top005_raw_regions.tsv', na = "")
#union_intervals <- union(chicken_pos_1, chicken_pos_2)
#write_tsv(union_intervals,'chicken_union_BW_pos_DCMS_top005_raw_regions.tsv', na = "")

####循环 pos
species_list <- c("chicken", "pig", "cattle", "sheep", "dog")
process_species_comparison <- function(species, base_dir = NULL) {
  if (is.null(base_dir)) {
    base_dir <- "~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS"
  }
  pos_file1 <- paste0(species, "/", species, "_pairwise_pos_top5_regions.tsv")
  pos_data1 <- read_tsv(pos_file1)[, c(1, 2, 3)]
  pos_file2 <- paste0(base_dir, "/", species, "/", species, "_BW_pos_DCMS_top_005_raw_regions.txt")
  pos_data2 <- read.table(pos_file2, sep = "\t")
  colnames(pos_data2) <- colnames(pos_data1)
  results <- list()
  results$intersect <- inner_join(pos_data1, pos_data2, by = colnames(pos_data1))
  results$union <- union(pos_data1, pos_data2)
  write_tsv(results$intersect, 
            paste0(species, "_intersect_BW_pos_DCMS_top005_raw_regions.tsv"), 
            na = "")
  write_tsv(results$union, 
            paste0(species, "_union_BW_pos_DCMS_top005_raw_regions.tsv"), 
            na = "")
  
  return(list(
    species = species,
    intersection_count = nrow(results$intersect),
    union_count = nrow(results$union)
  ))
}

results_summary <- lapply(species_list, process_species_comparison)

# 打印汇总信息
cat("\nProcessing Summary:\n")
for (result in results_summary) {
  cat(result$species, ": Intersection =", result$intersection_count, 
      ", Union =", result$union_count, "\n")
}


####neg
process_species_comparison <- function(species, base_dir = NULL) {
  if (is.null(base_dir)) {
    base_dir <- "~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS"
  }
  pos_file1 <- paste0(species, "/", species, "_pairwise_neg_top5_regions.tsv")
  pos_data1 <- read_tsv(pos_file1)[, c(1, 2, 3)]
  pos_file2 <- paste0(base_dir, "/", species, "/", species, "_BW_neg_DCMS_top_005_raw_regions.txt")
  pos_data2 <- read.table(pos_file2, sep = "\t")
  colnames(pos_data2) <- colnames(pos_data1)
  results <- list()
  results$intersect <- inner_join(pos_data1, pos_data2, by = colnames(pos_data1))
  results$union <- union(pos_data1, pos_data2)
  write_tsv(results$intersect, 
            paste0(species, "_intersect_BW_neg_DCMS_top005_raw_regions.tsv"), 
            na = "")
  write_tsv(results$union, 
            paste0(species, "_union_BW_neg_DCMS_top005_raw_regions.tsv"), 
            na = "")
  
  return(list(
    species = species,
    intersection_count = nrow(results$intersect),
    union_count = nrow(results$union)
  ))
}

results_summary <- lapply(species_list, process_species_comparison)

cat("\nProcessing Summary:\n")
for (result in results_summary) {
  cat(result$species, ": Intersection =", result$intersection_count, 
      ", Union =", result$union_count, "\n")
}
