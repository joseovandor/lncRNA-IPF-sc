# ============================================================
# 07_quantify_recurrent_ipf_loci.R
#
# lncRNA-IPF-sc | Step 7: quantification of loci present in 3/3 IPF
# Author: Jose A. Ovando-Ricardez
#
# Builds a master table of consensus loci detected in all three IPF
# samples, split into GENCODE_lncRNA, External_lncRNA and
# Unannotated_candidate, and summarises per-sample counts, positive
# cells and % positive cells. Loci are ranked descriptively
# (min IPF counts, min % positive cells, CV) to prioritise candidates.
#
# Samples: IPF VUILD53, VUILD63, VUILD64; Control VUHD67
#
# Note: the IPF/control log2 ratio is descriptive only, not a formal
# differential expression test.
#
# Inputs:  data/TAR_catalog/external_annotation/IPF_uTAR_consensus_external_annotated.csv
#          data/TAR_catalog/IPF_uTAR_sample_mapping.csv
# Outputs: data/TAR_catalog/quantification_3of3_IPF/,
#          data/plots/04_lncRNA/04_quantification_3of3_IPF/
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

CATALOG_FILE <- file.path(BASE_DIR, "data/TAR_catalog/external_annotation/IPF_uTAR_consensus_external_annotated.csv")
MAPPING_FILE <- file.path(BASE_DIR, "data/TAR_catalog/IPF_uTAR_sample_mapping.csv")
OUTPUT_DIR   <- file.path(BASE_DIR, "data/TAR_catalog/quantification_3of3_IPF")
PLOT_DIR     <- file.path(BASE_DIR, "data/plots/04_lncRNA/04_quantification_3of3_IPF")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Check input files
# ------------------------------------------------------------

if (!file.exists(CATALOG_FILE)) {
  stop(paste("File not found:", CATALOG_FILE))
}

if (!file.exists(MAPPING_FILE)) {
  stop(paste("File not found:", MAPPING_FILE))
}


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("QUANTIFICATION OF lncRNA RECURRENT IN 3/3 IPF\n")
cat("====================================================\n\n")

catalog <- read.csv(CATALOG_FILE, stringsAsFactors = FALSE, check.names = FALSE)

mapping <- read.csv(MAPPING_FILE, stringsAsFactors = FALSE, check.names = FALSE)

cat("Consensus loci:", nrow(catalog), "\n")
cat("Mapping rows:", nrow(mapping), "\n\n")


# ------------------------------------------------------------
# 4. Required columns
# ------------------------------------------------------------

required_catalog <- c("TAR_ID", "chr", "start", "end", "strand", "length", "annotation_class_external",
                      "n_samples", "n_IPF", "n_Control", "presence_class")

required_mapping <- c("TAR_ID", "sample", "condition", "original_feature", "total_counts",
                      "cells_detected", "pct_cells")

missing_catalog <- setdiff(required_catalog, colnames(catalog))

missing_mapping <- setdiff(required_mapping, colnames(mapping))

if (length(missing_catalog) > 0) {
  stop(paste("Missing columns in catalog:", paste(missing_catalog, collapse = ", ")))
}

if (length(missing_mapping) > 0) {
  stop(paste("Missing columns in mapping:", paste(missing_mapping, collapse = ", ")))
}


# ------------------------------------------------------------
# 5. Samples
# ------------------------------------------------------------

ipf_samples <- c("VUILD53", "VUILD63", "VUILD64")

control_samples <- c("VUHD67")

all_samples <- c(ipf_samples, control_samples)


# ------------------------------------------------------------
# 6. Samples present in the mapping
# ------------------------------------------------------------

cat("Samples found:\n")

print(sort(unique(mapping$sample)))

cat("\n")


# ------------------------------------------------------------
# 7. Aggregate original features to consensus level
#
# A consensus TAR can have more than one original_feature in a
# sample. Counts are summed. cells_detected is summed as well, but
# this can overestimate unique cells when several features of the
# same locus are detected in the same cell; it is used here as a
# provisional descriptive measure for the pilot.
# pct_cells is summarised conservatively as the maximum percentage
# observed among the original features of the locus.
# ------------------------------------------------------------

mapping_summary <- mapping %>%
  filter(sample %in% all_samples) %>%
  group_by(TAR_ID, sample, condition) %>%
  summarise(
    total_counts        = sum(total_counts, na.rm = TRUE),
    cells_detected_sum  = sum(cells_detected, na.rm = TRUE),
    pct_cells_max       = max(pct_cells, na.rm = TRUE),
    n_original_features = n_distinct(original_feature),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 8. Replace Inf in pct_cells_max (all-NA groups)
# ------------------------------------------------------------

mapping_summary$pct_cells_max[!is.finite(mapping_summary$pct_cells_max)] <- NA_real_


# ------------------------------------------------------------
# 9. Wide format
# ------------------------------------------------------------

wide_counts <- mapping_summary %>%
  select(TAR_ID, sample, total_counts) %>%
  pivot_wider(names_from = sample, values_from = total_counts, values_fill = 0, names_prefix = "counts_")

wide_cells <- mapping_summary %>%
  select(TAR_ID, sample, cells_detected_sum) %>%
  pivot_wider(names_from = sample, values_from = cells_detected_sum, values_fill = 0, names_prefix = "cells_")

wide_pct <- mapping_summary %>%
  select(TAR_ID, sample, pct_cells_max) %>%
  pivot_wider(names_from = sample, values_from = pct_cells_max, values_fill = 0, names_prefix = "pct_")

wide_features <- mapping_summary %>%
  select(TAR_ID, sample, n_original_features) %>%
  pivot_wider(names_from = sample, values_from = n_original_features, values_fill = 0,
              names_prefix = "features_")


# ------------------------------------------------------------
# 10. Join with catalogue
# ------------------------------------------------------------

master <- catalog %>%
  left_join(wide_counts, by = "TAR_ID") %>%
  left_join(wide_cells, by = "TAR_ID") %>%
  left_join(wide_pct, by = "TAR_ID") %>%
  left_join(wide_features, by = "TAR_ID")


# ------------------------------------------------------------
# 11. Replace numeric NA with 0
# ------------------------------------------------------------

numeric_prefixes <- c("counts_", "cells_", "pct_", "features_")

numeric_cols <- grep(paste0("^(", paste(numeric_prefixes, collapse = "|"), ")"), colnames(master), value = TRUE)

for (col in numeric_cols) {
  master[[col]][is.na(master[[col]])] <- 0
}


# ------------------------------------------------------------
# 12. Analysis groups
# ------------------------------------------------------------

master$analysis_group <- case_when(
  master$annotation_class_external == "GENCODE_lncRNA" ~ "Known_GENCODE",
  master$annotation_class_external == "External_lncRNA" ~ "Known_external",
  master$annotation_class_external == "Unannotated_candidate" ~ "Unannotated",
  TRUE ~ "Other"
)


# ------------------------------------------------------------
# 13. Loci present in 3/3 IPF
# ------------------------------------------------------------

master_3of3 <- master %>%
  filter(n_IPF == 3)

cat("Total loci present in 3/3 IPF:", nrow(master_3of3), "\n\n")


# ------------------------------------------------------------
# 14. Summary by class
# ------------------------------------------------------------

group_summary <- master_3of3 %>%
  count(analysis_group, annotation_class_external, name = "n_loci") %>%
  arrange(desc(n_loci))

cat("====================================================\n")
cat("3/3 IPF LOCI BY CLASS\n")
cat("====================================================\n")

print(group_summary)

cat("\n")


# ------------------------------------------------------------
# 15. Main subsets
# ------------------------------------------------------------

known_gencode_3of3 <- master_3of3 %>%
  filter(analysis_group == "Known_GENCODE")

known_external_3of3 <- master_3of3 %>%
  filter(analysis_group == "Known_external")

unannotated_3of3 <- master_3of3 %>%
  filter(analysis_group == "Unannotated", length >= 200)

unannotated_3of3_no_control <- unannotated_3of3 %>%
  filter(n_Control == 0)

unannotated_3of3_with_control <- unannotated_3of3 %>%
  filter(n_Control > 0)


# ------------------------------------------------------------
# 16. IPF expression metrics: required columns
# ------------------------------------------------------------

required_count_cols <- paste0("counts_", ipf_samples)

required_pct_cols <- paste0("pct_", ipf_samples)

missing_count_cols <- setdiff(required_count_cols, colnames(master_3of3))

missing_pct_cols <- setdiff(required_pct_cols, colnames(master_3of3))

if (length(missing_count_cols) > 0) {
  stop(paste("Missing count columns:", paste(missing_count_cols, collapse = ", ")))
}

if (length(missing_pct_cols) > 0) {
  stop(paste("Missing percentage columns:", paste(missing_pct_cols, collapse = ", ")))
}


# ------------------------------------------------------------
# 17. Mean / median / CV of counts across the 3 IPF samples
# ------------------------------------------------------------

ipf_count_matrix <- as.matrix(master_3of3[, required_count_cols, drop = FALSE])

master_3of3$IPF_mean_counts   <- rowMeans(ipf_count_matrix)
master_3of3$IPF_median_counts <- apply(ipf_count_matrix, 1, median)
master_3of3$IPF_min_counts    <- apply(ipf_count_matrix, 1, min)
master_3of3$IPF_max_counts    <- apply(ipf_count_matrix, 1, max)
master_3of3$IPF_sd_counts     <- apply(ipf_count_matrix, 1, sd)

master_3of3$IPF_CV_counts <- ifelse(master_3of3$IPF_mean_counts > 0,
                                    master_3of3$IPF_sd_counts / master_3of3$IPF_mean_counts, NA_real_)


# ------------------------------------------------------------
# 18. % positive cells in IPF
# ------------------------------------------------------------

ipf_pct_matrix <- as.matrix(master_3of3[, required_pct_cols, drop = FALSE])

master_3of3$IPF_mean_pct_cells <- rowMeans(ipf_pct_matrix)
master_3of3$IPF_min_pct_cells  <- apply(ipf_pct_matrix, 1, min)
master_3of3$IPF_max_pct_cells  <- apply(ipf_pct_matrix, 1, max)


# ------------------------------------------------------------
# 19. Control counts and %
# ------------------------------------------------------------

control_count_col <- "counts_VUHD67"

control_pct_col <- "pct_VUHD67"

if (!control_count_col %in% colnames(master_3of3)) {
  master_3of3[[control_count_col]] <- 0
}

if (!control_pct_col %in% colnames(master_3of3)) {
  master_3of3[[control_pct_col]] <- 0
}

master_3of3$control_detected <- master_3of3[[control_count_col]] > 0


# ------------------------------------------------------------
# 20. Descriptive IPF / control ratio
#     (descriptive only, not formal DE; pseudocount = 1)
# ------------------------------------------------------------

master_3of3$log2_ratio_IPF_vs_control <- log2((master_3of3$IPF_mean_counts + 1) /
                                                (master_3of3[[control_count_col]] + 1))


# ------------------------------------------------------------
# 21. Rebuild subsets with metrics
# ------------------------------------------------------------

known_gencode_3of3 <- master_3of3 %>%
  filter(analysis_group == "Known_GENCODE")

known_external_3of3 <- master_3of3 %>%
  filter(analysis_group == "Known_external")

unannotated_3of3 <- master_3of3 %>%
  filter(analysis_group == "Unannotated", length >= 200)

unannotated_3of3_no_control <- unannotated_3of3 %>%
  filter(n_Control == 0)


# ------------------------------------------------------------
# 22. Descriptive consistency ranking (not a statistical test)
#     Priority: higher min counts, higher min % positive cells,
#     lower CV across the 3 IPF samples
# ------------------------------------------------------------

unannotated_3of3_ranked <- unannotated_3of3 %>%
  arrange(desc(IPF_min_counts), desc(IPF_min_pct_cells), IPF_CV_counts)

unannotated_3of3_no_control_ranked <- unannotated_3of3_no_control %>%
  arrange(desc(IPF_min_counts), desc(IPF_min_pct_cells), IPF_CV_counts)

known_gencode_3of3_ranked <- known_gencode_3of3 %>%
  arrange(desc(IPF_min_counts), desc(IPF_min_pct_cells), IPF_CV_counts)

known_external_3of3_ranked <- known_external_3of3 %>%
  arrange(desc(IPF_min_counts), desc(IPF_min_pct_cells), IPF_CV_counts)


# ------------------------------------------------------------
# 23. Save master table
# ------------------------------------------------------------

write.csv(master_3of3, file.path(OUTPUT_DIR, "IPF_lncRNA_master_3of3_IPF.csv"), row.names = FALSE)
write.csv(group_summary, file.path(OUTPUT_DIR, "IPF_lncRNA_3of3_group_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 24. Save subsets
# ------------------------------------------------------------

write.csv(known_gencode_3of3_ranked, file.path(OUTPUT_DIR, "Known_GENCODE_lncRNA_3of3_IPF.csv"), row.names = FALSE)
write.csv(known_external_3of3_ranked, file.path(OUTPUT_DIR, "Known_external_lncRNA_3of3_IPF.csv"), row.names = FALSE)
write.csv(unannotated_3of3_ranked, file.path(OUTPUT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF.csv"),
          row.names = FALSE)
write.csv(unannotated_3of3_no_control_ranked,
          file.path(OUTPUT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF_0of1_Control.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 25. Results summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MAIN SUMMARY\n")
cat("====================================================\n\n")

cat("GENCODE lncRNA present in 3/3 IPF:", nrow(known_gencode_3of3), "\n")
cat("External lncRNA present in 3/3 IPF:", nrow(known_external_3of3), "\n")
cat("Unannotated >=200 bp present in 3/3 IPF:", nrow(unannotated_3of3), "\n")
cat("Unannotated >=200 bp in 3/3 IPF and 0/1 control:", nrow(unannotated_3of3_no_control), "\n\n")


# ------------------------------------------------------------
# 26. Expression of unannotated 3/3 candidates
# ------------------------------------------------------------

cat("====================================================\n")
cat("EXPRESSION - UNANNOTATED 3/3 IPF\n")
cat("====================================================\n")

cat("\nMean IPF counts:\n")
print(summary(unannotated_3of3$IPF_mean_counts))

cat("\nMinimum counts across the 3 IPF samples:\n")
print(summary(unannotated_3of3$IPF_min_counts))

cat("\nMean % positive cells in IPF:\n")
print(summary(unannotated_3of3$IPF_mean_pct_cells))

cat("\nMinimum % positive cells across the 3 IPF samples:\n")
print(summary(unannotated_3of3$IPF_min_pct_cells))

cat("\nCV of counts across IPF samples:\n")
print(summary(unannotated_3of3$IPF_CV_counts))


# ------------------------------------------------------------
# 27. Control
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("PRESENCE IN CONTROL\n")
cat("====================================================\n")

cat("Unannotated 3/3 also detected in control:", sum(unannotated_3of3$control_detected, na.rm = TRUE), "\n")
cat("Unannotated 3/3 not detected in control:", sum(!unannotated_3of3$control_detected, na.rm = TRUE), "\n")


# ------------------------------------------------------------
# 28. Top 30 candidates: 3/3 IPF and 0/1 control
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP 30 UNANNOTATED 3/3 IPF + 0/1 CONTROL\n")
cat("====================================================\n")

top_cols <- intersect(
  c("TAR_ID", "chr", "start", "end", "strand", "length",
    "counts_VUILD53", "counts_VUILD63", "counts_VUILD64", "counts_VUHD67",
    "pct_VUILD53", "pct_VUILD63", "pct_VUILD64", "pct_VUHD67",
    "IPF_mean_counts", "IPF_min_counts", "IPF_mean_pct_cells", "IPF_min_pct_cells", "IPF_CV_counts",
    "log2_ratio_IPF_vs_control"),
  colnames(unannotated_3of3_no_control_ranked)
)

print(head(unannotated_3of3_no_control_ranked[, top_cols, drop = FALSE], 30))


# ------------------------------------------------------------
# 29. Top 20 GENCODE lncRNA 3/3
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP 20 GENCODE lncRNA 3/3 IPF\n")
cat("====================================================\n")

print(head(known_gencode_3of3_ranked[, intersect(top_cols, colnames(known_gencode_3of3_ranked)), drop = FALSE], 20))


# ------------------------------------------------------------
# 30. Top 20 external lncRNA 3/3
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP 20 EXTERNAL lncRNA 3/3 IPF\n")
cat("====================================================\n")

external_cols <- intersect(
  c(top_cols, "LncBook_ID", "NONCODE_ID", "LNCipedia_ID", "FANTOM_ID", "n_external_lncRNA_sources"),
  colnames(known_external_3of3_ranked)
)

print(head(known_external_3of3_ranked[, external_cols, drop = FALSE], 20))


# ------------------------------------------------------------
# 31. Figure 1: number of 3/3 IPF loci by class
# ------------------------------------------------------------

plot_table <- table(master_3of3$analysis_group)

pdf(file.path(PLOT_DIR, "01_lncRNA_3of3_IPF_by_class.pdf"), width = 8, height = 6)
barplot(plot_table, las = 2, ylab = "Number of loci", main = "lncRNA / uTAR present in 3/3 IPF samples")
dev.off()

png(file.path(PLOT_DIR, "01_lncRNA_3of3_IPF_by_class.png"), width = 1800, height = 1400, res = 200)
barplot(plot_table, las = 2, ylab = "Number of loci", main = "lncRNA / uTAR present in 3/3 IPF samples")
dev.off()


# ------------------------------------------------------------
# 32. Figure 2: count distribution
# ------------------------------------------------------------

pdf(file.path(PLOT_DIR, "02_unannotated_3of3_IPF_mean_counts.pdf"), width = 8, height = 6)
hist(log10(unannotated_3of3$IPF_mean_counts + 1), breaks = 50, xlab = "log10(mean IPF counts + 1)",
     main = "Expression of unannotated candidates in 3/3 IPF")
dev.off()

png(file.path(PLOT_DIR, "02_unannotated_3of3_IPF_mean_counts.png"), width = 1800, height = 1400, res = 200)
hist(log10(unannotated_3of3$IPF_mean_counts + 1), breaks = 50, xlab = "log10(mean IPF counts + 1)",
     main = "Expression of unannotated candidates in 3/3 IPF")
dev.off()


# ------------------------------------------------------------
# 33. RDS
# ------------------------------------------------------------

saveRDS(master_3of3, file.path(OUTPUT_DIR, "IPF_lncRNA_master_3of3_IPF.rds"))


# ------------------------------------------------------------
# 34. Done
# ------------------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("3/3 IPF ANALYSIS FINISHED\n")
cat("====================================================\n\n")

cat("Tables saved to:\n", OUTPUT_DIR, "\n\n")

cat("Figures saved to:\n", PLOT_DIR, "\n")
