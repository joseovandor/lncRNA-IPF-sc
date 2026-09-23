# ============================================================
# 14_compartment_specificity.R
#
# lncRNA-IPF-sc | Step 14: TAR specificity across broad cell compartments
# Author: Jose A. Ovando-Ricardez
#
# Collapses the manual RNA annotation into 4 broad compartments
# (Epithelial, Immune, Mesenchymal, Endothelial) and, for every TAR in the
# 3/3 IPF sets (Known_GENCODE, Known_external, Unannotated), computes the
# per-sample and per-compartment detection, the dominant compartment and a
# descriptive IPF vs Control comparison. Also reports the previous priority
# loci (SMAD4, ADNP, DYRK1A, TCF4, PHLPP1, YES1, NUP50).
#
# Note: the IPF vs Control comparison is descriptive (3 IPF vs 1 Control).
#
# Inputs:  data/single_cell_RNA_Habermann/IPF_pilot_RNA_manual_compartments_res010.rds
#            (or IPF_pilot_manual_compartments_res010.rds)
#          data/single_cell_TAR/IPFTAR_consensus_single_cell_matrices.rds (Step 13)
#          data/TAR_catalog/quantification_3of3_IPF/{Known_GENCODE_lncRNA_3of3_IPF.csv,
#            Known_external_lncRNA_3of3_IPF.csv, Unannotated_lncRNA_candidates_3of3_IPF.csv,
#            Unannotated_lncRNA_candidates_3of3_IPF_0of1_Control.csv}
# Outputs: data/single_cell_TAR/broad_compartment_analysis/,
#          data/plots/07_TAR_broad_compartments/
# ============================================================


# ------------------------------------------------------------
# 1. Clean environment
# ------------------------------------------------------------

rm(list = ls())
gc()

# Sourced after rm(list = ls()) so the config variables are not removed
source("config/config.R")


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})


# ------------------------------------------------------------
# 3. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

DATA_DIR    <- file.path(BASE_DIR, "data")
RNA_DIR     <- file.path(DATA_DIR, "single_cell_RNA_Habermann")
TAR_DIR     <- file.path(DATA_DIR, "single_cell_TAR")
CATALOG_DIR <- file.path(DATA_DIR, "TAR_catalog")
QUANT_DIR   <- file.path(CATALOG_DIR, "quantification_3of3_IPF")
OUTPUT_DIR  <- file.path(TAR_DIR, "broad_compartment_analysis")
PLOT_DIR    <- file.path(DATA_DIR, "plots", "07_TAR_broad_compartments")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 4. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")
ipf_samples <- c("VUILD53", "VUILD63", "VUILD64")
control_sample <- "VUHD67"

broad_levels <- c("Epithelial", "Immune", "Mesenchymal", "Endothelial")


# ------------------------------------------------------------
# 5. Load RNA object
# ------------------------------------------------------------

rna_candidates <- c(
  file.path(RNA_DIR, "IPF_pilot_RNA_manual_compartments_res010.rds"),
  file.path(RNA_DIR, "IPF_pilot_manual_compartments_res010.rds")
)

rna_file <- rna_candidates[file.exists(rna_candidates)]

if (length(rna_file) == 0) {
  stop("Annotated RNA object not found.")
}

rna_file <- rna_file[1]

cat("\n")
cat("====================================================\n")
cat("LOADING RNA OBJECT\n")
cat("====================================================\n\n")

cat(rna_file, "\n")

combined <- readRDS(rna_file)

cat("Cells:", ncol(combined), "\n")


# ------------------------------------------------------------
# 6. Recover cell annotation
# ------------------------------------------------------------

metadata_columns <- colnames(combined@meta.data)

if (!"cell_compartment" %in% metadata_columns) {
  if ("cell_compartment_manual" %in% metadata_columns) {
    combined$cell_compartment <- combined$cell_compartment_manual
  } else {
    stop("The object does not contain cell_compartment.")
  }
}


# ------------------------------------------------------------
# 7. Collapse to 4 broad compartments
# ------------------------------------------------------------

old_compartment <- as.character(combined$cell_compartment)

combined$broad_compartment <- dplyr::case_when(
  old_compartment %in% c("Airway_Epithelial", "Alveolar_Epithelial") ~ "Epithelial",
  old_compartment %in% c("Myeloid", "Lymphoid") ~ "Immune",
  old_compartment == "Mesenchymal" ~ "Mesenchymal",
  old_compartment == "Endothelial" ~ "Endothelial",
  TRUE ~ NA_character_
)

combined$broad_compartment <- factor(combined$broad_compartment, levels = broad_levels)

if (any(is.na(combined$broad_compartment))) {
  stop("Some cells have no broad_compartment.")
}


# ------------------------------------------------------------
# 8. Cell summary
# ------------------------------------------------------------

broad_cell_counts <- combined@meta.data %>%
  dplyr::count(sample, condition, broad_compartment, name = "n_cells")

cat("\n")
cat("====================================================\n")
cat("CELLS IN THE 4 COMPARTMENTS\n")
cat("====================================================\n\n")

print(tibble::as_tibble(broad_cell_counts), n = Inf)

write.csv(broad_cell_counts, file.path(OUTPUT_DIR, "Broad_compartment_cell_counts.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 9. Load consensus TAR matrices (Step 13)
# ------------------------------------------------------------

matrix_file <- file.path(TAR_DIR, "IPFTAR_consensus_single_cell_matrices.rds")

if (!file.exists(matrix_file)) {
  stop(paste("File not found:", matrix_file))
}

tar_matrices <- readRDS(matrix_file)

cat("\n")
cat("Available matrices:\n")
print(names(tar_matrices))

missing_matrices <- setdiff(samples, names(tar_matrices))

if (length(missing_matrices) > 0) {
  stop(paste("Missing TAR matrices for:", paste(missing_matrices, collapse = ", ")))
}


# ------------------------------------------------------------
# 10. Load the 3 lncRNA groups
# ------------------------------------------------------------

known_gencode <- read.csv(file.path(QUANT_DIR, "Known_GENCODE_lncRNA_3of3_IPF.csv"),
                          stringsAsFactors = FALSE)

known_external <- read.csv(file.path(QUANT_DIR, "Known_external_lncRNA_3of3_IPF.csv"),
                           stringsAsFactors = FALSE)

unannotated <- read.csv(file.path(QUANT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF.csv"),
                        stringsAsFactors = FALSE)

unannotated_373 <- read.csv(file.path(QUANT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF_0of1_Control.csv"),
                            stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 11. Analysis universe
# ------------------------------------------------------------

candidate_metadata <- dplyr::bind_rows(
  known_gencode %>% dplyr::mutate(analysis_group = "Known_GENCODE"),
  known_external %>% dplyr::mutate(analysis_group = "Known_external"),
  unannotated %>% dplyr::mutate(analysis_group = "Unannotated")
)

candidate_metadata <- candidate_metadata %>%
  dplyr::distinct(TAR_ID, .keep_all = TRUE)

candidate_metadata$previous_0_control <- candidate_metadata$TAR_ID %in% unannotated_373$TAR_ID

candidate_ids <- candidate_metadata$TAR_ID

cat("\n")
cat("====================================================\n")
cat("lncRNA/TAR UNIVERSE\n")
cat("====================================================\n\n")

print(table(candidate_metadata$analysis_group))

cat("\nTotal:", length(candidate_ids), "\n")


# ------------------------------------------------------------
# 12. Check IDs in matrices
# ------------------------------------------------------------

for (sample_name in samples) {
  mat <- tar_matrices[[sample_name]]
  missing_ids <- setdiff(candidate_ids, rownames(mat))
  if (length(missing_ids) > 0) {
    stop(paste(sample_name, "is missing", length(missing_ids), "candidate_ids."))
  }
}


# ------------------------------------------------------------
# 13. TAR x sample x compartment metrics
# ------------------------------------------------------------

result_list <- list()

for (sample_name in samples) {
  cat("\n")
  cat("====================================================\n")
  cat("Processing:", sample_name, "\n")
  cat("====================================================\n")

  mat <- tar_matrices[[sample_name]]

  cat("Matrix:", nrow(mat), "TAR x", ncol(mat), "cells\n")

  # RNA metadata for this sample
  sample_meta <- combined@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(sample == sample_name)

  # Match matrix columns to metadata
  idx <- match(colnames(mat), sample_meta$cell)

  if (any(is.na(idx))) {
    stop(paste("TAR cells without RNA metadata in", sample_name))
  }

  sample_meta <- sample_meta[idx, , drop = FALSE]

  if (!identical(as.character(sample_meta$cell), as.character(colnames(mat)))) {
    stop(paste("Metadata/TAR order mismatch in", sample_name))
  }

  # Keep the 4,881 candidates
  mat <- mat[candidate_ids, , drop = FALSE]

  for (compartment_name in broad_levels) {
    cell_indices <- which(as.character(sample_meta$broad_compartment) == compartment_name)

    n_cells_compartment <- length(cell_indices)

    cat(compartment_name, ":", n_cells_compartment, "cells\n")

    if (n_cells_compartment == 0) {
      next
    }

    sub_mat <- mat[, cell_indices, drop = FALSE]

    total_counts <- Matrix::rowSums(sub_mat)
    cells_positive <- Matrix::rowSums(sub_mat > 0)
    pct_positive <- 100 * cells_positive / n_cells_compartment
    mean_counts_all <- total_counts / n_cells_compartment
    mean_counts_positive <- ifelse(cells_positive > 0, total_counts / cells_positive, 0)

    result_name <- paste(sample_name, compartment_name, sep = "__")

    result_list[[result_name]] <- data.frame(
      TAR_ID = candidate_ids,
      sample = sample_name,
      condition = ifelse(sample_name %in% ipf_samples, "IPF", "Control"),
      broad_compartment = compartment_name,
      n_cells = n_cells_compartment,
      cells_positive = as.numeric(cells_positive),
      pct_positive = as.numeric(pct_positive),
      total_counts = as.numeric(total_counts),
      mean_counts_all = as.numeric(mean_counts_all),
      mean_counts_positive = as.numeric(mean_counts_positive),
      stringsAsFactors = FALSE
    )
  }
}


# ------------------------------------------------------------
# 14. Combine results
# ------------------------------------------------------------

tar_broad_long <- dplyr::bind_rows(result_list)

tar_broad_long <- tar_broad_long %>%
  dplyr::left_join(candidate_metadata %>% dplyr::select(TAR_ID, analysis_group, previous_0_control),
                   by = "TAR_ID")

write.csv(tar_broad_long, file.path(OUTPUT_DIR, "IPFTAR_by_sample_4compartments.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 15. IPF summary
# ------------------------------------------------------------
# mean_pct_positive_IPF:   mean of the per-sample % positive cells (3 IPF samples)
# pooled_pct_positive_IPF: total positive cells / total cells
# The mean prevents the sample with the most cells from dominating.

ipf_summary <- tar_broad_long %>%
  dplyr::filter(condition == "IPF") %>%
  dplyr::group_by(TAR_ID, analysis_group, previous_0_control, broad_compartment) %>%
  dplyr::summarise(
    IPF_samples = dplyr::n_distinct(sample),
    IPF_samples_positive = sum(cells_positive > 0),
    total_cells_IPF = sum(n_cells),
    positive_cells_IPF = sum(cells_positive),
    pooled_pct_positive_IPF = 100 * positive_cells_IPF / total_cells_IPF,
    mean_pct_positive_IPF = mean(pct_positive, na.rm = TRUE),
    median_pct_positive_IPF = median(pct_positive, na.rm = TRUE),
    min_pct_positive_IPF = min(pct_positive, na.rm = TRUE),
    max_pct_positive_IPF = max(pct_positive, na.rm = TRUE),
    total_counts_IPF = sum(total_counts, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 16. Control summary
# ------------------------------------------------------------

control_summary <- tar_broad_long %>%
  dplyr::filter(condition == "Control") %>%
  dplyr::select(TAR_ID, broad_compartment, control_cells = n_cells, control_positive_cells = cells_positive,
                control_pct_positive = pct_positive, control_total_counts = total_counts,
                control_mean_counts_all = mean_counts_all)


# ------------------------------------------------------------
# 17. IPF vs Control
# ------------------------------------------------------------

comparison <- ipf_summary %>%
  dplyr::left_join(control_summary, by = c("TAR_ID", "broad_compartment")) %>%
  dplyr::mutate(
    delta_pct_IPF_vs_Control = mean_pct_positive_IPF - control_pct_positive,
    log2_prevalence_ratio = log2((mean_pct_positive_IPF + 0.1) / (control_pct_positive + 0.1)),
    present_3of3_in_compartment = IPF_samples_positive == 3,
    absent_control_in_compartment = control_positive_cells == 0
  )

write.csv(comparison, file.path(OUTPUT_DIR, "IPFTAR_4compartment_IPF_vs_Control.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 18. Dominant compartment
# ------------------------------------------------------------

specificity_ranked <- ipf_summary %>%
  dplyr::group_by(TAR_ID) %>%
  dplyr::arrange(dplyr::desc(mean_pct_positive_IPF), .by_group = TRUE) %>%
  dplyr::mutate(compartment_rank = dplyr::row_number())

specificity_top <- specificity_ranked %>%
  dplyr::filter(compartment_rank <= 2) %>%
  dplyr::select(TAR_ID, compartment_rank, broad_compartment, mean_pct_positive_IPF) %>%
  tidyr::pivot_wider(names_from = compartment_rank, values_from = c(broad_compartment, mean_pct_positive_IPF),
                     names_glue = "{.value}_{compartment_rank}")

specificity <- specificity_top %>%
  dplyr::transmute(
    TAR_ID = TAR_ID,
    dominant_compartment = broad_compartment_1,
    dominant_mean_pct = mean_pct_positive_IPF_1,
    second_compartment = broad_compartment_2,
    second_mean_pct = mean_pct_positive_IPF_2,
    specificity_margin = mean_pct_positive_IPF_1 - mean_pct_positive_IPF_2
  ) %>%
  dplyr::left_join(candidate_metadata %>% dplyr::select(TAR_ID, analysis_group, previous_0_control),
                   by = "TAR_ID")

write.csv(specificity, file.path(OUTPUT_DIR, "IPFTAR_4compartment_specificity.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 19. Dominant compartment summary
# ------------------------------------------------------------

dominant_summary <- specificity %>%
  dplyr::count(analysis_group, dominant_compartment, name = "n_TAR") %>%
  dplyr::group_by(analysis_group) %>%
  dplyr::mutate(percentage = 100 * n_TAR / sum(n_TAR)) %>%
  dplyr::ungroup()

write.csv(dominant_summary, file.path(OUTPUT_DIR, "Dominant_4compartment_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 20. Previous priority loci
# ------------------------------------------------------------

interesting_loci <- data.frame(
  TAR_ID = c("IPFTAR032288",
             "IPFTAR049318",
             "IPFTAR053240", "IPFTAR053241", "IPFTAR053242",
             "IPFTAR032366", "IPFTAR032367", "IPFTAR032368", "IPFTAR032369",
             "IPFTAR034005", "IPFTAR031406", "IPFTAR056555"),
  locus = c("SMAD4_near_intergenic",
            "ADNP_antisense",
            "DYRK1A_antisense", "DYRK1A_antisense", "DYRK1A_antisense",
            "TCF4_antisense", "TCF4_antisense", "TCF4_antisense", "TCF4_antisense",
            "PHLPP1_antisense", "YES1_antisense", "NUP50_antisense"),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 21. Check priority loci
# ------------------------------------------------------------

missing_interesting <- setdiff(interesting_loci$TAR_ID, candidate_ids)

if (length(missing_interesting) > 0) {
  warning(paste("These loci are not in candidate_ids:", paste(missing_interesting, collapse = ", ")))
}


# ------------------------------------------------------------
# 22. Priority loci results
# ------------------------------------------------------------

interesting_results <- comparison %>%
  dplyr::inner_join(interesting_loci, by = "TAR_ID") %>%
  dplyr::select(locus, TAR_ID, analysis_group, broad_compartment, IPF_samples_positive,
                total_cells_IPF, positive_cells_IPF, mean_pct_positive_IPF, pooled_pct_positive_IPF,
                control_cells, control_positive_cells, control_pct_positive, delta_pct_IPF_vs_Control,
                present_3of3_in_compartment, absent_control_in_compartment) %>%
  dplyr::arrange(locus, dplyr::desc(mean_pct_positive_IPF))

write.csv(interesting_results, file.path(OUTPUT_DIR, "Previous_priority_loci_4compartment.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 23. Main compartment of each priority locus
# ------------------------------------------------------------

interesting_dominant <- interesting_results %>%
  dplyr::group_by(locus, TAR_ID) %>%
  dplyr::slice_max(order_by = mean_pct_positive_IPF, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()


# ------------------------------------------------------------
# 24. Known lncRNAs present 3/3 in the same compartment
# ------------------------------------------------------------

known_3of3_compartment <- comparison %>%
  dplyr::filter(analysis_group %in% c("Known_GENCODE", "Known_external"), present_3of3_in_compartment)


# ------------------------------------------------------------
# 25. Known IPF > Control
# ------------------------------------------------------------

known_IPF_enriched <- known_3of3_compartment %>%
  dplyr::filter(delta_pct_IPF_vs_Control > 0) %>%
  dplyr::arrange(dplyr::desc(delta_pct_IPF_vs_Control))


# ------------------------------------------------------------
# 26. Known 3/3 IPF with 0 positive Control cells
# ------------------------------------------------------------

known_3of3_no_control <- known_3of3_compartment %>%
  dplyr::filter(absent_control_in_compartment) %>%
  dplyr::arrange(dplyr::desc(mean_pct_positive_IPF))


# ------------------------------------------------------------
# 27. Known summary
# ------------------------------------------------------------

known_summary <- known_3of3_compartment %>%
  dplyr::group_by(analysis_group, broad_compartment) %>%
  dplyr::summarise(
    n_3of3 = dplyr::n(),
    n_IPF_gt_control = sum(delta_pct_IPF_vs_Control > 0, na.rm = TRUE),
    n_control_zero = sum(absent_control_in_compartment, na.rm = TRUE),
    median_IPF_pct = median(mean_pct_positive_IPF, na.rm = TRUE),
    median_control_pct = median(control_pct_positive, na.rm = TRUE),
    median_delta = median(delta_pct_IPF_vs_Control, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(known_3of3_compartment, file.path(OUTPUT_DIR, "Known_lncRNA_3of3_IPF_by_compartment.csv"),
          row.names = FALSE)

write.csv(known_3of3_no_control, file.path(OUTPUT_DIR, "Known_lncRNA_3of3_IPF_0Control_same_compartment.csv"),
          row.names = FALSE)

write.csv(known_summary, file.path(OUTPUT_DIR, "Known_lncRNA_3of3_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 28. Top 10 known per compartment
# ------------------------------------------------------------

top10_known <- known_IPF_enriched %>%
  dplyr::group_by(broad_compartment) %>%
  dplyr::arrange(dplyr::desc(delta_pct_IPF_vs_Control), dplyr::desc(mean_pct_positive_IPF), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

write.csv(top10_known, file.path(OUTPUT_DIR, "Top10_known_lncRNA_IPF_enriched_per_compartment.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 29. Unannotated 3/3 in the same compartment
# ------------------------------------------------------------

unannotated_3of3_compartment <- comparison %>%
  dplyr::filter(analysis_group == "Unannotated", present_3of3_in_compartment)


# ------------------------------------------------------------
# 30. Top 10 unannotated per compartment
# ------------------------------------------------------------

top10_unannotated <- unannotated_3of3_compartment %>%
  dplyr::group_by(broad_compartment) %>%
  dplyr::arrange(dplyr::desc(delta_pct_IPF_vs_Control), dplyr::desc(mean_pct_positive_IPF), .by_group = TRUE) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

write.csv(top10_unannotated, file.path(OUTPUT_DIR, "Top10_unannotated_IPF_enriched_per_compartment.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 31. Subset of 373 (unannotated, 0/1 Control)
# ------------------------------------------------------------

subset_373 <- comparison %>%
  dplyr::filter(previous_0_control)

write.csv(subset_373, file.path(OUTPUT_DIR, "Unannotated_373_4compartment_context.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 32. Figure: cell numbers
# ------------------------------------------------------------

p_cells <- ggplot(broad_cell_counts, aes(x = sample, y = n_cells, fill = broad_compartment)) +
  geom_col() +
  theme_classic() +
  xlab(NULL) +
  ylab("Cells") +
  ggtitle("Four broad cell compartments")

ggsave(filename = file.path(PLOT_DIR, "01_Broad_compartment_cell_counts.png"), plot = p_cells,
       width = 9, height = 6, dpi = 300)


# ------------------------------------------------------------
# 33. Heatmap: priority loci
# ------------------------------------------------------------

p_interesting <- ggplot(interesting_results,
                        aes(x = broad_compartment, y = paste0(locus, " | ", TAR_ID), fill = mean_pct_positive_IPF)) +
  geom_tile() +
  theme_classic() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Previous candidate TARs - mean % positive cells in IPF")

ggsave(filename = file.path(PLOT_DIR, "02_Previous_candidates_4compartments.png"), plot = p_interesting,
       width = 10, height = 8, dpi = 300)


# ------------------------------------------------------------
# 34. Save RNA object with 4 compartments
# ------------------------------------------------------------

saveRDS(combined, file.path(OUTPUT_DIR, "IPF_pilot_RNA_4_broad_compartments.rds"))


# ------------------------------------------------------------
# 35. Final summary
# ------------------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("STEP 14 SUMMARY\n")
cat("====================================================\n\n")

# A. Cells
cat("1. CELLS IN 4 COMPARTMENTS\n\n")

print(tibble::as_tibble(broad_cell_counts), n = Inf)

# B. Dominance
cat("\n2. DOMINANT COMPARTMENT BY CLASS\n\n")

print(tibble::as_tibble(dominant_summary), n = Inf, width = Inf)

# C. Previous loci
cat("\n3. PREVIOUS PRIORITY LOCI\n\n")

print(
  interesting_dominant %>%
    dplyr::select(locus, TAR_ID, analysis_group, broad_compartment, IPF_samples_positive,
                  positive_cells_IPF, total_cells_IPF, mean_pct_positive_IPF, pooled_pct_positive_IPF,
                  control_positive_cells, control_cells, control_pct_positive, delta_pct_IPF_vs_Control) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

# D. Known summary
cat("\n4. KNOWN lncRNA 3/3 IPF - SUMMARY\n\n")

print(tibble::as_tibble(known_summary), n = Inf, width = Inf)

# E. Top known
cat("\n5. TOP 10 KNOWN lncRNA IPF > CONTROL\n\n")

print(
  top10_known %>%
    dplyr::select(TAR_ID, analysis_group, broad_compartment, IPF_samples_positive, positive_cells_IPF,
                  total_cells_IPF, mean_pct_positive_IPF, control_positive_cells, control_cells,
                  control_pct_positive, delta_pct_IPF_vs_Control) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

# F. Top unannotated
cat("\n6. TOP 10 UNANNOTATED 3/3 IPF\n\n")

print(
  top10_unannotated %>%
    dplyr::select(TAR_ID, broad_compartment, IPF_samples_positive, positive_cells_IPF, total_cells_IPF,
                  mean_pct_positive_IPF, control_positive_cells, control_cells, control_pct_positive,
                  delta_pct_IPF_vs_Control, previous_0_control) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

# G. Known with 0 Control
cat("\n7. KNOWN 3/3 IPF + 0 CONTROL IN THE SAME COMPARTMENT\n\n")

cat("Total:", nrow(known_3of3_no_control), "\n\n")

print(
  known_3of3_no_control %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::select(TAR_ID, analysis_group, broad_compartment, IPF_samples_positive, positive_cells_IPF,
                  total_cells_IPF, mean_pct_positive_IPF, control_positive_cells, control_cells,
                  delta_pct_IPF_vs_Control) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

# H. Strongest specificity
cat("\n8. HIGHEST SPECIFICITY ACROSS COMPARTMENTS\n\n")

print(
  specificity %>%
    dplyr::arrange(dplyr::desc(specificity_margin)) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::select(TAR_ID, analysis_group, dominant_compartment, dominant_mean_pct, second_compartment,
                  second_mean_pct, specificity_margin, previous_0_control) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\n")
cat("====================================================\n")
cat("STEP 14 FINISHED\n")
cat("====================================================\n")
