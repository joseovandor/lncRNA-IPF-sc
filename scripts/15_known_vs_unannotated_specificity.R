# ============================================================
# 15_known_vs_unannotated_specificity.R
#
# lncRNA-IPF-sc | Step 15: compartment specificity, known vs unannotated TARs
# Author: Jose A. Ovando-Ricardez
#
# Post-hoc analysis of the broad-compartment results (Step 14):
#   - recomputes compartment specificity per TAR (top vs second compartment
#     by mean % positive IPF cells),
#   - focuses on the Epithelial and Mesenchymal compartments,
#   - re-checks previously prioritised loci (SMAD4, ADNP, DYRK1A, TCF4,
#     PHLPP1, YES1, NUP50),
#   - ranks compartment-specific unannotated TARs and known lncRNAs
#     (recovering the known lncRNA name/ID), and builds a combined
#     known vs unannotated table and a descriptive shortlist.
#
# Note: the shortlist is descriptive, NOT a final candidate selection.
#
# Inputs:  data/single_cell_TAR/broad_compartment_analysis/IPFTAR_4compartment_IPF_vs_Control.csv
#          data/TAR_catalog/quantification_3of3_IPF/{Known_GENCODE_lncRNA_3of3_IPF.csv,
#            Known_external_lncRNA_3of3_IPF.csv, Unannotated_lncRNA_candidates_3of3_IPF.csv,
#            Unannotated_lncRNA_candidates_3of3_IPF_0of1_Control.csv}
# Outputs: data/single_cell_TAR/broad_compartment_analysis/posthoc_specificity/,
#          data/plots/08_TAR_posthoc_specificity/ (created; no plots written)
# ============================================================

# ------------------------------------------------------------
# 1. Clean environment
# ------------------------------------------------------------

rm(list = ls())
gc()

# Loaded after rm() so that PROJECT_DIR is not cleared
source("config/config.R")


# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
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
TAR_DIR     <- file.path(DATA_DIR, "single_cell_TAR")
BROAD_DIR   <- file.path(TAR_DIR, "broad_compartment_analysis")
CATALOG_DIR <- file.path(DATA_DIR, "TAR_catalog")
QUANT_DIR   <- file.path(CATALOG_DIR, "quantification_3of3_IPF")
OUTPUT_DIR  <- file.path(BROAD_DIR, "posthoc_specificity")
PLOT_DIR    <- file.path(DATA_DIR, "plots", "08_TAR_posthoc_specificity")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 4. Inputs from the broad-compartment analysis
# ------------------------------------------------------------

comparison_file  <- file.path(BROAD_DIR, "IPFTAR_4compartment_IPF_vs_Control.csv")
specificity_file <- file.path(BROAD_DIR, "IPFTAR_4compartment_specificity.csv")

if (!file.exists(comparison_file)) {
  stop(paste("File not found:", comparison_file))
}

comparison <- read.csv(comparison_file, stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 5. Rebuild ipf_summary from comparison
# ------------------------------------------------------------

ipf_summary <- comparison %>%
  dplyr::select(
    TAR_ID, analysis_group, previous_0_control, broad_compartment, IPF_samples,
    IPF_samples_positive, total_cells_IPF, positive_cells_IPF, pooled_pct_positive_IPF,
    mean_pct_positive_IPF, median_pct_positive_IPF, min_pct_positive_IPF,
    max_pct_positive_IPF, total_counts_IPF
  )


# ------------------------------------------------------------
# 6. Candidate metadata
# ------------------------------------------------------------

known_gencode <- read.csv(
  file.path(QUANT_DIR, "Known_GENCODE_lncRNA_3of3_IPF.csv"), stringsAsFactors = FALSE
)

known_external <- read.csv(
  file.path(QUANT_DIR, "Known_external_lncRNA_3of3_IPF.csv"), stringsAsFactors = FALSE
)

unannotated <- read.csv(
  file.path(QUANT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF.csv"), stringsAsFactors = FALSE
)

unannotated_373 <- read.csv(
  file.path(QUANT_DIR, "Unannotated_lncRNA_candidates_3of3_IPF_0of1_Control.csv"),
  stringsAsFactors = FALSE
)

candidate_metadata <- dplyr::bind_rows(
  known_gencode %>% dplyr::mutate(analysis_group = "Known_GENCODE"),
  known_external %>% dplyr::mutate(analysis_group = "Known_external"),
  unannotated %>% dplyr::mutate(analysis_group = "Unannotated")
) %>%
  dplyr::distinct(TAR_ID, .keep_all = TRUE)

candidate_metadata$previous_0_control <- candidate_metadata$TAR_ID %in% unannotated_373$TAR_ID


# ------------------------------------------------------------
# 7. Compartment specificity
# ------------------------------------------------------------
# Per TAR: compartment 1 = highest mean % positive IPF, compartment 2 = second highest.
# specificity_margin = pct_compartment_1 - pct_compartment_2

specificity_fixed <- ipf_summary %>%
  dplyr::ungroup() %>%
  dplyr::group_by(TAR_ID) %>%
  dplyr::arrange(dplyr::desc(mean_pct_positive_IPF), .by_group = TRUE) %>%
  dplyr::mutate(rank_compartment = dplyr::row_number()) %>%
  dplyr::filter(rank_compartment <= 2) %>%
  dplyr::select(
    TAR_ID, analysis_group, previous_0_control, rank_compartment, broad_compartment,
    mean_pct_positive_IPF
  ) %>%
  tidyr::pivot_wider(
    names_from = rank_compartment,
    values_from = c(broad_compartment, mean_pct_positive_IPF),
    names_glue = "{.value}_{rank_compartment}"
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(specificity_margin = mean_pct_positive_IPF_1 - mean_pct_positive_IPF_2)

write.csv(
  specificity_fixed, file.path(OUTPUT_DIR, "IPFTAR_specificity_fixed.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 8. Dominant compartment distribution
# ------------------------------------------------------------

dominant_summary_fixed <- specificity_fixed %>%
  dplyr::count(analysis_group, broad_compartment_1, name = "n_TAR") %>%
  dplyr::group_by(analysis_group) %>%
  dplyr::mutate(percentage = 100 * n_TAR / sum(n_TAR)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(dominant_compartment = broad_compartment_1)

write.csv(
  dominant_summary_fixed, file.path(OUTPUT_DIR, "Dominant_compartment_summary_fixed.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Previously prioritised loci
# ------------------------------------------------------------

interesting_loci <- data.frame(
  TAR_ID = c(
    "IPFTAR032288",
    "IPFTAR049318",
    "IPFTAR053240", "IPFTAR053241", "IPFTAR053242",
    "IPFTAR032366", "IPFTAR032367", "IPFTAR032368", "IPFTAR032369",
    "IPFTAR034005", "IPFTAR031406", "IPFTAR056555"
  ),
  locus = c(
    "SMAD4_near_intergenic",
    "ADNP_antisense",
    "DYRK1A_antisense", "DYRK1A_antisense", "DYRK1A_antisense",
    "TCF4_antisense", "TCF4_antisense", "TCF4_antisense", "TCF4_antisense",
    "PHLPP1_antisense", "YES1_antisense", "NUP50_antisense"
  ),
  stringsAsFactors = FALSE
)

interesting_specificity <- specificity_fixed %>%
  dplyr::filter(TAR_ID %in% interesting_loci$TAR_ID) %>%
  dplyr::left_join(interesting_loci, by = "TAR_ID") %>%
  dplyr::arrange(dplyr::desc(specificity_margin))

write.csv(
  interesting_specificity,
  file.path(OUTPUT_DIR, "Previous_priority_loci_specificity_fixed.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 10. Epithelial / Mesenchymal results
# ------------------------------------------------------------

epi_mes_comparison <- comparison %>%
  dplyr::filter(broad_compartment %in% c("Epithelial", "Mesenchymal"))

write.csv(
  epi_mes_comparison,
  file.path(OUTPUT_DIR, "IPFTAR_Epithelial_Mesenchymal_IPF_vs_Control.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 11. Unannotated 3/3 IPF - Epithelial
# ------------------------------------------------------------

unannotated_epi <- comparison %>%
  dplyr::filter(
    analysis_group == "Unannotated", broad_compartment == "Epithelial", IPF_samples_positive == 3
  ) %>%
  dplyr::arrange(dplyr::desc(delta_pct_IPF_vs_Control), dplyr::desc(mean_pct_positive_IPF))

write.csv(
  unannotated_epi, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_Epithelial.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 12. Unannotated 3/3 IPF - Mesenchymal
# ------------------------------------------------------------

unannotated_mes <- comparison %>%
  dplyr::filter(
    analysis_group == "Unannotated", broad_compartment == "Mesenchymal", IPF_samples_positive == 3
  ) %>%
  dplyr::arrange(dplyr::desc(delta_pct_IPF_vs_Control), dplyr::desc(mean_pct_positive_IPF))

write.csv(
  unannotated_mes, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_Mesenchymal.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 13. Unannotated, Epithelial-dominant
# ------------------------------------------------------------

unannotated_epi_specific <- specificity_fixed %>%
  dplyr::filter(analysis_group == "Unannotated", broad_compartment_1 == "Epithelial") %>%
  dplyr::arrange(dplyr::desc(specificity_margin))

write.csv(
  unannotated_epi_specific,
  file.path(OUTPUT_DIR, "Unannotated_specificity_Epithelial.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 14. Unannotated, Mesenchymal-dominant
# ------------------------------------------------------------

unannotated_mes_specific <- specificity_fixed %>%
  dplyr::filter(analysis_group == "Unannotated", broad_compartment_1 == "Mesenchymal") %>%
  dplyr::arrange(dplyr::desc(specificity_margin))

write.csv(
  unannotated_mes_specific,
  file.path(OUTPUT_DIR, "Unannotated_specificity_Mesenchymal.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 15. Known lncRNA metadata
# ------------------------------------------------------------

known_metadata <- candidate_metadata %>%
  dplyr::filter(analysis_group %in% c("Known_GENCODE", "Known_external"))


# ------------------------------------------------------------
# 16. Join specificity with known metadata
# ------------------------------------------------------------

known_specificity <- specificity_fixed %>%
  dplyr::filter(analysis_group %in% c("Known_GENCODE", "Known_external")) %>%
  dplyr::left_join(known_metadata, by = c("TAR_ID", "analysis_group", "previous_0_control"))


# ------------------------------------------------------------
# 17. Add IPF vs Control stats of the dominant compartment
# ------------------------------------------------------------

known_specificity <- known_specificity %>%
  dplyr::left_join(
    comparison %>%
      dplyr::select(
        TAR_ID, broad_compartment, IPF_samples_positive, positive_cells_IPF, total_cells_IPF,
        mean_pct_positive_IPF, pooled_pct_positive_IPF, control_positive_cells, control_cells,
        control_pct_positive, delta_pct_IPF_vs_Control, absent_control_in_compartment
      ),
    by = c("TAR_ID" = "TAR_ID", "broad_compartment_1" = "broad_compartment")
  )


# ------------------------------------------------------------
# 18. Known lncRNA name
# ------------------------------------------------------------
# Priority: GENCODE gene_name > LncBook ID > NONCODE ID > LNCipedia ID > FANTOM ID

known_named <- known_specificity %>%
  dplyr::mutate(
    known_name = dplyr::case_when(
      !is.na(GENCODE_lncRNA_gene_name) & GENCODE_lncRNA_gene_name != "" ~ GENCODE_lncRNA_gene_name,
      !is.na(LncBook_ID) & LncBook_ID != "" ~ LncBook_ID,
      !is.na(NONCODE_ID) & NONCODE_ID != "" ~ NONCODE_ID,
      !is.na(LNCipedia_ID) & LNCipedia_ID != "" ~ LNCipedia_ID,
      !is.na(FANTOM_ID) & FANTOM_ID != "" ~ FANTOM_ID,
      TRUE ~ NA_character_
    ),
    known_source = dplyr::case_when(
      !is.na(GENCODE_lncRNA_gene_name) & GENCODE_lncRNA_gene_name != "" ~ "GENCODE",
      !is.na(LncBook_ID) & LncBook_ID != "" ~ "LncBook",
      !is.na(NONCODE_ID) & NONCODE_ID != "" ~ "NONCODE",
      !is.na(LNCipedia_ID) & LNCipedia_ID != "" ~ "LNCipedia",
      !is.na(FANTOM_ID) & FANTOM_ID != "" ~ "FANTOM",
      TRUE ~ NA_character_
    )
  )


# ------------------------------------------------------------
# 19-22. Known lncRNAs by dominant compartment
# ------------------------------------------------------------

known_epithelial <- known_named %>%
  dplyr::filter(broad_compartment_1 == "Epithelial") %>%
  dplyr::arrange(dplyr::desc(specificity_margin))

known_mesenchymal <- known_named %>%
  dplyr::filter(broad_compartment_1 == "Mesenchymal") %>%
  dplyr::arrange(dplyr::desc(specificity_margin))

# 3/3 IPF and 0 positive control cells
known_epithelial_3of3_0control <- known_epithelial %>%
  dplyr::filter(IPF_samples_positive == 3, control_positive_cells == 0) %>%
  dplyr::arrange(dplyr::desc(specificity_margin), dplyr::desc(mean_pct_positive_IPF_1))

known_mesenchymal_3of3_0control <- known_mesenchymal %>%
  dplyr::filter(IPF_samples_positive == 3, control_positive_cells == 0) %>%
  dplyr::arrange(dplyr::desc(specificity_margin), dplyr::desc(mean_pct_positive_IPF_1))


# ------------------------------------------------------------
# 23. Save known lncRNA tables
# ------------------------------------------------------------

write.csv(
  known_epithelial, file.path(OUTPUT_DIR, "Known_lncRNA_specificity_Epithelial.csv"),
  row.names = FALSE
)
write.csv(
  known_mesenchymal, file.path(OUTPUT_DIR, "Known_lncRNA_specificity_Mesenchymal.csv"),
  row.names = FALSE
)
write.csv(
  known_epithelial_3of3_0control,
  file.path(OUTPUT_DIR, "Known_lncRNA_Epithelial_3of3_IPF_0Control.csv"), row.names = FALSE
)
write.csv(
  known_mesenchymal_3of3_0control,
  file.path(OUTPUT_DIR, "Known_lncRNA_Mesenchymal_3of3_IPF_0Control.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 24. Combined known vs unannotated table (Epithelial / Mesenchymal)
# ------------------------------------------------------------

known_final <- known_named %>%
  dplyr::filter(broad_compartment_1 %in% c("Epithelial", "Mesenchymal")) %>%
  dplyr::select(
    TAR_ID, known_name, known_source, analysis_group, chr, start, end, strand,
    broad_compartment_1, mean_pct_positive_IPF_1, broad_compartment_2, mean_pct_positive_IPF_2,
    specificity_margin, IPF_samples_positive, positive_cells_IPF, total_cells_IPF,
    control_positive_cells, control_cells, control_pct_positive, delta_pct_IPF_vs_Control
  )

unannotated_final <- specificity_fixed %>%
  dplyr::filter(
    analysis_group == "Unannotated", broad_compartment_1 %in% c("Epithelial", "Mesenchymal")
  ) %>%
  dplyr::left_join(
    comparison %>%
      dplyr::select(
        TAR_ID, broad_compartment, IPF_samples_positive, positive_cells_IPF, total_cells_IPF,
        control_positive_cells, control_cells, control_pct_positive, delta_pct_IPF_vs_Control
      ),
    by = c("TAR_ID" = "TAR_ID", "broad_compartment_1" = "broad_compartment")
  ) %>%
  dplyr::left_join(
    candidate_metadata %>% dplyr::select(TAR_ID, chr, start, end, strand),
    by = "TAR_ID"
  ) %>%
  dplyr::mutate(known_name = NA_character_, known_source = NA_character_) %>%
  dplyr::select(
    TAR_ID, known_name, known_source, analysis_group, chr, start, end, strand,
    broad_compartment_1, mean_pct_positive_IPF_1, broad_compartment_2, mean_pct_positive_IPF_2,
    specificity_margin, IPF_samples_positive, positive_cells_IPF, total_cells_IPF,
    control_positive_cells, control_cells, control_pct_positive, delta_pct_IPF_vs_Control
  )

final_comparison_table <- dplyr::bind_rows(known_final, unannotated_final)

write.csv(
  final_comparison_table,
  file.path(OUTPUT_DIR, "Known_vs_Unannotated_Epithelial_Mesenchymal.csv"), row.names = FALSE
)


# ------------------------------------------------------------
# 25. Descriptive shortlist (not a final selection)
# ------------------------------------------------------------
# 3/3 IPF, Epithelial or Mesenchymal dominant, delta IPF - Control > 0

shortlist <- final_comparison_table %>%
  dplyr::filter(IPF_samples_positive == 3, delta_pct_IPF_vs_Control > 0) %>%
  dplyr::arrange(
    broad_compartment_1, dplyr::desc(specificity_margin), dplyr::desc(delta_pct_IPF_vs_Control)
  )

write.csv(
  shortlist, file.path(OUTPUT_DIR, "Shortlist_Epithelial_Mesenchymal_3of3_IPF.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 15 SUMMARY\n")
cat("====================================================\n\n")

cat("Known epithelial-dominant:", nrow(known_epithelial), "\n")
cat("Known mesenchymal-dominant:", nrow(known_mesenchymal), "\n")
cat("Known epithelial 3/3 IPF + 0 control:", nrow(known_epithelial_3of3_0control), "\n")
cat("Known mesenchymal 3/3 IPF + 0 control:", nrow(known_mesenchymal_3of3_0control), "\n\n")

cat("Unannotated epithelial-dominant:", nrow(unannotated_epi_specific), "\n")
cat("Unannotated mesenchymal-dominant:", nrow(unannotated_mes_specific), "\n\n")

cat("TOP 20 KNOWN EPITHELIAL 3/3 IPF + 0 CONTROL\n\n")

print(
  known_epithelial_3of3_0control %>%
    dplyr::select(
      TAR_ID, known_name, known_source, chr, start, end, strand, mean_pct_positive_IPF_1,
      broad_compartment_2, mean_pct_positive_IPF_2, specificity_margin, positive_cells_IPF,
      control_pct_positive
    ) %>%
    dplyr::slice_head(n = 20) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\nTOP KNOWN MESENCHYMAL 3/3 IPF + 0 CONTROL\n\n")

print(
  known_mesenchymal_3of3_0control %>%
    dplyr::select(
      TAR_ID, known_name, known_source, chr, start, end, strand, mean_pct_positive_IPF_1,
      broad_compartment_2, mean_pct_positive_IPF_2, specificity_margin, positive_cells_IPF,
      control_pct_positive
    ) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\nTOP 20 UNANNOTATED EPITHELIAL\n\n")

print(
  shortlist %>%
    dplyr::filter(analysis_group == "Unannotated", broad_compartment_1 == "Epithelial") %>%
    dplyr::slice_head(n = 20) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\nTOP 20 UNANNOTATED MESENCHYMAL\n\n")

print(
  shortlist %>%
    dplyr::filter(analysis_group == "Unannotated", broad_compartment_1 == "Mesenchymal") %>%
    dplyr::slice_head(n = 20) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\nPREVIOUS LOCI OF INTEREST\n\n")

print(
  interesting_specificity %>%
    dplyr::select(
      locus, TAR_ID, broad_compartment_1, mean_pct_positive_IPF_1, broad_compartment_2,
      mean_pct_positive_IPF_2, specificity_margin, previous_0_control
    ) %>%
    tibble::as_tibble(),
  n = Inf, width = Inf
)

cat("\n")
cat("====================================================\n")
cat("STEP 15 FINISHED\n")
cat("====================================================\n")
