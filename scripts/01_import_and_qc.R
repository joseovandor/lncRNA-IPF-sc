# ============================================================
# 01_import_and_qc.R
#
# lncRNA-IPF-sc | Step 1: import Cell Ranger and TAR matrices and run QC
# Author: Jose A. Ovando-Ricardez
#
# - Reads the conventional Cell Ranger matrix and the TAR-scRNA-seq
#   feature-barcode matrix of every sample.
# - Checks barcode correspondence between both matrices.
# - Computes per-sample and per-cell RNA/TAR statistics and TAR sparsity.
# - Preliminary TAR classification from the feature-name suffix produced
#   by the TAR pipeline: "_0" = unannotated (uTAR candidate),
#   "_1" = annotated (aTAR). A "_0" feature is NOT yet a novel lncRNA.
# - No TAR is removed at this stage.
#
# Outputs: data/TAR_QC_summary.csv, data/TAR_cell_QC.csv,
#          data/plots/01_QC/*.png|pdf
# ============================================================

source("config/config.R")

library(Seurat)
library(Matrix)
library(ggplot2)


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

data_dir <- DATA_DIR

plots_dir            <- file.path(data_dir, "plots")
qc_plot_dir          <- file.path(plots_dir, "01_QC")
tar_catalog_plot_dir <- file.path(plots_dir, "02_TAR_catalog")
single_cell_plot_dir <- file.path(plots_dir, "03_single_cell")
lncrna_plot_dir      <- file.path(plots_dir, "04_lncRNA")

dir.create(qc_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tar_catalog_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(single_cell_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(lncrna_plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nDirectories created:\n")
cat(plots_dir, "\n")
cat(qc_plot_dir, "\n")
cat(tar_catalog_plot_dir, "\n")
cat(single_cell_plot_dir, "\n")
cat(lncrna_plot_dir, "\n")


# ------------------------------------------------------------
# 2. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")

condition <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)

sample_levels <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")


# ------------------------------------------------------------
# 3. Result containers
# ------------------------------------------------------------

qc_summary <- data.frame()
cell_qc_list <- list()


# ------------------------------------------------------------
# 4. Per-sample processing
# ------------------------------------------------------------

for (s in samples) {

  cat("\n")
  cat("====================================================\n")
  cat("Sample:", s, "\n")
  cat("Condition:", condition[s], "\n")
  cat("====================================================\n")

  # 4.1 Conventional Cell Ranger matrix
  cr <- Read10X(data.dir = file.path(data_dir, s, "cellranger"))

  # 4.2 TAR matrix
  tar <- Read10X(data.dir = file.path(data_dir, s, "TAR", "TAR_feature_bc_matrix"))

  # 4.3 Dimensions
  cat("\nCell Ranger:\n")
  cat("  Features:", nrow(cr), "\n")
  cat("  Cells   :", ncol(cr), "\n")

  cat("\nTAR:\n")
  cat("  Features:", nrow(tar), "\n")
  cat("  Cells   :", ncol(tar), "\n")

  # 4.4 Barcode correspondence
  shared_barcodes <- intersect(colnames(cr), colnames(tar))

  pct_cr_in_tar <- length(shared_barcodes) / ncol(cr) * 100
  pct_tar_in_cr <- length(shared_barcodes) / ncol(tar) * 100

  cat("\nShared barcodes:", length(shared_barcodes), "\n")
  cat("CR barcodes present in TAR:", round(pct_cr_in_tar, 2), "%\n")
  cat("TAR barcodes present in CR:", round(pct_tar_in_cr, 2), "%\n")

  # 4.5 Inspect TAR feature names
  cat("\nFirst 20 TAR features:\n")
  print(head(rownames(tar), 20))

  cat("\nLast 20 TAR features:\n")
  print(tail(rownames(tar), 20))

  # 4.6 Global counts per TAR
  tar_row_counts <- Matrix::rowSums(tar)

  tar_detected_gt0 <- sum(tar_row_counts > 0)
  tar_ge10         <- sum(tar_row_counts >= 10)
  tar_ge100        <- sum(tar_row_counts >= 100)

  cat("\nTAR detected (>0 counts):", tar_detected_gt0, "\n")
  cat("TAR with >= 10 counts:", tar_ge10, "\n")
  cat("TAR with >= 100 counts:", tar_ge100, "\n")

  # 4.7 RNA QC per cell
  rna_counts_cell   <- Matrix::colSums(cr)
  rna_features_cell <- Matrix::colSums(cr > 0)

  cat("\nRNA counts per cell:\n")
  print(summary(rna_counts_cell))

  cat("\nRNA features per cell:\n")
  print(summary(rna_features_cell))

  # 4.8 TAR QC per cell
  tar_counts_cell   <- Matrix::colSums(tar)
  tar_features_cell <- Matrix::colSums(tar > 0)

  cat("\nTAR counts per cell:\n")
  print(summary(tar_counts_cell))

  cat("\nTAR detected per cell:\n")
  print(summary(tar_features_cell))

  # 4.9 TAR sparsity
  nnz_tar             <- length(tar@x)
  total_positions_tar <- nrow(tar) * ncol(tar)
  sparsity_tar        <- 1 - (nnz_tar / total_positions_tar)

  cat("\nTAR sparsity:", round(sparsity_tar, 6), "\n")

  # 4.10 Preliminary classification from the feature suffix
  #      "_0" = no associated annotation (uTAR candidate)
  #      "_1" = associated annotation (aTAR)
  is_utar          <- grepl("_0$", rownames(tar))
  is_annotated_tar <- grepl("_1$", rownames(tar))

  n_utar          <- sum(is_utar)
  n_annotated_tar <- sum(is_annotated_tar)
  n_other         <- nrow(tar) - n_utar - n_annotated_tar

  cat("\nPreliminary classification:\n")
  cat("  uTAR candidates:", n_utar, "\n")
  cat("  Annotated TAR:", n_annotated_tar, "\n")
  cat("  Other/unclassified:", n_other, "\n")

  # 4.11 Per-cell table for plots
  cell_qc_list[[s]] <- data.frame(
    sample       = s,
    condition    = condition[s],
    barcode      = colnames(cr),
    RNA_counts   = as.numeric(rna_counts_cell),
    RNA_features = as.numeric(rna_features_cell),
    TAR_counts   = as.numeric(tar_counts_cell),
    TAR_features = as.numeric(tar_features_cell),
    stringsAsFactors = FALSE
  )

  # 4.12 Per-sample summary
  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample                   = s,
      condition                = condition[s],
      cells                    = ncol(cr),
      shared_barcodes          = length(shared_barcodes),
      pct_CR_in_TAR            = pct_cr_in_tar,
      pct_TAR_in_CR            = pct_tar_in_cr,
      RNA_features             = nrow(cr),
      RNA_total_counts         = sum(cr),
      RNA_median_counts_cell   = median(rna_counts_cell),
      RNA_mean_counts_cell     = mean(rna_counts_cell),
      RNA_median_features_cell = median(rna_features_cell),
      RNA_mean_features_cell   = mean(rna_features_cell),
      TAR_features             = nrow(tar),
      TAR_total_counts         = sum(tar),
      TAR_median_counts_cell   = median(tar_counts_cell),
      TAR_mean_counts_cell     = mean(tar_counts_cell),
      TAR_median_features_cell = median(tar_features_cell),
      TAR_mean_features_cell   = mean(tar_features_cell),
      TAR_detected_gt0         = tar_detected_gt0,
      TAR_ge10                 = tar_ge10,
      TAR_ge100                = tar_ge100,
      candidate_uTAR           = n_utar,
      annotated_TAR            = n_annotated_tar,
      other_TAR                = n_other,
      TAR_sparsity             = sparsity_tar,
      stringsAsFactors = FALSE
    )
  )

  # 4.13 Free memory
  rm(cr, tar, shared_barcodes, tar_row_counts, rna_counts_cell,
     rna_features_cell, tar_counts_cell, tar_features_cell)
  gc()
}


# ------------------------------------------------------------
# 5. Combine per-cell information
# ------------------------------------------------------------

cell_qc <- do.call(rbind, cell_qc_list)
rownames(cell_qc) <- NULL

cell_qc$sample    <- factor(cell_qc$sample, levels = sample_levels)
qc_summary$sample <- factor(qc_summary$sample, levels = sample_levels)


# ------------------------------------------------------------
# 6. Overall summary
# ------------------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("OVERALL SUMMARY\n")
cat("====================================================\n\n")

print(qc_summary)


# ------------------------------------------------------------
# 7. Save tables
# ------------------------------------------------------------

write.csv(qc_summary, file = file.path(data_dir, "TAR_QC_summary.csv"), row.names = FALSE)
write.csv(cell_qc, file = file.path(data_dir, "TAR_cell_QC.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 8. Helper: save PNG + PDF
# ------------------------------------------------------------

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(filename = file.path(qc_plot_dir, paste0(filename, ".png")),
         plot = plot_object, width = width, height = height, dpi = 300)
  ggsave(filename = file.path(qc_plot_dir, paste0(filename, ".pdf")),
         plot = plot_object, width = width, height = height)
}


# ------------------------------------------------------------
# 9. Plot 1: total TAR per sample
# ------------------------------------------------------------

p1 <- ggplot(qc_summary, aes(x = sample, y = TAR_features)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(TAR_features, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Total TAR detected per sample", x = "Sample", y = "Number of TAR") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)
save_plot(p1, "01_total_TAR_per_sample")


# ------------------------------------------------------------
# 10. Plot 2: uTAR vs annotated TAR
# ------------------------------------------------------------

tar_class_df <- rbind(
  data.frame(sample = qc_summary$sample, class = "uTAR candidate", count = qc_summary$candidate_uTAR),
  data.frame(sample = qc_summary$sample, class = "Annotated TAR", count = qc_summary$annotated_TAR)
)

tar_class_df$sample <- factor(tar_class_df$sample, levels = sample_levels)

p2 <- ggplot(tar_class_df, aes(x = sample, y = count, fill = class)) +
  geom_col(position = "stack", width = 0.7) +
  labs(title = "Preliminary composition of TAR regions", x = "Sample",
       y = "Number of TAR", fill = "Class") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
save_plot(p2, "02_uTAR_vs_annotated_TAR")


# ------------------------------------------------------------
# 11. Plot 3: median TAR counts per cell
# ------------------------------------------------------------

p3 <- ggplot(qc_summary, aes(x = sample, y = TAR_median_counts_cell)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(TAR_median_counts_cell, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Median TAR counts per cell", x = "Sample", y = "Median TAR counts") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)
save_plot(p3, "03_median_TAR_counts_per_cell")


# ------------------------------------------------------------
# 12. Plot 4: median TAR features per cell
# ------------------------------------------------------------

p4 <- ggplot(qc_summary, aes(x = sample, y = TAR_median_features_cell)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(TAR_median_features_cell, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Median TAR detected per cell", x = "Sample", y = "TAR detected per cell") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p4)
save_plot(p4, "04_median_TAR_features_per_cell")


# ------------------------------------------------------------
# 13. Plot 5: distribution of TAR counts per cell (log10, large
#     differences between samples)
# ------------------------------------------------------------

p5 <- ggplot(cell_qc, aes(x = sample, y = TAR_counts)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Distribution of TAR counts per cell", x = "Sample",
       y = "TAR counts per cell (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p5)
save_plot(p5, "05_distribution_TAR_counts_per_cell")


# ------------------------------------------------------------
# 14. Plot 6: distribution of TAR features per cell
# ------------------------------------------------------------

p6 <- ggplot(cell_qc, aes(x = sample, y = TAR_features)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Distribution of TAR detected per cell", x = "Sample",
       y = "TAR per cell (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p6)
save_plot(p6, "06_distribution_TAR_features_per_cell")


# ------------------------------------------------------------
# 15. Plot 7: RNA counts per cell
# ------------------------------------------------------------

p7 <- ggplot(cell_qc, aes(x = sample, y = RNA_counts)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Distribution of RNA counts per cell", x = "Sample",
       y = "RNA counts per cell (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p7)
save_plot(p7, "07_distribution_RNA_counts_per_cell")


# ------------------------------------------------------------
# 16. Plot 8: RNA features per cell
# ------------------------------------------------------------

p8 <- ggplot(cell_qc, aes(x = sample, y = RNA_features)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Distribution of genes detected per cell", x = "Sample",
       y = "Genes detected per cell (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p8)
save_plot(p8, "08_distribution_RNA_features_per_cell")


# ------------------------------------------------------------
# 17. Plot 9: TAR matrix sparsity (% of zero entries)
# ------------------------------------------------------------

qc_summary$TAR_sparsity_percent <- qc_summary$TAR_sparsity * 100

p9 <- ggplot(qc_summary, aes(x = sample, y = TAR_sparsity_percent)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(TAR_sparsity_percent, 2), "%")), vjust = -0.4, size = 4) +
  labs(title = "Sparsity of the TAR matrices", x = "Sample", y = "Zero entries (%)") +
  ylim(0, 105) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p9)
save_plot(p9, "09_TAR_sparsity")


# ------------------------------------------------------------
# 18. Plot 10: TAR with >=10 and >=100 counts
# ------------------------------------------------------------

tar_threshold_df <- rbind(
  data.frame(sample = qc_summary$sample, threshold = ">=10 counts", count = qc_summary$TAR_ge10),
  data.frame(sample = qc_summary$sample, threshold = ">=100 counts", count = qc_summary$TAR_ge100)
)

tar_threshold_df$sample <- factor(tar_threshold_df$sample, levels = sample_levels)

p10 <- ggplot(tar_threshold_df, aes(x = sample, y = count, fill = threshold)) +
  geom_col(position = "dodge", width = 0.7) +
  labs(title = "Expression support of TAR regions", x = "Sample",
       y = "Number of TAR", fill = "Threshold") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p10)
save_plot(p10, "10_TAR_expression_thresholds")


# ------------------------------------------------------------
# 19. Plot 11: median RNA vs TAR counts per cell
# ------------------------------------------------------------

rna_tar_median_df <- rbind(
  data.frame(sample = qc_summary$sample, assay = "RNA", median_counts = qc_summary$RNA_median_counts_cell),
  data.frame(sample = qc_summary$sample, assay = "TAR", median_counts = qc_summary$TAR_median_counts_cell)
)

rna_tar_median_df$sample <- factor(rna_tar_median_df$sample, levels = sample_levels)

p11 <- ggplot(rna_tar_median_df, aes(x = sample, y = median_counts, fill = assay)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_y_log10() +
  labs(title = "RNA and TAR depth per sample", x = "Sample",
       y = "Median counts per cell (log10)", fill = "Matrix") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p11)
save_plot(p11, "11_RNA_vs_TAR_median_counts")


# ------------------------------------------------------------
# 20. Plot 12: RNA vs TAR depth per cell
# ------------------------------------------------------------

p12 <- ggplot(cell_qc, aes(x = RNA_counts, y = TAR_counts)) +
  geom_point(alpha = 0.2, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ sample, scales = "free") +
  labs(title = "RNA vs TAR depth per cell", x = "RNA counts per cell (log10)",
       y = "TAR counts per cell (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p12)
save_plot(p12, "12_RNA_vs_TAR_cell_depth", width = 10, height = 7)


# ------------------------------------------------------------
# 21. Save updated summary (adds TAR_sparsity_percent)
# ------------------------------------------------------------

write.csv(qc_summary, file = file.path(data_dir, "TAR_QC_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 22. Done
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("QC FINISHED\n")
cat("====================================================\n")

cat("\nSummary table:\n", file.path(data_dir, "TAR_QC_summary.csv"), "\n")
cat("\nPer-cell QC:\n", file.path(data_dir, "TAR_cell_QC.csv"), "\n")
cat("\nFigures saved in:\n", qc_plot_dir, "\n")

cat("\nFiles written:\n")
print(list.files(qc_plot_dir))
