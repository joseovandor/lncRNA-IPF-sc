# ============================================================
# 04_known_ipf_lncrna_detection.R
#
# lncRNA-IPF-sc | Step 4: screen known IPF/fibrosis-related lncRNAs
# Author: Jose A. Ovando-Ricardez
#
# 1. Looks for known lncRNAs in the Cell Ranger RNA matrix.
# 2. Computes total counts, positive cells, % positive cells and
#    mean expression.
# 3. Looks for the same lncRNAs in the annotated TAR catalogue.
# 4. Integrates RNA + TAR into one summary table.
# 5. Produces PNG and PDF figures.
#
# Outputs: data/known_IPF_lncRNA/, data/plots/04_lncRNA/02_known_IPF_lncRNA/
# ============================================================

source("config/config.R")

library(Seurat)
library(Matrix)
library(ggplot2)


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

data_dir <- DATA_DIR

catalog_dir    <- file.path(data_dir, "TAR_catalog")
annotation_dir <- file.path(catalog_dir, "annotation")
plots_dir      <- file.path(data_dir, "plots", "04_lncRNA", "02_known_IPF_lncRNA")
output_dir     <- file.path(data_dir, "known_IPF_lncRNA")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nResults directory:\n")
cat(output_dir, "\n")
cat("\nFigure directory:\n")
cat(plots_dir, "\n")


# ------------------------------------------------------------
# 2. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")
sample_levels <- samples

condition <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)


# ------------------------------------------------------------
# 3. lncRNAs to evaluate (editable list; some have direct evidence
#    in IPF, others have been linked to pulmonary fibrosis or
#    fibrotic regulation)
# ------------------------------------------------------------

known_ipf_lncrnas <- c(
  "MEG3", "FENDRR", "MALAT1", "NEAT1", "H19", "PVT1",
  "GAS5", "HOTAIR", "TUG1", "CDKN2B-AS1", "LINC00470", "UCA1"
)

cat("\nlncRNAs to evaluate:\n")
print(known_ipf_lncrnas)


# ------------------------------------------------------------
# 4. Screening in the Cell Ranger RNA matrix
# ------------------------------------------------------------

rna_results_list <- list()
rna_cell_list <- list()

for (s in samples) {

  cat("\n")
  cat("====================================================\n")
  cat("RNA - Sample:", s, "\n")
  cat("====================================================\n")

  # 4.1 Read Cell Ranger matrix
  cr <- Read10X(data.dir = file.path(data_dir, s, "cellranger"))

  genes_available <- rownames(cr)

  # 4.2 Look for each lncRNA
  for (gene in known_ipf_lncrnas) {

    idx <- which(toupper(genes_available) == toupper(gene))

    # Not found
    if (length(idx) == 0) {
      rna_results_list[[paste(s, gene, sep = "_")]] <- data.frame(
        sample                         = s,
        condition                      = unname(condition[s]),
        lncRNA                         = gene,
        RNA_present                    = FALSE,
        RNA_total_counts               = 0,
        RNA_cells_detected             = 0,
        RNA_pct_cells                  = 0,
        RNA_mean_counts_all_cells      = 0,
        RNA_mean_counts_positive_cells = 0,
        stringsAsFactors = FALSE
      )
      next
    }

    # Found
    gene_matrix <- cr[idx, , drop = FALSE]
    gene_counts <- Matrix::colSums(gene_matrix)

    cells_detected <- sum(gene_counts > 0)
    pct_cells      <- cells_detected / ncol(cr) * 100
    total_counts   <- sum(gene_counts)
    mean_all       <- mean(gene_counts)

    if (cells_detected > 0) {
      mean_positive <- mean(gene_counts[gene_counts > 0])
    } else {
      mean_positive <- 0
    }

    rna_results_list[[paste(s, gene, sep = "_")]] <- data.frame(
      sample                         = s,
      condition                      = unname(condition[s]),
      lncRNA                         = gene,
      RNA_present                    = TRUE,
      RNA_total_counts               = total_counts,
      RNA_cells_detected             = cells_detected,
      RNA_pct_cells                  = pct_cells,
      RNA_mean_counts_all_cells      = mean_all,
      RNA_mean_counts_positive_cells = mean_positive,
      stringsAsFactors = FALSE
    )

    # Per-cell expression
    rna_cell_list[[paste(s, gene, sep = "_")]] <- data.frame(
      sample    = s,
      condition = unname(condition[s]),
      barcode   = colnames(cr),
      lncRNA    = gene,
      counts    = as.numeric(gene_counts),
      stringsAsFactors = FALSE
    )
  }

  rm(cr)
  gc()
}


# ------------------------------------------------------------
# 5. Combine RNA results
# ------------------------------------------------------------

rna_summary <- do.call(rbind, rna_results_list)
rownames(rna_summary) <- NULL

rna_cell <- do.call(rbind, rna_cell_list)
rownames(rna_cell) <- NULL


# ------------------------------------------------------------
# 6. Print RNA results
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RESULTS IN THE RNA MATRIX\n")
cat("====================================================\n\n")

print(rna_summary)


# ------------------------------------------------------------
# 7. Read annotated TAR catalogue
# ------------------------------------------------------------

annotated_catalog_file <- file.path(annotation_dir, "IPF_uTAR_consensus_annotated.csv")
mapping_file           <- file.path(catalog_dir, "IPF_uTAR_sample_mapping.csv")

if (!file.exists(annotated_catalog_file)) {
  stop(paste("File not found:", annotated_catalog_file))
}

if (!file.exists(mapping_file)) {
  stop(paste("File not found:", mapping_file))
}

annotated_catalog <- read.csv(annotated_catalog_file, stringsAsFactors = FALSE)
mapping           <- read.csv(mapping_file, stringsAsFactors = FALSE)


# ------------------------------------------------------------
# 8. TARs associated with the lncRNAs of interest
# ------------------------------------------------------------

tar_hits_list <- list()

for (gene in known_ipf_lncrnas) {

  hit <- grepl(paste0("(^|;)", gene, "(;|$)"),
               annotated_catalog$GENCODE_lncRNA_gene_name, ignore.case = TRUE)
  hit[is.na(hit)] <- FALSE

  tmp <- annotated_catalog[hit, , drop = FALSE]

  if (nrow(tmp) > 0) {
    tmp$query_lncRNA <- gene
    tar_hits_list[[gene]] <- tmp
  }
}


# ------------------------------------------------------------
# 9. Combine TAR hits
# ------------------------------------------------------------

if (length(tar_hits_list) > 0) {
  tar_hits <- do.call(rbind, tar_hits_list)
  rownames(tar_hits) <- NULL
} else {
  tar_hits <- data.frame()
}

cat("\n")
cat("====================================================\n")
cat("TARs ASSOCIATED WITH KNOWN lncRNAs\n")
cat("====================================================\n\n")

if (nrow(tar_hits) > 0) {
  print(tar_hits[, intersect(c("TAR_ID", "chr", "start", "end", "strand", "query_lncRNA",
                               "GENCODE_lncRNA_gene_name", "n_samples", "n_IPF", "n_Control",
                               "presence_class"),
                             colnames(tar_hits)), drop = FALSE])
} else {
  cat("No TARs associated with the selected lncRNAs were found.\n")
}


# ------------------------------------------------------------
# 10. TAR summary per lncRNA and sample
# ------------------------------------------------------------

tar_summary_list <- list()

for (gene in known_ipf_lncrnas) {

  if (nrow(tar_hits) == 0 || !(gene %in% tar_hits$query_lncRNA)) {
    for (s in samples) {
      tar_summary_list[[paste(gene, s, sep = "_")]] <- data.frame(
        sample                 = s,
        condition              = unname(condition[s]),
        lncRNA                 = gene,
        TAR_present            = FALSE,
        TAR_loci               = 0,
        TAR_total_counts       = 0,
        TAR_max_cells_detected = 0,
        TAR_max_pct_cells      = 0,
        stringsAsFactors = FALSE
      )
    }
    next
  }

  gene_tar_ids <- unique(tar_hits$TAR_ID[tar_hits$query_lncRNA == gene])

  for (s in samples) {

    gene_mapping <- mapping[mapping$TAR_ID %in% gene_tar_ids & mapping$sample == s, , drop = FALSE]

    if (nrow(gene_mapping) == 0) {
      tar_summary_list[[paste(gene, s, sep = "_")]] <- data.frame(
        sample                 = s,
        condition              = unname(condition[s]),
        lncRNA                 = gene,
        TAR_present            = FALSE,
        TAR_loci               = 0,
        TAR_total_counts       = 0,
        TAR_max_cells_detected = 0,
        TAR_max_pct_cells      = 0,
        stringsAsFactors = FALSE
      )
    } else {
      tar_summary_list[[paste(gene, s, sep = "_")]] <- data.frame(
        sample                 = s,
        condition              = unname(condition[s]),
        lncRNA                 = gene,
        TAR_present            = TRUE,
        TAR_loci               = length(unique(gene_mapping$TAR_ID)),
        TAR_total_counts       = sum(gene_mapping$total_counts, na.rm = TRUE),
        TAR_max_cells_detected = max(gene_mapping$cells_detected, na.rm = TRUE),
        TAR_max_pct_cells      = max(gene_mapping$pct_cells, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}


# ------------------------------------------------------------
# 11. Combine TAR summary
# ------------------------------------------------------------

tar_summary <- do.call(rbind, tar_summary_list)
rownames(tar_summary) <- NULL


# ------------------------------------------------------------
# 12. Integrate RNA + TAR
# ------------------------------------------------------------

combined_summary <- merge(rna_summary, tar_summary,
                          by = c("sample", "condition", "lncRNA"), all = TRUE)

combined_summary$sample <- factor(combined_summary$sample, levels = sample_levels)
combined_summary$lncRNA <- factor(combined_summary$lncRNA, levels = known_ipf_lncrnas)


# ------------------------------------------------------------
# 13. Overall detection source
# ------------------------------------------------------------

combined_summary$detected_any <- combined_summary$RNA_present | combined_summary$TAR_present

combined_summary$detection_source <- "Not detected"

combined_summary$detection_source[
  combined_summary$RNA_present & !combined_summary$TAR_present
] <- "RNA only"

combined_summary$detection_source[
  !combined_summary$RNA_present & combined_summary$TAR_present
] <- "TAR only"

combined_summary$detection_source[
  combined_summary$RNA_present & combined_summary$TAR_present
] <- "RNA + TAR"


# ------------------------------------------------------------
# 14. Save tables
# ------------------------------------------------------------

write.csv(rna_summary, file.path(output_dir, "known_IPF_lncRNA_RNA_summary.csv"), row.names = FALSE)
write.csv(tar_summary, file.path(output_dir, "known_IPF_lncRNA_TAR_summary.csv"), row.names = FALSE)
write.csv(combined_summary, file.path(output_dir, "known_IPF_lncRNA_combined_summary.csv"), row.names = FALSE)
write.csv(tar_hits, file.path(output_dir, "known_IPF_lncRNA_TAR_loci.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 15. Helper: save PNG + PDF
# ------------------------------------------------------------

save_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(filename = file.path(plots_dir, paste0(filename, ".png")),
         plot = plot_object, width = width, height = height, dpi = 300)
  ggsave(filename = file.path(plots_dir, paste0(filename, ".pdf")),
         plot = plot_object, width = width, height = height)
}


# ------------------------------------------------------------
# 16. Figure 1: % RNA+ cells
# ------------------------------------------------------------

p1 <- ggplot(combined_summary, aes(x = sample, y = lncRNA, size = RNA_pct_cells)) +
  geom_point() +
  scale_size_continuous(range = c(0, 9)) +
  labs(title = "Detection of known lncRNAs in the RNA matrix", x = "Sample", y = "lncRNA",
       size = "% positive cells") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)
save_plot(p1, "01_known_IPF_lncRNA_RNA_detection")


# ------------------------------------------------------------
# 17. Figure 2: RNA counts
# ------------------------------------------------------------

plot_df <- combined_summary
plot_df$RNA_total_counts_plot <- plot_df$RNA_total_counts + 1

p2 <- ggplot(plot_df, aes(x = sample, y = lncRNA, size = RNA_pct_cells,
                          fill = log10(RNA_total_counts_plot))) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(0, 9)) +
  labs(title = "Expression of known lncRNAs in Cell Ranger", x = "Sample", y = "lncRNA",
       size = "% positive cells", fill = "log10(counts + 1)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
save_plot(p2, "02_known_IPF_lncRNA_RNA_dotplot")


# ------------------------------------------------------------
# 18. Figure 3: RNA + TAR evidence
# ------------------------------------------------------------

p3 <- ggplot(combined_summary, aes(x = sample, y = lncRNA, shape = detection_source)) +
  geom_point(size = 4) +
  labs(title = "Detection sources of known lncRNAs", x = "Sample", y = "lncRNA", shape = "Source") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)
save_plot(p3, "03_known_IPF_lncRNA_RNA_TAR_detection")


# ------------------------------------------------------------
# 19. Figure 4: TAR counts
# ------------------------------------------------------------

tar_plot_df <- combined_summary
tar_plot_df$TAR_total_counts_plot <- tar_plot_df$TAR_total_counts + 1

p4 <- ggplot(tar_plot_df, aes(x = sample, y = lncRNA, size = TAR_max_pct_cells,
                              fill = log10(TAR_total_counts_plot))) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(0, 9)) +
  labs(title = "TAR signal associated with known lncRNAs", x = "Sample", y = "lncRNA",
       size = "Max. % TAR+ cells", fill = "log10(TAR counts + 1)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p4)
save_plot(p4, "04_known_IPF_lncRNA_TAR_dotplot")


# ------------------------------------------------------------
# 20. Figure 5: number of samples with detection per lncRNA
# ------------------------------------------------------------

presence_gene <- aggregate(detected_any ~ lncRNA, data = combined_summary, FUN = sum)
colnames(presence_gene)[2] <- "n_samples_detected"

p5 <- ggplot(presence_gene, aes(x = lncRNA, y = n_samples_detected)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_samples_detected), vjust = -0.4, size = 4) +
  ylim(0, 4.5) +
  labs(title = "Number of samples with detection of each lncRNA", x = "lncRNA",
       y = "Samples with detection") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p5)
save_plot(p5, "05_known_IPF_lncRNA_samples_detected", width = 9, height = 6)


# ------------------------------------------------------------
# 21. Figure 6: per-cell distribution (lncRNAs found in RNA only;
#     log1p to show both high and low expression)
# ------------------------------------------------------------

if (nrow(rna_cell) > 0) {

  p6 <- ggplot(rna_cell, aes(x = sample, y = log1p(counts))) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA) +
    facet_wrap(~ lncRNA, scales = "free_y") +
    labs(title = "Cellular expression of known lncRNAs", x = "Sample", y = "log1p(counts)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))

  print(p6)
  save_plot(p6, "06_known_IPF_lncRNA_cell_distribution", width = 12, height = 10)
}


# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("KNOWN lncRNA SCREENING FINISHED\n")
cat("====================================================\n\n")

cat("lncRNAs evaluated:", length(known_ipf_lncrnas), "\n\n")

cat("Presence per sample:\n\n")
print(combined_summary[, c("lncRNA", "sample", "RNA_present", "RNA_total_counts",
                           "RNA_cells_detected", "RNA_pct_cells", "TAR_present", "TAR_loci",
                           "TAR_total_counts", "TAR_max_pct_cells", "detection_source")])

cat("\n")
cat("====================================================\n")
cat("NUMBER OF SAMPLES WITH DETECTION\n")
cat("====================================================\n\n")
print(presence_gene)

cat("\n")
cat("====================================================\n")
cat("FILES WRITTEN\n")
cat("====================================================\n\n")
print(list.files(output_dir))

cat("\n")
cat("====================================================\n")
cat("FIGURES\n")
cat("====================================================\n\n")
print(list.files(plots_dir))
