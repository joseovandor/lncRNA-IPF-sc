# ============================================================
# 05_known_ipf_lncrna_exploratory.R
#
# lncRNA-IPF-sc | Step 5: exploratory positive control
# Author: Jose A. Ovando-Ricardez
# Known lncRNAs, IPF vs control (whole-sample CPM).
#
# This is descriptive only (3 IPF vs 1 control); it is NOT a formal
# differential-expression test.
#
# Outputs: data/known_IPF_lncRNA/known_IPF_lncRNA_IPF_vs_Control_exploratory.csv,
#          data/plots/04_lncRNA/02_known_IPF_lncRNA/07_*, 08_*
# ============================================================

source("config/config.R")

library(Seurat)
library(Matrix)
library(ggplot2)

data_dir <- DATA_DIR

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")

condition <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)

known_ipf_lncrnas <- c(
  "MEG3", "FENDRR", "MALAT1", "NEAT1", "H19", "PVT1",
  "GAS5", "HOTAIR", "TUG1", "CDKN2B-AS1", "LINC00470", "UCA1"
)

results_list <- list()

for (s in samples) {

  cat("\nProcessing:", s, "\n")

  cr <- Read10X(data.dir = file.path(data_dir, s, "cellranger"))

  # Total library size of the sample
  library_size <- sum(cr)

  for (gene in known_ipf_lncrnas) {

    idx <- which(toupper(rownames(cr)) == toupper(gene))

    # Gene not present in the matrix
    if (length(idx) == 0) {
      results_list[[paste(s, gene, sep = "_")]] <- data.frame(
        sample          = s,
        condition       = unname(condition[s]),
        lncRNA          = gene,
        feature_present = FALSE,
        raw_counts      = 0,
        cells_detected  = 0,
        pct_cells       = 0,
        CPM             = 0,
        stringsAsFactors = FALSE
      )
      next
    }

    gene_counts <- Matrix::colSums(cr[idx, , drop = FALSE])

    raw_counts     <- sum(gene_counts)
    cells_detected <- sum(gene_counts > 0)
    pct_cells      <- cells_detected / ncol(cr) * 100

    # Simple library-size normalisation
    CPM <- raw_counts / library_size * 1e6

    results_list[[paste(s, gene, sep = "_")]] <- data.frame(
      sample          = s,
      condition       = unname(condition[s]),
      lncRNA          = gene,
      feature_present = TRUE,
      raw_counts      = raw_counts,
      cells_detected  = cells_detected,
      pct_cells       = pct_cells,
      CPM             = CPM,
      stringsAsFactors = FALSE
    )
  }

  rm(cr)
  gc()
}


# ------------------------------------------------------------
# Combine results
# ------------------------------------------------------------

expr_df <- do.call(rbind, results_list)
rownames(expr_df) <- NULL

expr_df$detected <- expr_df$raw_counts > 0


# ------------------------------------------------------------
# IPF vs control summary
# ------------------------------------------------------------

ipf_df     <- expr_df[expr_df$condition == "IPF", ]
control_df <- expr_df[expr_df$condition == "Control", ]

ipf_summary <- aggregate(cbind(CPM, pct_cells) ~ lncRNA, data = ipf_df, FUN = mean)
colnames(ipf_summary) <- c("lncRNA", "IPF_mean_CPM", "IPF_mean_pct_cells")

control_summary <- control_df[, c("lncRNA", "CPM", "pct_cells")]
colnames(control_summary) <- c("lncRNA", "Control_CPM", "Control_pct_cells")

comparison <- merge(ipf_summary, control_summary, by = "lncRNA")


# ------------------------------------------------------------
# Exploratory log2 fold change (small pseudocount avoids division
# by zero). Descriptive only, NOT a formal DE test.
# ------------------------------------------------------------

comparison$log2FC_IPF_vs_Control <- log2(
  (comparison$IPF_mean_CPM + 0.1) / (comparison$Control_CPM + 0.1)
)

comparison$direction <- ifelse(
  comparison$log2FC_IPF_vs_Control > 0,
  "Higher in IPF",
  ifelse(comparison$log2FC_IPF_vs_Control < 0, "Lower in IPF", "Similar")
)

comparison <- comparison[order(comparison$log2FC_IPF_vs_Control, decreasing = TRUE), ]


# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("POSITIVE CONTROL: lncRNA IPF vs CONTROL\n")
cat("====================================================\n\n")

print(comparison, row.names = FALSE)

cat("\n")
cat("====================================================\n")
cat("EXPRESSION PER SAMPLE\n")
cat("====================================================\n\n")

print(expr_df[, c("lncRNA", "sample", "condition", "raw_counts", "cells_detected",
                  "pct_cells", "CPM", "detected")], row.names = FALSE)


# ------------------------------------------------------------
# Save table
# ------------------------------------------------------------

output_dir <- file.path(data_dir, "known_IPF_lncRNA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(comparison,
          file.path(output_dir, "known_IPF_lncRNA_IPF_vs_Control_exploratory.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# Figure 1: log2FC IPF vs control
# ------------------------------------------------------------

plot_dir <- file.path(data_dir, "plots", "04_lncRNA", "02_known_IPF_lncRNA")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

p1 <- ggplot(comparison, aes(x = reorder(lncRNA, log2FC_IPF_vs_Control), y = log2FC_IPF_vs_Control)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Known lncRNAs: IPF vs control", subtitle = "Exploratory CPM-based comparison",
       x = "lncRNA", y = "log2FC (IPF mean / control)") +
  theme_bw(base_size = 12)

print(p1)

ggsave(file.path(plot_dir, "07_known_IPF_lncRNA_log2FC_IPF_vs_Control.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "07_known_IPF_lncRNA_log2FC_IPF_vs_Control.pdf"), p1, width = 8, height = 6)


# ------------------------------------------------------------
# Figure 2: CPM per sample
# ------------------------------------------------------------

expr_df$sample <- factor(expr_df$sample, levels = samples)

p2 <- ggplot(expr_df, aes(x = sample, y = CPM + 0.1, group = lncRNA)) +
  geom_point(size = 2) +
  facet_wrap(~ lncRNA, scales = "free_y") +
  scale_y_log10() +
  labs(title = "Expression of known lncRNAs per sample", x = "Sample", y = "CPM + 0.1 (log10)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p2)

ggsave(file.path(plot_dir, "08_known_IPF_lncRNA_CPM_per_sample.png"), p2, width = 12, height = 9, dpi = 300)
ggsave(file.path(plot_dir, "08_known_IPF_lncRNA_CPM_per_sample.pdf"), p2, width = 12, height = 9)

cat("\n")
cat("====================================================\n")
cat("DONE\n")
cat("====================================================\n")
