# ============================================================
# 02_build_consensus_utar_catalog.R
#
# lncRNA-IPF-sc | Step 2: cross-sample consensus uTAR catalogue
# Author: Jose A. Ovando-Ricardez
#
# 1. Reads the TAR matrices of the four samples.
# 2. Extracts preliminary uTARs (features ending in "_0").
# 3. Recovers genomic coordinates, strand and expression.
# 4. Converts uTARs to GRanges.
# 5. Merges OVERLAPPING regions across samples, respecting strand.
# 6. Assigns identifiers IPFTAR000001, IPFTAR000002, ...
# 7. Builds the original uTAR -> consensus locus mapping table.
# 8. Builds a per-sample presence table.
# 9. Saves the catalogue as CSV and BED.
# 10. Produces descriptive figures of the catalogue.
#
# Notes:
# - The 500-bp merge was already applied within each sample by the
#   TAR-scRNA-seq pipeline; it is NOT applied again here.
# - Only truly overlapping regions on the same strand are merged.
# - "IPF_only_in_pilot" describes presence in this pilot only; it does
#   NOT imply biological specificity for IPF.
# - A uTAR must not be called a "novel lncRNA" at this stage.
#
# Outputs: data/TAR_catalog/IPF_uTAR_consensus_catalog.{csv,bed},
#          IPF_uTAR_sample_mapping.csv, IPF_uTAR_sample_summary.csv,
#          *_GRanges.rds, data/plots/02_TAR_catalog/
# ============================================================

source("config/config.R")

library(Seurat)
library(Matrix)
library(GenomicRanges)
library(IRanges)
library(ggplot2)


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

data_dir <- DATA_DIR

plots_dir        <- file.path(data_dir, "plots")
catalog_plot_dir <- file.path(plots_dir, "02_TAR_catalog")
catalog_dir      <- file.path(data_dir, "TAR_catalog")

dir.create(catalog_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(catalog_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nCatalogue directory:\n")
cat(catalog_dir, "\n")
cat("\nFigure directory:\n")
cat(catalog_plot_dir, "\n")


# ------------------------------------------------------------
# 2. Samples and conditions
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")
sample_levels <- samples

condition <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)

ipf_samples     <- c("VUILD53", "VUILD63", "VUILD64")
control_samples <- c("VUHD67")


# ------------------------------------------------------------
# 3. Result containers
# ------------------------------------------------------------

utar_list <- list()
sample_summary <- data.frame()


# ------------------------------------------------------------
# 4. Extract uTARs from each sample
# ------------------------------------------------------------

for (s in samples) {

  cat("\n")
  cat("====================================================\n")
  cat("Processing sample:", s, "\n")
  cat("Condition:", unname(condition[s]), "\n")
  cat("====================================================\n")

  # 4.1 Read TAR matrix
  tar <- Read10X(data.dir = file.path(data_dir, s, "TAR", "TAR_feature_bc_matrix"))

  # 4.2 Select preliminary uTARs ("_0" = no associated annotation in the
  #     feature names produced by the pipeline)
  utar_names <- rownames(tar)[grepl("_0$", rownames(tar))]

  cat("\nPreliminary uTARs:", length(utar_names), "\n")

  # 4.3 Parse coordinates with a regular expression.
  #     Example: chrX_94354499_94356249_-_66_0
  #     Parsed from the end to tolerate contig names containing "_".
  pattern <- "^(.+)_([0-9]+)_([0-9]+)_([+-])_([0-9]+)_0$"

  matched <- regexec(pattern, utar_names)
  parsed  <- regmatches(utar_names, matched)
  valid   <- lengths(parsed) == 6

  if (any(!valid)) {
    cat("WARNING:", sum(!valid), "uTARs could not be parsed and will be excluded.\n")
    writeLines(utar_names[!valid],
               con = file.path(catalog_dir, paste0(s, "_uTAR_unparsed.txt")))
  }

  parsed      <- parsed[valid]
  valid_names <- utar_names[valid]

  utar_df <- data.frame(
    feature         = valid_names,
    chr             = vapply(parsed, `[`, character(1), 2),
    start           = as.numeric(vapply(parsed, `[`, character(1), 3)),
    end             = as.numeric(vapply(parsed, `[`, character(1), 4)),
    strand          = vapply(parsed, `[`, character(1), 5),
    feature_support = as.numeric(vapply(parsed, `[`, character(1), 6)),
    sample          = s,
    condition       = unname(condition[s]),
    stringsAsFactors = FALSE
  )

  # 4.4 Expression of each uTAR
  utar_matrix <- tar[valid_names, , drop = FALSE]

  utar_df$total_counts   <- as.numeric(Matrix::rowSums(utar_matrix))
  utar_df$cells_detected <- as.numeric(Matrix::rowSums(utar_matrix > 0))
  utar_df$pct_cells      <- utar_df$cells_detected / ncol(utar_matrix) * 100

  # 4.5 Length of each original uTAR
  utar_df$length <- utar_df$end - utar_df$start + 1

  # 4.6 Convert to GRanges
  gr <- GRanges(
    seqnames = utar_df$chr,
    ranges   = IRanges(start = utar_df$start, end = utar_df$end),
    strand   = utar_df$strand
  )

  mcols(gr)$sample          <- utar_df$sample
  mcols(gr)$condition       <- utar_df$condition
  mcols(gr)$feature         <- utar_df$feature
  mcols(gr)$feature_support <- utar_df$feature_support
  mcols(gr)$total_counts    <- utar_df$total_counts
  mcols(gr)$cells_detected  <- utar_df$cells_detected
  mcols(gr)$pct_cells       <- utar_df$pct_cells
  mcols(gr)$original_length <- utar_df$length

  # 4.7 Store GRanges
  utar_list[[s]] <- gr

  # 4.8 Per-sample summary
  sample_summary <- rbind(
    sample_summary,
    data.frame(
      sample                = s,
      condition             = unname(condition[s]),
      cells                 = ncol(tar),
      total_TAR             = nrow(tar),
      original_uTAR         = length(valid_names),
      median_uTAR_counts    = median(utar_df$total_counts),
      median_cells_detected = median(utar_df$cells_detected),
      median_pct_cells      = median(utar_df$pct_cells),
      median_uTAR_length    = median(utar_df$length),
      stringsAsFactors = FALSE
    )
  )

  rm(tar, utar_matrix, utar_df, gr, parsed, matched)
  gc()
}


# ------------------------------------------------------------
# 5. Combine all uTARs
# ------------------------------------------------------------

all_utar <- do.call(c, unname(utar_list))

cat("\n")
cat("====================================================\n")
cat("uTARs BEFORE BUILDING THE CONSENSUS CATALOGUE\n")
cat("====================================================\n")
cat("\nTotal original uTARs:", length(all_utar), "\n")


# ------------------------------------------------------------
# 6. Consensus catalogue
#    reduce() merges only overlapping regions:
#    - min.gapwidth = 1: regions must truly overlap to be merged.
#    - ignore.strand = FALSE: + and - are kept separate.
# ------------------------------------------------------------

consensus_utar <- reduce(all_utar, ignore.strand = FALSE, min.gapwidth = 1)

cat("Consensus uTARs:", length(consensus_utar), "\n")


# ------------------------------------------------------------
# 7. Unique identifiers
# ------------------------------------------------------------

mcols(consensus_utar)$TAR_ID <- sprintf("IPFTAR%06d", seq_along(consensus_utar))


# ------------------------------------------------------------
# 8. Original uTAR -> consensus mapping
# ------------------------------------------------------------

hits <- findOverlaps(consensus_utar, all_utar, ignore.strand = FALSE)

mapping <- data.frame(
  TAR_ID           = mcols(consensus_utar)$TAR_ID[queryHits(hits)],
  consensus_chr    = as.character(seqnames(consensus_utar[queryHits(hits)])),
  consensus_start  = start(consensus_utar[queryHits(hits)]),
  consensus_end    = end(consensus_utar[queryHits(hits)]),
  consensus_strand = as.character(strand(consensus_utar[queryHits(hits)])),
  sample           = mcols(all_utar)$sample[subjectHits(hits)],
  condition        = mcols(all_utar)$condition[subjectHits(hits)],
  original_feature = mcols(all_utar)$feature[subjectHits(hits)],
  feature_support  = mcols(all_utar)$feature_support[subjectHits(hits)],
  total_counts     = mcols(all_utar)$total_counts[subjectHits(hits)],
  cells_detected   = mcols(all_utar)$cells_detected[subjectHits(hits)],
  pct_cells        = mcols(all_utar)$pct_cells[subjectHits(hits)],
  original_length  = mcols(all_utar)$original_length[subjectHits(hits)],
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. Validate mapping
# ------------------------------------------------------------

cat("\nOriginal uTARs mapped:", nrow(mapping), "\n")
cat("Original uTARs expected:", length(all_utar), "\n")

if (nrow(mapping) != length(all_utar)) {
  warning("The number of mapping rows does not match the number of original uTARs.")
}


# ------------------------------------------------------------
# 10. Per-sample presence table
# ------------------------------------------------------------

presence_table <- table(mapping$TAR_ID, mapping$sample)
presence_table <- (presence_table > 0) * 1

presence_df <- as.data.frame.matrix(presence_table)
presence_df$TAR_ID <- rownames(presence_df)
rownames(presence_df) <- NULL

# Make sure every sample exists as a column
for (s in samples) {
  if (!(s %in% colnames(presence_df))) {
    presence_df[[s]] <- 0
  }
}


# ------------------------------------------------------------
# 11. Master catalogue table
# ------------------------------------------------------------

consensus_df <- data.frame(
  TAR_ID = mcols(consensus_utar)$TAR_ID,
  chr    = as.character(seqnames(consensus_utar)),
  start  = start(consensus_utar),
  end    = end(consensus_utar),
  strand = as.character(strand(consensus_utar)),
  length = width(consensus_utar),
  stringsAsFactors = FALSE
)

consensus_df <- merge(consensus_df, presence_df, by = "TAR_ID", all.x = TRUE)


# ------------------------------------------------------------
# 12-14. Number of samples / IPF samples / controls with detection
# ------------------------------------------------------------

consensus_df$n_samples <- rowSums(consensus_df[, samples, drop = FALSE])
consensus_df$n_IPF     <- rowSums(consensus_df[, ipf_samples, drop = FALSE])
consensus_df$n_Control <- rowSums(consensus_df[, control_samples, drop = FALSE])


# ------------------------------------------------------------
# 15. Descriptive presence class
#     "IPF_only_in_pilot" does NOT mean "IPF-specific": it only means
#     detected in >=1 IPF sample and not detected in VUHD67.
# ------------------------------------------------------------

consensus_df$presence_class <- "Other"

consensus_df$presence_class[
  consensus_df$n_IPF > 0 & consensus_df$n_Control == 0
] <- "IPF_only_in_pilot"

consensus_df$presence_class[
  consensus_df$n_IPF == 0 & consensus_df$n_Control > 0
] <- "Control_only_in_pilot"

consensus_df$presence_class[
  consensus_df$n_IPF > 0 & consensus_df$n_Control > 0
] <- "Shared_IPF_Control"


# ------------------------------------------------------------
# 16. Exact presence pattern
# ------------------------------------------------------------

consensus_df$presence_pattern <- apply(
  consensus_df[, samples, drop = FALSE],
  1,
  function(x) {
    present <- samples[as.numeric(x) == 1]
    if (length(present) == 0) {
      return("None")
    }
    paste(present, collapse = " + ")
  }
)


# ------------------------------------------------------------
# 17. Length >= 200 bp
# ------------------------------------------------------------

consensus_df$length_ge_200 <- consensus_df$length >= 200


# ------------------------------------------------------------
# 18. Expression summary per consensus locus (median over the
#     original uTARs mapped to it)
# ------------------------------------------------------------

mapping_summary <- aggregate(
  cbind(total_counts, cells_detected, pct_cells) ~ TAR_ID,
  data = mapping,
  FUN = median
)

colnames(mapping_summary)[colnames(mapping_summary) == "total_counts"]   <- "median_total_counts_original"
colnames(mapping_summary)[colnames(mapping_summary) == "cells_detected"] <- "median_cells_detected_original"
colnames(mapping_summary)[colnames(mapping_summary) == "pct_cells"]      <- "median_pct_cells_original"

consensus_df <- merge(consensus_df, mapping_summary, by = "TAR_ID", all.x = TRUE)


# ------------------------------------------------------------
# 19. Sort catalogue
# ------------------------------------------------------------

consensus_df <- consensus_df[order(consensus_df$chr, consensus_df$start, consensus_df$end), ]
rownames(consensus_df) <- NULL


# ------------------------------------------------------------
# 20. Save tables
# ------------------------------------------------------------

write.csv(consensus_df, file = file.path(catalog_dir, "IPF_uTAR_consensus_catalog.csv"), row.names = FALSE)
write.csv(mapping, file = file.path(catalog_dir, "IPF_uTAR_sample_mapping.csv"), row.names = FALSE)
write.csv(sample_summary, file = file.path(catalog_dir, "IPF_uTAR_sample_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 21. Save BED (0-based start: BED start = R start - 1, BED end = R end)
# ------------------------------------------------------------

bed_df <- data.frame(
  chr    = consensus_df$chr,
  start  = consensus_df$start - 1,
  end    = consensus_df$end,
  name   = consensus_df$TAR_ID,
  score  = 0,
  strand = consensus_df$strand
)

write.table(bed_df, file = file.path(catalog_dir, "IPF_uTAR_consensus_catalog.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 22. Save R objects
# ------------------------------------------------------------

saveRDS(consensus_utar, file = file.path(catalog_dir, "IPF_uTAR_consensus_GRanges.rds"))
saveRDS(utar_list, file = file.path(catalog_dir, "IPF_uTAR_per_sample_GRanges.rds"))


# ------------------------------------------------------------
# 23. Helper: save PNG + PDF
# ------------------------------------------------------------

save_catalog_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(filename = file.path(catalog_plot_dir, paste0(filename, ".png")),
         plot = plot_object, width = width, height = height, dpi = 300)
  ggsave(filename = file.path(catalog_plot_dir, paste0(filename, ".pdf")),
         plot = plot_object, width = width, height = height)
}


# ------------------------------------------------------------
# 24. Factors
# ------------------------------------------------------------

sample_summary$sample <- factor(sample_summary$sample, levels = sample_levels)
mapping$sample        <- factor(mapping$sample, levels = sample_levels)


# ------------------------------------------------------------
# 25. Figure 1: original uTARs per sample
# ------------------------------------------------------------

p1 <- ggplot(sample_summary, aes(x = sample, y = original_uTAR)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(original_uTAR, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Preliminary uTARs detected per sample", x = "Sample", y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)
save_catalog_plot(p1, "01_original_uTAR_per_sample")


# ------------------------------------------------------------
# 26. Figure 2: redundancy reduction
# ------------------------------------------------------------

reduction_df <- data.frame(
  category = factor(c("Original uTARs", "Consensus uTARs"),
                    levels = c("Original uTARs", "Consensus uTARs")),
  count = c(length(all_utar), nrow(consensus_df))
)

p2 <- ggplot(reduction_df, aes(x = category, y = count)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Redundancy reduction when building the catalogue", x = NULL, y = "Number of regions") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
save_catalog_plot(p2, "02_original_vs_consensus_uTAR")


# ------------------------------------------------------------
# 27. Figure 3: number of samples in which each locus is detected
# ------------------------------------------------------------

recurrence_df <- as.data.frame(table(consensus_df$n_samples))
colnames(recurrence_df) <- c("n_samples", "count")
recurrence_df$n_samples <- as.numeric(as.character(recurrence_df$n_samples))

p3 <- ggplot(recurrence_df, aes(x = factor(n_samples), y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Recurrence of consensus uTARs across samples", x = "Number of samples",
       y = "Number of consensus uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)
save_catalog_plot(p3, "03_consensus_uTAR_recurrence")


# ------------------------------------------------------------
# 28. Figure 4: IPF / control presence class
# ------------------------------------------------------------

presence_class_df <- as.data.frame(table(consensus_df$presence_class))
colnames(presence_class_df) <- c("presence_class", "count")

p4 <- ggplot(presence_class_df, aes(x = presence_class, y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Descriptive distribution of uTARs by condition", x = "Overall presence pattern",
       y = "Number of consensus uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p4)
save_catalog_plot(p4, "04_consensus_uTAR_presence_class")


# ------------------------------------------------------------
# 29. Figure 5: length distribution
# ------------------------------------------------------------

p5 <- ggplot(consensus_df, aes(x = length)) +
  geom_histogram(bins = 80) +
  scale_x_log10() +
  labs(title = "Length distribution of consensus uTARs", x = "uTAR length (bp, log10)",
       y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p5)
save_catalog_plot(p5, "05_consensus_uTAR_length_distribution")


# ------------------------------------------------------------
# 30. Figure 6: >=200 bp vs <200 bp
# ------------------------------------------------------------

length_class_df <- as.data.frame(table(ifelse(consensus_df$length_ge_200, ">=200 bp", "<200 bp")))
colnames(length_class_df) <- c("length_class", "count")

p6 <- ggplot(length_class_df, aes(x = length_class, y = count)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Length of consensus uTARs", x = "Length class", y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p6)
save_catalog_plot(p6, "06_consensus_uTAR_length_200bp")


# ------------------------------------------------------------
# 31. Figure 7: distribution by chromosome
# ------------------------------------------------------------

chr_df <- as.data.frame(table(consensus_df$chr))
colnames(chr_df) <- c("chr", "count")
chr_df <- chr_df[order(chr_df$count, decreasing = TRUE), ]
chr_df$chr <- factor(chr_df$chr, levels = rev(chr_df$chr))

p7 <- ggplot(chr_df, aes(x = chr, y = count)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(title = "Genomic distribution of consensus uTARs", x = "Chromosome / contig",
       y = "Number of uTARs") +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5))

print(p7)
save_catalog_plot(p7, "07_consensus_uTAR_by_chromosome",
                  width = 9, height = max(7, 0.22 * nrow(chr_df)))


# ------------------------------------------------------------
# 32. Figure 8: exact presence patterns
# ------------------------------------------------------------

pattern_df <- as.data.frame(table(consensus_df$presence_pattern))
colnames(pattern_df) <- c("pattern", "count")
pattern_df <- pattern_df[order(pattern_df$count, decreasing = TRUE), ]
pattern_df$pattern <- factor(pattern_df$pattern, levels = rev(pattern_df$pattern))

p8 <- ggplot(pattern_df, aes(x = pattern, y = count)) +
  geom_col(width = 0.75) +
  coord_flip() +
  labs(title = "uTAR presence combinations across samples", x = "Samples with detection",
       y = "Number of consensus uTARs") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5))

print(p8)
save_catalog_plot(p8, "08_consensus_uTAR_presence_patterns", width = 10, height = 7)


# ------------------------------------------------------------
# 33. Figure 9: counts of original uTARs per sample
# ------------------------------------------------------------

p9 <- ggplot(mapping, aes(x = sample, y = total_counts)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Distribution of uTAR counts per sample", x = "Sample",
       y = "Counts per uTAR (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p9)
save_catalog_plot(p9, "09_original_uTAR_counts_distribution")


# ------------------------------------------------------------
# 34. Figure 10: cells detecting each uTAR
# ------------------------------------------------------------

p10 <- ggplot(mapping, aes(x = sample, y = cells_detected)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "Number of cells detecting each uTAR", x = "Sample",
       y = "Cells expressing the uTAR (log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p10)
save_catalog_plot(p10, "10_original_uTAR_cells_detected")


# ------------------------------------------------------------
# 35. Figure 11: percentage of cells detecting each uTAR
# ------------------------------------------------------------

p11 <- ggplot(mapping, aes(x = sample, y = pct_cells)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(title = "Fraction of cells expressing each uTAR", x = "Sample", y = "Positive cells (%)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p11)
save_catalog_plot(p11, "11_original_uTAR_percent_cells")


# ------------------------------------------------------------
# 36. Figure 12: recurrence across IPF samples
# ------------------------------------------------------------

ipf_recurrence_df <- as.data.frame(table(consensus_df$n_IPF))
colnames(ipf_recurrence_df) <- c("n_IPF", "count")
ipf_recurrence_df$n_IPF <- as.numeric(as.character(ipf_recurrence_df$n_IPF))

p12 <- ggplot(ipf_recurrence_df, aes(x = factor(n_IPF), y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "uTAR recurrence across the three IPF samples", x = "Number of IPF samples",
       y = "Number of consensus uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p12)
save_catalog_plot(p12, "12_consensus_uTAR_IPF_recurrence")


# ------------------------------------------------------------
# 37. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CONSENSUS CATALOGUE FINISHED\n")
cat("====================================================\n")

cat("\nOriginal uTARs:", length(all_utar), "\n")
cat("Consensus uTARs:", nrow(consensus_df), "\n")
cat("\nReduction:", round((1 - nrow(consensus_df) / length(all_utar)) * 100, 2), "%\n")
cat("\nConsensus uTARs >=200 bp:", sum(consensus_df$length_ge_200), "of", nrow(consensus_df),
    "(", round(mean(consensus_df$length_ge_200) * 100, 2), "% )\n")

cat("\nDistribution by number of samples:\n")
print(table(consensus_df$n_samples))

cat("\nDescriptive distribution by condition:\n")
print(table(consensus_df$presence_class))

cat("\nRecurrence in IPF samples:\n")
print(table(consensus_df$n_IPF))

cat("\nLength summary:\n")
print(summary(consensus_df$length))

cat("\nCatalogue files:\n")
print(list.files(catalog_dir))

cat("\nFigures:\n")
print(list.files(catalog_plot_dir))

cat("\n====================================================\n")
cat("NEXT STEP\n")
cat("====================================================\n")
cat(
  "\nAnnotate the consensus uTARs against\n",
  "GENCODE and lncRNA-specific databases to separate:\n",
  "1) known lncRNAs\n",
  "2) uTARs with external evidence\n",
  "3) potentially novel candidates\n",
  "without calling a uTAR a novel lncRNA only because it is unannotated.\n"
)
