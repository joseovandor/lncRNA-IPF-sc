# ============================================================
# 11_priorityA_structural_refinement.R
#
# lncRNA-IPF-sc | Step 11: structural refinement of Priority A candidates
# Author: Jose A. Ovando-Ricardez
#
# For each of the 48 Priority A candidates, finds the best-overlapping
# per-sample TAR (TAR_reads.bed.gz) in VUILD53, VUILD63, VUILD64 and VUHD67,
# computes overlap bp, candidate/local-TAR fraction covered, Jaccard and
# boundary differences, summarises reproducibility across the 3 IPF samples,
# checks structural signal in the control and accounts for cluster complexity.
# Candidates are classified as stable_single_locus, needs_BAM_refinement,
# fragmented_complex_locus or poorly_reproducible.
#
# Note: this classification is STRUCTURAL only. It does NOT imply differential
# expression, IPF specificity, a validated functional lncRNA or a confirmed
# mature transcript. TAR_reads.bed.gz carries no splice-junction information.
#
# Inputs:  data/TAR_catalog/genomic_context/Unannotated_3of3_IPF_genomic_context.csv
#          data/<sample>/TAR/TAR_reads.bed.gz
# Outputs: data/TAR_catalog/priorityA_systematic_refinement/ (incl. IGV/)
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(GenomicRanges)
  library(IRanges)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

INPUT_CONTEXT <- file.path(BASE_DIR, "data/TAR_catalog/genomic_context/Unannotated_3of3_IPF_genomic_context.csv")
OUTPUT_DIR    <- file.path(BASE_DIR, "data/TAR_catalog/priorityA_systematic_refinement")
IGV_DIR       <- file.path(OUTPUT_DIR, "IGV")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(IGV_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")

conditions <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)

IPF_samples <- c("VUILD53", "VUILD63", "VUILD64")

CONTROL_sample <- "VUHD67"

tar_bed_paths <- setNames(file.path(BASE_DIR, "data", samples, "TAR", "TAR_reads.bed.gz"), samples)


# ------------------------------------------------------------
# 3. Check inputs
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INPUT CHECK\n")
cat("====================================================\n\n")

if (!file.exists(INPUT_CONTEXT)) {
  stop(paste("File not found:", INPUT_CONTEXT))
}

missing_tar_bed <- names(tar_bed_paths)[!file.exists(tar_bed_paths)]

if (length(missing_tar_bed) > 0) {
  stop(paste("Missing TAR_reads.bed.gz for:", paste(missing_tar_bed, collapse = ", ")))
}

input_check <- data.frame(
  sample    = samples,
  condition = unname(conditions[samples]),
  TAR_bed   = unname(tar_bed_paths[samples]),
  exists    = file.exists(tar_bed_paths[samples]),
  stringsAsFactors = FALSE
)

print(input_check)


# ------------------------------------------------------------
# 4. Load genomic-context catalogue
# ------------------------------------------------------------

result <- read.csv(INPUT_CONTEXT, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c(
  "TAR_ID", "chr", "start", "end", "strand", "length",
  "priority_group",
  "cluster_id", "cluster_n_TAR", "cluster_span_bp",
  "genomic_context",
  "overlapping_gene_names", "overlapping_gene_ids",
  "nearest_gene", "nearest_gene_id", "nearest_gene_type", "nearest_gene_distance",
  "nearest_gene_chr", "nearest_gene_start", "nearest_gene_end", "nearest_gene_strand",
  "genes_within_10kb", "genes_within_50kb", "genes_within_100kb",
  "IPF_min_counts", "IPF_mean_counts", "IPF_CV_counts", "IPF_min_pct_cells"
)

missing_columns <- setdiff(required_columns, colnames(result))

if (length(missing_columns) > 0) {
  cat("\nAvailable columns:\n")
  print(colnames(result))
  stop(paste("Missing columns:", paste(missing_columns, collapse = ", ")))
}

cat("\nCatalogue loaded:", nrow(result), "candidates\n")


# ------------------------------------------------------------
# 5. Extract Priority A
# ------------------------------------------------------------

priority_A <- result %>%
  dplyr::filter(.data$priority_group == "Priority_A") %>%
  dplyr::arrange(.data$chr, .data$start)

cat("\n")
cat("====================================================\n")
cat("PRIORITY A\n")
cat("====================================================\n\n")

cat("Priority A:", nrow(priority_A), "\n")

if (nrow(priority_A) != 48) {
  warning(paste("Expected 48 Priority A candidates but found", nrow(priority_A)))
}

if (nrow(priority_A) == 0) {
  stop("No Priority A candidates found.")
}


# ------------------------------------------------------------
# 6. Cluster information
# ------------------------------------------------------------

cluster_summary <- result %>%
  dplyr::group_by(.data$cluster_id) %>%
  dplyr::summarise(
    n_candidate_TAR_in_cluster = dplyr::n(),
    cluster_start = min(.data$start, na.rm = TRUE),
    cluster_end   = max(.data$end, na.rm = TRUE),
    cluster_span  = .data$cluster_end - .data$cluster_start + 1,
    .groups = "drop"
  )

priorityA_cluster_count <- priority_A %>%
  dplyr::count(.data$cluster_id, name = "n_PriorityA_in_cluster")

priority_A <- priority_A %>%
  dplyr::left_join(cluster_summary, by = "cluster_id") %>%
  dplyr::left_join(priorityA_cluster_count, by = "cluster_id")


# ------------------------------------------------------------
# 7. Priority A GRanges
# ------------------------------------------------------------

candidate_gr <- GRanges(
  seqnames = priority_A$chr,
  ranges   = IRanges(start = priority_A$start, end = priority_A$end),
  strand   = priority_A$strand
)

names(candidate_gr) <- priority_A$TAR_ID

mcols(candidate_gr)$cluster_id <- priority_A$cluster_id


# ------------------------------------------------------------
# 8. Read per-sample TAR BED files
# ------------------------------------------------------------

read_tar_bed <- function(path, sample_name, condition_name) {

  x <- readr::read_tsv(
    path,
    col_names = c("chr", "bed_start", "bed_end", "col4", "col5", "strand", "bed_support"),
    col_types = readr::cols(
      chr         = readr::col_character(),
      bed_start   = readr::col_double(),
      bed_end     = readr::col_double(),
      col4        = readr::col_character(),
      col5        = readr::col_character(),
      strand      = readr::col_character(),
      bed_support = readr::col_double()
    ),
    progress = FALSE
  )

  # BED is 0-based half-open; convert to 1-based closed coordinates
  x <- x %>%
    dplyr::mutate(
      sample          = sample_name,
      condition       = condition_name,
      start           = .data$bed_start + 1,
      end             = .data$bed_end,
      interval_length = .data$end - .data$start + 1
    )

  x
}

cat("\n")
cat("====================================================\n")
cat("READING TAR BED\n")
cat("====================================================\n\n")

sample_tar_list <- list()

for (sample_name in samples) {

  cat("Reading:", sample_name, "\n")

  tmp <- read_tar_bed(
    path           = tar_bed_paths[sample_name],
    sample_name    = sample_name,
    condition_name = unname(conditions[sample_name])
  )

  cat("  TAR:", nrow(tmp), "\n")

  sample_tar_list[[sample_name]] <- tmp
}


# ------------------------------------------------------------
# 9. TAR BED summary
# ------------------------------------------------------------

sample_summary_list <- list()

for (sample_name in samples) {

  tmp <- sample_tar_list[[sample_name]]

  sample_summary_list[[sample_name]] <- data.frame(
    sample         = sample_name,
    condition      = unname(conditions[sample_name]),
    n_TAR          = nrow(tmp),
    total_support  = sum(tmp$bed_support, na.rm = TRUE),
    median_support = median(tmp$bed_support, na.rm = TRUE),
    median_length  = median(tmp$interval_length, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

sample_summary <- dplyr::bind_rows(sample_summary_list)

cat("\n")
cat("====================================================\n")
cat("TAR BED SUMMARY\n")
cat("====================================================\n\n")

print(sample_summary)

write.csv(sample_summary, file.path(OUTPUT_DIR, "TAR_BED_sample_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 10. Per-sample overlaps
# ------------------------------------------------------------

overlap_list <- list()

overlap_counter <- 1

cat("\n")
cat("====================================================\n")
cat("COMPUTING PRIORITY A OVERLAPS\n")
cat("====================================================\n\n")

for (sample_name in samples) {

  cat("Processing:", sample_name, "\n")

  sample_df <- sample_tar_list[[sample_name]]

  sample_gr <- GRanges(
    seqnames = sample_df$chr,
    ranges   = IRanges(start = sample_df$start, end = sample_df$end),
    strand   = sample_df$strand
  )

  mcols(sample_gr)$bed_support <- sample_df$bed_support

  hits <- findOverlaps(candidate_gr, sample_gr, ignore.strand = FALSE)

  if (length(hits) == 0) {
    next
  }

  for (i in seq_along(hits)) {

    q <- queryHits(hits)[i]
    s <- subjectHits(hits)[i]

    candidate  <- candidate_gr[q]
    sample_tar <- sample_gr[s]

    overlap_start <- max(start(candidate), start(sample_tar))
    overlap_end   <- min(end(candidate), end(sample_tar))
    overlap_bp    <- max(0, overlap_end - overlap_start + 1)

    candidate_width <- as.numeric(width(candidate))
    sample_width    <- as.numeric(width(sample_tar))

    union_bp <- candidate_width + sample_width - overlap_bp

    overlap_list[[overlap_counter]] <- data.frame(
      TAR_ID                      = names(candidate),
      cluster_id                  = as.character(mcols(candidate)$cluster_id),
      sample                      = sample_name,
      condition                   = unname(conditions[sample_name]),
      chr                         = as.character(seqnames(candidate)),
      strand                      = as.character(strand(candidate)),
      candidate_start             = as.numeric(start(candidate)),
      candidate_end               = as.numeric(end(candidate)),
      candidate_width             = candidate_width,
      sample_TAR_start            = as.numeric(start(sample_tar)),
      sample_TAR_end              = as.numeric(end(sample_tar)),
      sample_TAR_width            = sample_width,
      bed_support                 = as.numeric(mcols(sample_tar)$bed_support),
      overlap_bp                  = overlap_bp,
      candidate_fraction_covered  = overlap_bp / candidate_width,
      sample_TAR_fraction_covered = overlap_bp / sample_width,
      jaccard                     = overlap_bp / union_bp,
      start_difference            = as.numeric(start(sample_tar) - start(candidate)),
      end_difference              = as.numeric(end(sample_tar) - end(candidate)),
      abs_start_difference        = abs(as.numeric(start(sample_tar) - start(candidate))),
      abs_end_difference          = abs(as.numeric(end(sample_tar) - end(candidate))),
      stringsAsFactors = FALSE
    )

    overlap_counter <- overlap_counter + 1
  }
}

if (length(overlap_list) == 0) {
  stop("No Priority A / TAR_reads overlaps found.")
}

candidate_overlap <- dplyr::bind_rows(overlap_list)

write.csv(candidate_overlap, file.path(OUTPUT_DIR, "PriorityA_all_overlaps.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 11. Best match per candidate / sample
# ------------------------------------------------------------

best_overlap <- candidate_overlap %>%
  dplyr::group_by(.data$TAR_ID, .data$sample) %>%
  dplyr::arrange(dplyr::desc(.data$jaccard), dplyr::desc(.data$candidate_fraction_covered),
                 dplyr::desc(.data$bed_support), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

write.csv(best_overlap, file.path(OUTPUT_DIR, "PriorityA_best_overlap_observed.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 12. Complete 48 x 4 candidate-by-sample matrix
# ------------------------------------------------------------

candidate_sample_grid <- tidyr::expand_grid(TAR_ID = priority_A$TAR_ID, sample = samples) %>%
  dplyr::mutate(condition = unname(conditions[.data$sample]))

candidate_sample_complete <- candidate_sample_grid %>%
  dplyr::left_join(best_overlap %>% dplyr::select(-condition), by = c("TAR_ID", "sample")) %>%
  dplyr::mutate(
    detected                    = ifelse(is.na(.data$jaccard), 0L, 1L),
    candidate_fraction_covered  = tidyr::replace_na(.data$candidate_fraction_covered, 0),
    sample_TAR_fraction_covered = tidyr::replace_na(.data$sample_TAR_fraction_covered, 0),
    jaccard                     = tidyr::replace_na(.data$jaccard, 0),
    overlap_bp                  = tidyr::replace_na(.data$overlap_bp, 0),
    bed_support                 = tidyr::replace_na(.data$bed_support, 0)
  )

write.csv(candidate_sample_complete, file.path(OUTPUT_DIR, "PriorityA_best_match_48x4.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 13. IPF summary
# ------------------------------------------------------------

IPF_summary <- candidate_sample_complete %>%
  dplyr::filter(.data$condition == "IPF") %>%
  dplyr::group_by(.data$TAR_ID) %>%
  dplyr::summarise(
    n_IPF_detected                = sum(.data$detected),
    min_IPF_candidate_fraction    = min(.data$candidate_fraction_covered),
    mean_IPF_candidate_fraction   = mean(.data$candidate_fraction_covered),
    median_IPF_candidate_fraction = median(.data$candidate_fraction_covered),
    min_IPF_jaccard               = min(.data$jaccard),
    mean_IPF_jaccard              = mean(.data$jaccard),
    median_IPF_jaccard            = median(.data$jaccard),
    n_IPF_jaccard_ge_0.50         = sum(.data$jaccard >= 0.50),
    n_IPF_jaccard_ge_0.80         = sum(.data$jaccard >= 0.80),
    start_range_IPF = ifelse(
      sum(.data$detected) >= 2,
      max(.data$sample_TAR_start[.data$detected == 1], na.rm = TRUE) -
        min(.data$sample_TAR_start[.data$detected == 1], na.rm = TRUE),
      NA_real_
    ),
    end_range_IPF = ifelse(
      sum(.data$detected) >= 2,
      max(.data$sample_TAR_end[.data$detected == 1], na.rm = TRUE) -
        min(.data$sample_TAR_end[.data$detected == 1], na.rm = TRUE),
      NA_real_
    ),
    mean_abs_start_difference_IPF = ifelse(
      any(.data$detected == 1),
      mean(.data$abs_start_difference[.data$detected == 1], na.rm = TRUE),
      NA_real_
    ),
    mean_abs_end_difference_IPF = ifelse(
      any(.data$detected == 1),
      mean(.data$abs_end_difference[.data$detected == 1], na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 14. Control summary
# ------------------------------------------------------------

control_summary <- candidate_sample_complete %>%
  dplyr::filter(.data$condition == "Control") %>%
  dplyr::transmute(
    TAR_ID                     = .data$TAR_ID,
    control_detected           = .data$detected,
    control_candidate_fraction = .data$candidate_fraction_covered,
    control_jaccard            = .data$jaccard,
    control_overlap_bp         = .data$overlap_bp,
    control_bed_support        = .data$bed_support,
    control_TAR_start          = .data$sample_TAR_start,
    control_TAR_end            = .data$sample_TAR_end
  )


# ------------------------------------------------------------
# 15. Assemble refinement table
# ------------------------------------------------------------

candidate_refinement <- priority_A %>%
  dplyr::select(
    TAR_ID, chr, start, end, strand, length,
    cluster_id, cluster_n_TAR, n_candidate_TAR_in_cluster, n_PriorityA_in_cluster,
    cluster_start, cluster_end, cluster_span,
    genomic_context,
    overlapping_gene_names, overlapping_gene_ids,
    nearest_gene, nearest_gene_id, nearest_gene_type, nearest_gene_distance,
    nearest_gene_chr, nearest_gene_start, nearest_gene_end, nearest_gene_strand,
    genes_within_10kb, genes_within_50kb, genes_within_100kb,
    IPF_min_counts, IPF_mean_counts, IPF_CV_counts, IPF_min_pct_cells
  ) %>%
  dplyr::left_join(IPF_summary, by = "TAR_ID") %>%
  dplyr::left_join(control_summary, by = "TAR_ID")


# ------------------------------------------------------------
# 16. Structural reproducibility across IPF samples
# ------------------------------------------------------------

candidate_refinement <- candidate_refinement %>%
  dplyr::mutate(
    structural_reproducibility = dplyr::case_when(
      .data$n_IPF_detected == 3 & .data$min_IPF_candidate_fraction >= 0.50 &
        .data$median_IPF_jaccard >= 0.50 ~ "High",
      .data$n_IPF_detected == 3 & .data$min_IPF_candidate_fraction >= 0.25 &
        .data$median_IPF_jaccard >= 0.25 ~ "Moderate",
      TRUE ~ "Low"
    )
  )


# ------------------------------------------------------------
# 17. Structural signal in the control
# ------------------------------------------------------------

candidate_refinement <- candidate_refinement %>%
  dplyr::mutate(
    control_structure = dplyr::case_when(
      .data$control_detected == 0 ~ "No_overlap_detected",
      .data$control_candidate_fraction >= 0.50 ~ "Strong_overlap",
      .data$control_candidate_fraction >= 0.25 ~ "Moderate_overlap",
      TRUE ~ "Weak_overlap"
    )
  )


# ------------------------------------------------------------
# 18. Architecture classification
# ------------------------------------------------------------
# stable_single_locus      : single-TAR cluster, High reproducibility
# needs_BAM_refinement     : multi-TAR cluster, High reproducibility
# fragmented_complex_locus : multi-TAR cluster, Moderate/Low reproducibility
# poorly_reproducible      : single-TAR cluster, Moderate/Low reproducibility

candidate_refinement <- candidate_refinement %>%
  dplyr::mutate(
    structural_class = dplyr::case_when(
      .data$n_candidate_TAR_in_cluster == 1 & .data$structural_reproducibility == "High" ~
        "stable_single_locus",
      .data$n_candidate_TAR_in_cluster > 1 & .data$structural_reproducibility == "High" ~
        "needs_BAM_refinement",
      .data$n_candidate_TAR_in_cluster > 1 & .data$structural_reproducibility != "High" ~
        "fragmented_complex_locus",
      TRUE ~ "poorly_reproducible"
    )
  )


# ------------------------------------------------------------
# 19. Recommended next step
# ------------------------------------------------------------

candidate_refinement <- candidate_refinement %>%
  dplyr::mutate(
    next_step = dplyr::case_when(
      .data$structural_class == "stable_single_locus" ~
        "Structural annotation before coding-potential interpretation",
      .data$structural_class == "needs_BAM_refinement" ~
        "Recover BAM and resolve splice junctions and transcript boundaries",
      .data$structural_class == "fragmented_complex_locus" ~
        "Refine locus architecture before treating TARs as independent transcripts",
      TRUE ~
        "Lower structural priority; revisit after BAM recovery or additional samples"
    )
  )


# ------------------------------------------------------------
# 20. Sort and save full table
# ------------------------------------------------------------

candidate_refinement <- candidate_refinement %>%
  dplyr::arrange(
    factor(.data$structural_reproducibility, levels = c("High", "Moderate", "Low")),
    dplyr::desc(.data$min_IPF_candidate_fraction),
    dplyr::desc(.data$median_IPF_jaccard)
  )

write.csv(candidate_refinement, file.path(OUTPUT_DIR, "PriorityA_structural_refinement_complete.csv"),
          row.names = FALSE)

saveRDS(candidate_refinement, file.path(OUTPUT_DIR, "PriorityA_structural_refinement_complete.rds"))


# ------------------------------------------------------------
# 21. Subsets
# ------------------------------------------------------------

high_reproducibility     <- candidate_refinement %>% dplyr::filter(.data$structural_reproducibility == "High")
moderate_reproducibility <- candidate_refinement %>%
  dplyr::filter(.data$structural_reproducibility == "Moderate")
low_reproducibility      <- candidate_refinement %>% dplyr::filter(.data$structural_reproducibility == "Low")

stable_single       <- candidate_refinement %>% dplyr::filter(.data$structural_class == "stable_single_locus")
needs_BAM           <- candidate_refinement %>% dplyr::filter(.data$structural_class == "needs_BAM_refinement")
fragmented          <- candidate_refinement %>%
  dplyr::filter(.data$structural_class == "fragmented_complex_locus")
poorly_reproducible <- candidate_refinement %>% dplyr::filter(.data$structural_class == "poorly_reproducible")

no_control_overlap <- candidate_refinement %>% dplyr::filter(.data$control_detected == 0)

write.csv(high_reproducibility, file.path(OUTPUT_DIR, "PriorityA_high_structural_reproducibility.csv"),
          row.names = FALSE)
write.csv(moderate_reproducibility, file.path(OUTPUT_DIR, "PriorityA_moderate_structural_reproducibility.csv"),
          row.names = FALSE)
write.csv(low_reproducibility, file.path(OUTPUT_DIR, "PriorityA_low_structural_reproducibility.csv"),
          row.names = FALSE)
write.csv(stable_single, file.path(OUTPUT_DIR, "PriorityA_stable_single_locus.csv"), row.names = FALSE)
write.csv(needs_BAM, file.path(OUTPUT_DIR, "PriorityA_needs_BAM_refinement.csv"), row.names = FALSE)
write.csv(fragmented, file.path(OUTPUT_DIR, "PriorityA_fragmented_complex_locus.csv"), row.names = FALSE)
write.csv(poorly_reproducible, file.path(OUTPUT_DIR, "PriorityA_poorly_reproducible.csv"), row.names = FALSE)
write.csv(no_control_overlap, file.path(OUTPUT_DIR, "PriorityA_no_structural_overlap_control.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 22. Candidate x sample matrices (Jaccard, coverage, BED support)
# ------------------------------------------------------------

jaccard_matrix <- candidate_sample_complete %>%
  dplyr::select(TAR_ID, sample, jaccard) %>%
  tidyr::pivot_wider(names_from = sample, values_from = jaccard, values_fill = 0, names_prefix = "jaccard_")

write.csv(jaccard_matrix, file.path(OUTPUT_DIR, "PriorityA_Jaccard_matrix.csv"), row.names = FALSE)

coverage_matrix <- candidate_sample_complete %>%
  dplyr::select(TAR_ID, sample, candidate_fraction_covered) %>%
  tidyr::pivot_wider(names_from = sample, values_from = candidate_fraction_covered, values_fill = 0,
                     names_prefix = "coverage_")

write.csv(coverage_matrix, file.path(OUTPUT_DIR, "PriorityA_candidate_coverage_matrix.csv"), row.names = FALSE)

support_matrix <- candidate_sample_complete %>%
  dplyr::select(TAR_ID, sample, bed_support) %>%
  tidyr::pivot_wider(names_from = sample, values_from = bed_support, values_fill = 0, names_prefix = "support_")

write.csv(support_matrix, file.path(OUTPUT_DIR, "PriorityA_BED_support_matrix.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 23. IGV BED tracks (0-based starts)
# ------------------------------------------------------------

# All 48 Priority A candidates
priorityA_bed <- priority_A %>%
  dplyr::transmute(
    chr       = .data$chr,
    bed_start = pmax(.data$start - 1, 0),
    end       = .data$end,
    name      = paste0(.data$TAR_ID, "|", .data$cluster_id),
    score     = pmin(1000, round(.data$IPF_min_counts)),
    strand    = .data$strand
  )

write.table(priorityA_bed, file.path(IGV_DIR, "PriorityA_48_candidates.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# High reproducibility
high_bed <- high_reproducibility %>%
  dplyr::transmute(
    chr       = .data$chr,
    bed_start = pmax(.data$start - 1, 0),
    end       = .data$end,
    name      = paste0(.data$TAR_ID, "|", .data$structural_class),
    score     = pmin(1000, round(1000 * .data$median_IPF_jaccard)),
    strand    = .data$strand
  )

write.table(high_bed, file.path(IGV_DIR, "PriorityA_high_reproducibility.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# No structural overlap in the control
no_control_bed <- no_control_overlap %>%
  dplyr::transmute(
    chr       = .data$chr,
    bed_start = pmax(.data$start - 1, 0),
    end       = .data$end,
    name      = paste0(.data$TAR_ID, "|no_control_overlap"),
    score     = pmin(1000, round(1000 * .data$median_IPF_jaccard)),
    strand    = .data$strand
  )

write.table(no_control_bed, file.path(IGV_DIR, "PriorityA_no_structural_overlap_control.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 24. Summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 11 SUMMARY\n")
cat("====================================================\n\n")

cat("Priority A analysed:", nrow(candidate_refinement), "\n")
cat("High reproducibility:", nrow(high_reproducibility), "\n")
cat("Moderate reproducibility:", nrow(moderate_reproducibility), "\n")
cat("Low reproducibility:", nrow(low_reproducibility), "\n\n")
cat("stable_single_locus:", nrow(stable_single), "\n")
cat("needs_BAM_refinement:", nrow(needs_BAM), "\n")
cat("fragmented_complex_locus:", nrow(fragmented), "\n")
cat("poorly_reproducible:", nrow(poorly_reproducible), "\n\n")
cat("No structural overlap in control:", nrow(no_control_overlap), "\n\n")

cat("====================================================\n")
cat("CONTROL STRUCTURE\n")
cat("====================================================\n\n")

control_distribution <- candidate_refinement %>%
  dplyr::count(.data$control_structure) %>%
  dplyr::arrange(dplyr::desc(.data$n))

print(tibble::as_tibble(control_distribution), n = Inf, width = Inf)

cat("\n")
cat("====================================================\n")
cat("STRUCTURAL CLASSES\n")
cat("====================================================\n\n")

structural_distribution <- candidate_refinement %>%
  dplyr::count(.data$structural_reproducibility, .data$structural_class) %>%
  dplyr::arrange(.data$structural_reproducibility, .data$structural_class)

print(tibble::as_tibble(structural_distribution), n = Inf, width = Inf)

cat("\n")
cat("====================================================\n")
cat("HIGH STRUCTURAL REPRODUCIBILITY\n")
cat("====================================================\n\n")

high_table <- high_reproducibility %>%
  dplyr::select(
    TAR_ID, chr, start, end, strand,
    genomic_context, overlapping_gene_names,
    nearest_gene, nearest_gene_distance,
    cluster_id, n_candidate_TAR_in_cluster, n_PriorityA_in_cluster,
    min_IPF_candidate_fraction, median_IPF_candidate_fraction,
    min_IPF_jaccard, median_IPF_jaccard,
    start_range_IPF, end_range_IPF,
    control_candidate_fraction, control_jaccard,
    control_structure, structural_class
  ) %>%
  tibble::as_tibble()

print(high_table, n = Inf, width = Inf)

cat("\n")
cat("====================================================\n")
cat("NO STRUCTURAL OVERLAP IN CONTROL\n")
cat("====================================================\n\n")

no_control_table <- no_control_overlap %>%
  dplyr::select(
    TAR_ID, chr, start, end, strand,
    genomic_context, overlapping_gene_names,
    nearest_gene, nearest_gene_distance,
    min_IPF_candidate_fraction, median_IPF_jaccard,
    structural_reproducibility, structural_class
  ) %>%
  tibble::as_tibble()

print(no_control_table, n = Inf, width = Inf)

cat("\n")
cat("====================================================\n")
cat("PRIORITY A COMPACT TABLE\n")
cat("====================================================\n\n")

compact_table <- candidate_refinement %>%
  dplyr::select(
    TAR_ID, chr, start, end, strand,
    genomic_context, overlapping_gene_names,
    nearest_gene, nearest_gene_distance,
    cluster_id, n_candidate_TAR_in_cluster,
    n_IPF_detected,
    min_IPF_candidate_fraction, median_IPF_jaccard,
    control_detected, control_candidate_fraction,
    structural_reproducibility, structural_class
  ) %>%
  tibble::as_tibble()

print(compact_table, n = Inf, width = Inf)

cat("\n")
cat("====================================================\n")
cat("OUTPUT FILES\n")
cat("====================================================\n\n")

print(list.files(OUTPUT_DIR, recursive = TRUE, full.names = FALSE))


# ------------------------------------------------------------
# 25. Interpretation notes
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INTERPRETATION\n")
cat("====================================================\n\n")

cat(
  "High/Moderate/Low refers exclusively to the structural\n",
  "reproducibility of the TAR interval across\n",
  "the three IPF samples.\n\n",
  "Overlap in VUHD67 indicates that structural TAR\n",
  "signal exists at the same locus in the control.\n",
  "It does not imply equal expression between conditions.\n\n",
  "bed_support values must NOT be used directly\n",
  "to compare IPF vs control expression because of\n",
  "sequencing-depth differences between samples.\n\n",
  "Highly reproducible multi-TAR loci are kept\n",
  "as needs_BAM_refinement because splicing information\n",
  "is still needed to resolve transcript architecture.\n"
)

cat("\n")
cat("====================================================\n")
cat("STEP 11 COMPLETE\n")
cat("====================================================\n\n")

cat("Results saved to:\n", OUTPUT_DIR, "\n")
