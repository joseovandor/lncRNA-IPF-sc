# ============================================================
# 10_refine_focus_loci_tar_bed.R
#
# lncRNA-IPF-sc | Step 10: locus refinement of Priority A loci from TAR BED
# Author: Jose A. Ovando-Ricardez
#
# Refines the genomic architecture of Priority A loci using the per-sample
# TAR_reads.bed.gz files (no BAM required). Focus loci: ADNP antisense,
# DYRK1A antisense, TCF4 antisense and an intergenic region near SMAD4.
#
# TAR_reads.bed.gz has 7 columns (chr, start [0-based], end, ., ., strand,
# support) and carries no CIGAR, splice junctions or per-read cell barcodes.
# The script therefore evaluates spatial recurrence, TAR support, boundary
# reproducibility, cross-sample overlap and approximate locus continuity.
#
# Note: it does NOT reconstruct mature transcripts or isoforms.
#
# Inputs:  data/TAR_catalog/genomic_context/Unannotated_3of3_IPF_genomic_context.csv
#          data/<sample>/TAR/TAR_reads.bed.gz
# Outputs: data/TAR_catalog/locus_refinement_TARbed/ (CSV tables + IGV/ BED tracks)
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)

  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

INPUT_CONTEXT <- file.path(BASE_DIR, "data/TAR_catalog/genomic_context/Unannotated_3of3_IPF_genomic_context.csv")
OUTPUT_DIR    <- file.path(BASE_DIR, "data/TAR_catalog/locus_refinement_TARbed")
IGV_DIR       <- file.path(OUTPUT_DIR, "IGV")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(IGV_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")

conditions <- c(VUILD53 = "IPF", VUILD63 = "IPF", VUILD64 = "IPF", VUHD67 = "Control")


# ------------------------------------------------------------
# 3. TAR_reads.bed.gz paths
# ------------------------------------------------------------

tar_bed_paths <- setNames(file.path(BASE_DIR, "data", samples, "TAR", "TAR_reads.bed.gz"), samples)


# ------------------------------------------------------------
# 4. Check input files
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INPUT FILE CHECK\n")
cat("====================================================\n\n")

file_check <- data.frame(
  sample    = samples,
  condition = unname(conditions[samples]),
  TAR_bed   = unname(tar_bed_paths[samples]),
  exists    = file.exists(tar_bed_paths[samples]),
  stringsAsFactors = FALSE
)

print(file_check)

if (any(!file_check$exists)) {
  stop(paste("Missing TAR_reads.bed.gz for:", paste(file_check$sample[!file_check$exists], collapse = ", ")))
}

if (!file.exists(INPUT_CONTEXT)) {
  stop(paste("Genomic context file not found:", INPUT_CONTEXT))
}


# ------------------------------------------------------------
# 5. Load genomic context results (previous step)
# ------------------------------------------------------------

result <- read.csv(INPUT_CONTEXT, stringsAsFactors = FALSE, check.names = FALSE)

cat("\n")
cat("====================================================\n")
cat("INPUT\n")
cat("====================================================\n\n")

cat("Total candidates:", nrow(result), "\n")
cat("Priority A:", sum(result$priority_group == "Priority_A", na.rm = TRUE), "\n\n")


# ------------------------------------------------------------
# 6. Required columns
# ------------------------------------------------------------

required_columns <- c("TAR_ID", "chr", "start", "end", "strand", "length", "priority_group",
                      "cluster_id", "cluster_n_TAR", "genomic_context", "overlapping_gene_names",
                      "IPF_min_counts", "IPF_CV_counts", "IPF_min_pct_cells")

missing_columns <- setdiff(required_columns, colnames(result))

if (length(missing_columns) > 0) {
  stop(paste("Missing columns in INPUT_CONTEXT:", paste(missing_columns, collapse = ", ")))
}


# ------------------------------------------------------------
# 7. Priority A
# ------------------------------------------------------------

priority_A <- result %>%
  dplyr::filter(.data$priority_group == "Priority_A")

cat("Priority A loaded:", nrow(priority_A), "\n")

if (nrow(priority_A) == 0) {
  stop("No Priority A candidates found.")
}


# ------------------------------------------------------------
# 8. Priority A locus catalogue
# ------------------------------------------------------------

priority_A_loci <- priority_A %>%
  dplyr::group_by(.data$cluster_id) %>%
  dplyr::summarise(
    chr                = dplyr::first(.data$chr),
    locus_start        = min(.data$start, na.rm = TRUE),
    locus_end          = max(.data$end, na.rm = TRUE),
    strand             = dplyr::first(.data$strand),
    n_priorityA_TAR    = dplyr::n(),
    priorityA_members  = paste(.data$TAR_ID, collapse = ";"),
    genomic_context    = paste(unique(.data$genomic_context), collapse = ";"),
    overlapping_genes  = paste(unique(na.omit(.data$overlapping_gene_names)), collapse = ";"),
    max_IPF_min_counts = max(.data$IPF_min_counts, na.rm = TRUE),
    min_CV             = min(.data$IPF_CV_counts, na.rm = TRUE),
    max_min_pct_cells  = max(.data$IPF_min_pct_cells, na.rm = TRUE),
    .groups = "drop"
  )

priority_A_loci$cluster_total_TAR <- result$cluster_n_TAR[match(priority_A_loci$cluster_id, result$cluster_id)]

cat("Unique Priority A loci:", nrow(priority_A_loci), "\n\n")

write.csv(priority_A_loci, file.path(OUTPUT_DIR, "Priority_A_locus_catalog.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 9. Focus loci (seed TARs)
# ------------------------------------------------------------

focus_seeds <- data.frame(
  locus_name = c("ADNP_antisense", "DYRK1A_antisense", "TCF4_antisense", "SMAD4_near_intergenic"),
  seed_TAR   = c("IPFTAR049318", "IPFTAR053241", "IPFTAR032367", "IPFTAR032288"),
  stringsAsFactors = FALSE
)

missing_seeds <- setdiff(focus_seeds$seed_TAR, result$TAR_ID)

if (length(missing_seeds) > 0) {
  stop(paste("Seed TARs not found:", paste(missing_seeds, collapse = ", ")))
}


# ------------------------------------------------------------
# 10. Assign cluster to each locus
# ------------------------------------------------------------

focus_seeds$cluster_id <- result$cluster_id[match(focus_seeds$seed_TAR, result$TAR_ID)]

if (any(is.na(focus_seeds$cluster_id))) {
  stop("At least one seed TAR has no cluster_id.")
}

cat("\n")
cat("====================================================\n")
cat("FOCUS LOCI\n")
cat("====================================================\n\n")

print(focus_seeds)


# ------------------------------------------------------------
# 11. Retrieve all TARs of each cluster
# ------------------------------------------------------------

focus_members_list <- list()

for (i in seq_len(nrow(focus_seeds))) {
  current_cluster <- focus_seeds$cluster_id[i]
  current_locus <- focus_seeds$locus_name[i]

  tmp <- result %>%
    dplyr::filter(.data$cluster_id == current_cluster)

  tmp$locus_name <- current_locus

  focus_members_list[[i]] <- tmp
}

focus_members <- dplyr::bind_rows(focus_members_list)

focus_members <- focus_members %>%
  dplyr::arrange(.data$locus_name, .data$start)

if (nrow(focus_members) == 0) {
  stop("No TARs retrieved from the selected clusters.")
}


# ------------------------------------------------------------
# 12. Show TARs of the focus loci
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CLUSTER TARs\n")
cat("====================================================\n\n")

print(
  focus_members %>%
    dplyr::select(locus_name, TAR_ID, chr, start, end, strand, length, priority_group, genomic_context,
                  overlapping_gene_names, IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells)
)

write.csv(focus_members, file.path(OUTPUT_DIR, "Focus_loci_consensus_TAR.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 13. Locus boundaries
# ------------------------------------------------------------

focus_loci <- focus_members %>%
  dplyr::group_by(.data$locus_name, .data$cluster_id, .data$chr, .data$strand) %>%
  dplyr::summarise(
    locus_start     = min(.data$start, na.rm = TRUE),
    locus_end       = max(.data$end, na.rm = TRUE),
    n_consensus_TAR = dplyr::n(),
    consensus_TAR   = paste(.data$TAR_ID, collapse = ";"),
    .groups = "drop"
  )

print(focus_loci)


# ------------------------------------------------------------
# 14. Locus GRanges
# ------------------------------------------------------------

focus_gr <- GRanges(
  seqnames = focus_loci$chr,
  ranges   = IRanges(start = focus_loci$locus_start, end = focus_loci$locus_end),
  strand   = focus_loci$strand
)

names(focus_gr) <- focus_loci$locus_name


# ------------------------------------------------------------
# 15. +/-25 kb windows
# ------------------------------------------------------------

focus_window <- resize(focus_gr, width = width(focus_gr) + 50000, fix = "center")

start(focus_window) <- pmax(start(focus_window), 1)


# ------------------------------------------------------------
# 16. Consensus TAR GRanges
# ------------------------------------------------------------

candidate_gr <- GRanges(
  seqnames = focus_members$chr,
  ranges   = IRanges(start = focus_members$start, end = focus_members$end),
  strand   = focus_members$strand
)

names(candidate_gr) <- focus_members$TAR_ID

mcols(candidate_gr)$locus_name <- focus_members$locus_name


# ------------------------------------------------------------
# 17. Reader for TAR_reads.bed.gz
#     BED start is 0-based, GRanges start is 1-based:
#     start = bed_start + 1; end = bed_end
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


# ------------------------------------------------------------
# 18. Read TAR BED for all samples
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("READING TAR_reads.bed.gz\n")
cat("====================================================\n\n")

sample_tar_list <- list()

for (sample_name in samples) {
  cat("Reading:", sample_name, "\n")

  tmp <- read_tar_bed(path = tar_bed_paths[sample_name], sample_name = sample_name,
                      condition_name = unname(conditions[sample_name]))

  cat("  intervals:", nrow(tmp), "\n")

  sample_tar_list[[sample_name]] <- tmp
}


# ------------------------------------------------------------
# 19. Per-sample TAR BED summary
# ------------------------------------------------------------

sample_summary_list <- list()

for (sample_name in samples) {
  tmp <- sample_tar_list[[sample_name]]

  sample_summary_list[[sample_name]] <- data.frame(
    sample         = sample_name,
    condition      = unname(conditions[sample_name]),
    n_intervals    = nrow(tmp),
    total_support  = sum(tmp$bed_support, na.rm = TRUE),
    median_support = median(tmp$bed_support, na.rm = TRUE),
    median_length  = median(tmp$interval_length, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

sample_summary <- dplyr::bind_rows(sample_summary_list)

cat("\n")
cat("====================================================\n")
cat("TAR BED SUMMARY PER SAMPLE\n")
cat("====================================================\n\n")

print(sample_summary)

write.csv(sample_summary, file.path(OUTPUT_DIR, "TAR_reads_BED_sample_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 20. Extract local TARs within +/-25 kb
# ------------------------------------------------------------

local_tar_list <- list()
local_counter <- 1

for (sample_name in samples) {
  sample_df <- sample_tar_list[[sample_name]]

  sample_gr <- GRanges(
    seqnames = sample_df$chr,
    ranges   = IRanges(start = sample_df$start, end = sample_df$end),
    strand   = sample_df$strand
  )

  mcols(sample_gr)$row_id <- seq_len(nrow(sample_df))

  for (current_locus in names(focus_window)) {
    window_gr <- focus_window[current_locus]

    hits <- findOverlaps(sample_gr, window_gr, ignore.strand = FALSE)

    if (length(hits) > 0) {
      idx <- unique(queryHits(hits))

      tmp <- sample_df[idx, , drop = FALSE]

      tmp$locus_name <- current_locus

      local_tar_list[[local_counter]] <- tmp

      local_counter <- local_counter + 1
    }
  }
}


# ------------------------------------------------------------
# 21. Combine local TARs
# ------------------------------------------------------------

if (length(local_tar_list) == 0) {
  stop("No TARs recovered within the +/-25 kb windows.")
}

local_tar <- dplyr::bind_rows(local_tar_list)

cat("\n")
cat("====================================================\n")
cat("LOCAL TARs +/-25 kb\n")
cat("====================================================\n\n")

cat("Local intervals:", nrow(local_tar), "\n\n")

print(local_tar %>% dplyr::count(locus_name, sample, condition, name = "n_TAR"))

write.csv(local_tar, file.path(OUTPUT_DIR, "Focus_loci_sample_TAR_25kb.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 22. Local summary per locus / sample
# ------------------------------------------------------------

local_summary <- local_tar %>%
  dplyr::group_by(.data$locus_name, .data$sample, .data$condition) %>%
  dplyr::summarise(
    n_TAR          = dplyr::n(),
    total_support  = sum(.data$bed_support, na.rm = TRUE),
    max_support    = max(.data$bed_support, na.rm = TRUE),
    median_support = median(.data$bed_support, na.rm = TRUE),
    local_start    = min(.data$start, na.rm = TRUE),
    local_end      = max(.data$end, na.rm = TRUE),
    local_span     = .data$local_end - .data$local_start + 1,
    .groups = "drop"
  )

write.csv(local_summary, file.path(OUTPUT_DIR, "Focus_loci_sample_summary_25kb.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 23. Local TARs as GRanges
# ------------------------------------------------------------

local_gr <- GRanges(
  seqnames = local_tar$chr,
  ranges   = IRanges(start = local_tar$start, end = local_tar$end),
  strand   = local_tar$strand
)

mcols(local_gr)$sample <- local_tar$sample
mcols(local_gr)$condition <- local_tar$condition
mcols(local_gr)$bed_support <- local_tar$bed_support
mcols(local_gr)$locus_name <- local_tar$locus_name
mcols(local_gr)$local_row <- seq_len(nrow(local_tar))


# ------------------------------------------------------------
# 24. Overlap: consensus TARs vs sample TARs
# ------------------------------------------------------------

candidate_hits <- findOverlaps(candidate_gr, local_gr, ignore.strand = FALSE)

if (length(candidate_hits) == 0) {
  stop("No overlaps between consensus TARs and local TARs.")
}


# ------------------------------------------------------------
# 25. Keep overlaps within the same locus
# ------------------------------------------------------------

query_index <- queryHits(candidate_hits)
subject_index <- subjectHits(candidate_hits)

same_locus <- as.character(mcols(candidate_gr)$locus_name[query_index]) ==
  as.character(mcols(local_gr)$locus_name[subject_index])

candidate_hits <- candidate_hits[same_locus]

if (length(candidate_hits) == 0) {
  stop("No overlaps left after restricting to the same locus.")
}


# ------------------------------------------------------------
# 26. Overlap metrics
# ------------------------------------------------------------

overlap_list <- list()

for (i in seq_along(candidate_hits)) {
  q <- queryHits(candidate_hits)[i]
  s <- subjectHits(candidate_hits)[i]

  candidate <- candidate_gr[q]
  sample_tar <- local_gr[s]

  intersection_start <- max(start(candidate), start(sample_tar))
  intersection_end <- min(end(candidate), end(sample_tar))

  overlap_bp <- max(0, intersection_end - intersection_start + 1)

  candidate_width <- as.numeric(width(candidate))
  sample_width <- as.numeric(width(sample_tar))

  union_bp <- candidate_width + sample_width - overlap_bp

  overlap_list[[i]] <- data.frame(
    locus_name                  = as.character(mcols(candidate)$locus_name),
    TAR_ID                      = names(candidate),
    sample                      = as.character(mcols(sample_tar)$sample),
    condition                   = as.character(mcols(sample_tar)$condition),
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
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 27. Combine overlaps
# ------------------------------------------------------------

candidate_overlap <- dplyr::bind_rows(overlap_list)

if (nrow(candidate_overlap) == 0) {
  stop("candidate_overlap is empty.")
}

write.csv(candidate_overlap, file.path(OUTPUT_DIR, "Consensus_TAR_vs_sample_TAR_overlaps.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 28. Support summary per consensus TAR / sample
# ------------------------------------------------------------

candidate_sample_summary <- candidate_overlap %>%
  dplyr::group_by(.data$locus_name, .data$TAR_ID, .data$sample, .data$condition) %>%
  dplyr::summarise(
    n_overlapping_sample_TAR        = dplyr::n(),
    total_bed_support               = sum(.data$bed_support, na.rm = TRUE),
    max_bed_support                 = max(.data$bed_support, na.rm = TRUE),
    max_overlap_bp                  = max(.data$overlap_bp, na.rm = TRUE),
    max_candidate_fraction_covered  = max(.data$candidate_fraction_covered, na.rm = TRUE),
    max_sample_TAR_fraction_covered = max(.data$sample_TAR_fraction_covered, na.rm = TRUE),
    max_jaccard                     = max(.data$jaccard, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(candidate_sample_summary, file.path(OUTPUT_DIR, "Consensus_TAR_sample_support_summary.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 29. Support matrix
# ------------------------------------------------------------

support_wide <- candidate_sample_summary %>%
  dplyr::select(locus_name, TAR_ID, sample, total_bed_support) %>%
  tidyr::pivot_wider(names_from = sample, values_from = total_bed_support, values_fill = 0,
                     names_prefix = "support_")

write.csv(support_wide, file.path(OUTPUT_DIR, "Consensus_TAR_support_matrix.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 30. Best overlap per candidate / sample
# ------------------------------------------------------------

best_overlap <- candidate_overlap %>%
  dplyr::group_by(.data$locus_name, .data$TAR_ID, .data$sample, .data$condition) %>%
  dplyr::arrange(dplyr::desc(.data$jaccard), dplyr::desc(.data$overlap_bp), .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

write.csv(best_overlap, file.path(OUTPUT_DIR, "Consensus_TAR_best_sample_match.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 31. Recurrent local segments
#     Only overlapping TARs are merged (min.gapwidth = 1);
#     merely nearby TARs are not joined.
# ------------------------------------------------------------

local_segments_list <- list()
segment_counter <- 1

locus_names <- unique(local_tar$locus_name)

for (current_locus in locus_names) {
  locus_df <- local_tar %>%
    dplyr::filter(.data$locus_name == current_locus)

  if (nrow(locus_df) == 0) {
    next
  }

  locus_gr <- GRanges(
    seqnames = locus_df$chr,
    ranges   = IRanges(start = locus_df$start, end = locus_df$end),
    strand   = locus_df$strand
  )

  locus_reduced <- reduce(locus_gr, ignore.strand = FALSE, min.gapwidth = 1)

  segment_names <- paste0(current_locus, "_SEG", sprintf("%03d", seq_along(locus_reduced)))

  names(locus_reduced) <- segment_names

  hits <- findOverlaps(locus_reduced, locus_gr, ignore.strand = FALSE)

  for (i in seq_along(locus_reduced)) {
    source_idx <- subjectHits(hits)[queryHits(hits) == i]

    source_samples <- locus_df$sample[source_idx]
    source_conditions <- locus_df$condition[source_idx]
    source_support <- locus_df$bed_support[source_idx]

    ipf_samples_here <- unique(source_samples[source_conditions == "IPF"])
    control_samples_here <- unique(source_samples[source_conditions == "Control"])

    local_segments_list[[segment_counter]] <- data.frame(
      segment_id    = names(locus_reduced)[i],
      locus_name    = current_locus,
      chr           = as.character(seqnames(locus_reduced[i])),
      start         = as.numeric(start(locus_reduced[i])),
      end           = as.numeric(end(locus_reduced[i])),
      strand        = as.character(strand(locus_reduced[i])),
      length        = as.numeric(width(locus_reduced[i])),
      n_source_TAR  = length(source_idx),
      n_samples     = length(unique(source_samples)),
      n_IPF         = length(ipf_samples_here),
      n_Control     = length(control_samples_here),
      samples       = paste(sort(unique(source_samples)), collapse = ";"),
      total_support = sum(source_support, na.rm = TRUE),
      stringsAsFactors = FALSE
    )

    segment_counter <- segment_counter + 1
  }
}


# ------------------------------------------------------------
# 32. Combine segments
# ------------------------------------------------------------

if (length(local_segments_list) == 0) {
  stop("Could not build local segments.")
}

local_segments <- dplyr::bind_rows(local_segments_list)

local_segments <- local_segments %>%
  dplyr::mutate(
    recurrence_class = dplyr::case_when(
      .data$n_IPF == 3 & .data$n_Control == 0 ~ "3of3_IPF_0of1_Control",
      .data$n_IPF == 3 & .data$n_Control == 1 ~ "3of3_IPF_and_Control",
      .data$n_IPF > 0 & .data$n_Control == 0 ~ "IPF_only_pilot",
      TRUE ~ "Other"
    )
  )

write.csv(local_segments, file.path(OUTPUT_DIR, "Focus_loci_reduced_local_segments.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 33. Segments recurrent in 3/3 IPF
# ------------------------------------------------------------

segments_3of3 <- local_segments %>%
  dplyr::filter(.data$n_IPF == 3)

segments_3of3_no_control <- local_segments %>%
  dplyr::filter(.data$n_IPF == 3, .data$n_Control == 0)

write.csv(segments_3of3, file.path(OUTPUT_DIR, "Focus_loci_segments_3of3_IPF.csv"), row.names = FALSE)

write.csv(segments_3of3_no_control, file.path(OUTPUT_DIR, "Focus_loci_segments_3of3_IPF_0of1_control.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 34. Gaps between consensus TARs of the same locus
# ------------------------------------------------------------

gap_list <- list()
gap_counter <- 1

focus_locus_names <- unique(focus_members$locus_name)

for (current_locus in focus_locus_names) {
  tmp <- focus_members %>%
    dplyr::filter(.data$locus_name == current_locus) %>%
    dplyr::arrange(.data$start)

  if (nrow(tmp) < 2) {
    next
  }

  for (i in seq_len(nrow(tmp) - 1)) {
    gap_bp <- tmp$start[i + 1] - tmp$end[i] - 1

    gap_list[[gap_counter]] <- data.frame(
      locus_name  = current_locus,
      left_TAR    = tmp$TAR_ID[i],
      right_TAR   = tmp$TAR_ID[i + 1],
      left_end    = tmp$end[i],
      right_start = tmp$start[i + 1],
      gap_bp      = gap_bp,
      stringsAsFactors = FALSE
    )

    gap_counter <- gap_counter + 1
  }
}

if (length(gap_list) > 0) {
  candidate_gaps <- dplyr::bind_rows(gap_list)
} else {
  candidate_gaps <- data.frame(locus_name = character(), left_TAR = character(), right_TAR = character(),
                               left_end = numeric(), right_start = numeric(), gap_bp = numeric())
}

write.csv(candidate_gaps, file.path(OUTPUT_DIR, "Focus_loci_consensus_TAR_gaps.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 35. BED: consensus TARs
# ------------------------------------------------------------

consensus_bed <- focus_members %>%
  dplyr::transmute(
    chr       = .data$chr,
    bed_start = pmax(.data$start - 1, 0),
    end       = .data$end,
    name      = paste0(.data$locus_name, "|", .data$TAR_ID),
    score     = pmin(1000, round(.data$IPF_min_counts)),
    strand    = .data$strand
  )

write.table(consensus_bed, file = file.path(IGV_DIR, "Focus_consensus_TAR.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 36. BED: local TARs per sample
# ------------------------------------------------------------

for (sample_name in samples) {
  tmp <- local_tar %>%
    dplyr::filter(.data$sample == sample_name) %>%
    dplyr::transmute(
      chr       = .data$chr,
      bed_start = pmax(.data$start - 1, 0),
      end       = .data$end,
      name      = paste0(.data$locus_name, "|", sample_name, "|support=", .data$bed_support),
      score     = pmin(1000, round(.data$bed_support)),
      strand    = .data$strand
    )

  write.table(tmp, file = file.path(IGV_DIR, paste0(sample_name, "_Focus_local_TAR.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# ------------------------------------------------------------
# 37. BED: recurrent segments
# ------------------------------------------------------------

segments_bed <- local_segments %>%
  dplyr::transmute(
    chr       = .data$chr,
    bed_start = pmax(.data$start - 1, 0),
    end       = .data$end,
    name      = paste0(.data$segment_id, "|", .data$recurrence_class),
    score     = pmin(1000, round(.data$total_support)),
    strand    = .data$strand
  )

write.table(segments_bed, file = file.path(IGV_DIR, "Focus_reduced_local_segments.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 38. BED: +/-25 kb windows
# ------------------------------------------------------------

window_bed <- data.frame(
  chr    = as.character(seqnames(focus_window)),
  start  = pmax(start(focus_window) - 1, 0),
  end    = end(focus_window),
  name   = names(focus_window),
  score  = 0,
  strand = as.character(strand(focus_window)),
  stringsAsFactors = FALSE
)

write.table(window_bed, file = file.path(IGV_DIR, "Focus_loci_windows_25kb.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 39. Summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUMMARY\n")
cat("====================================================\n\n")

cat("Priority A TAR:", nrow(priority_A), "\n")
cat("Priority A loci:", nrow(priority_A_loci), "\n")
cat("Focus loci analysed:", nrow(focus_loci), "\n")
cat("Consensus TARs in those loci:", nrow(focus_members), "\n")
cat("Local TARs recovered:", nrow(local_tar), "\n")
cat("Consensus/sample overlaps:", nrow(candidate_overlap), "\n")
cat("Reduced local segments:", nrow(local_segments), "\n")
cat("Segments present in 3/3 IPF:", nrow(segments_3of3), "\n")
cat("Segments 3/3 IPF + 0/1 control:", nrow(segments_3of3_no_control), "\n\n")


# ------------------------------------------------------------
# 40. Local support per sample
# ------------------------------------------------------------

cat("====================================================\n")
cat("LOCAL SUPPORT PER SAMPLE\n")
cat("====================================================\n\n")

print(local_summary %>% dplyr::arrange(.data$locus_name, .data$sample))


# ------------------------------------------------------------
# 41. Consensus TAR support
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CONSENSUS TAR SUPPORT PER SAMPLE\n")
cat("====================================================\n\n")

print(candidate_sample_summary %>% dplyr::arrange(.data$locus_name, .data$TAR_ID, .data$sample))


# ------------------------------------------------------------
# 42. Best matches
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("BEST SAMPLE TAR FOR EACH CONSENSUS TAR\n")
cat("====================================================\n\n")

print(
  best_overlap %>%
    dplyr::select(locus_name, TAR_ID, sample, condition, candidate_start, candidate_end, sample_TAR_start,
                  sample_TAR_end, bed_support, overlap_bp, candidate_fraction_covered,
                  sample_TAR_fraction_covered, jaccard) %>%
    dplyr::arrange(.data$locus_name, .data$TAR_ID, .data$sample)
)


# ------------------------------------------------------------
# 43. Recurrent segments
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RECURRENT LOCAL SEGMENTS\n")
cat("====================================================\n\n")

print(local_segments %>%
        dplyr::arrange(.data$locus_name, dplyr::desc(.data$n_IPF), dplyr::desc(.data$total_support)))


# ------------------------------------------------------------
# 44. Gaps between consensus TARs
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("DISTANCE BETWEEN CONSENSUS TARs OF THE SAME LOCUS\n")
cat("====================================================\n\n")

if (nrow(candidate_gaps) > 0) {
  print(candidate_gaps)
} else {
  cat("No loci with more than one TAR.\n")
}


# ------------------------------------------------------------
# 45. Per-locus detail
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("ADNP ANTISENSE\n")
cat("====================================================\n\n")

print(best_overlap %>% dplyr::filter(.data$locus_name == "ADNP_antisense") %>%
        dplyr::arrange(.data$TAR_ID, .data$sample))

cat("\n")
cat("====================================================\n")
cat("DYRK1A ANTISENSE\n")
cat("====================================================\n\n")

print(best_overlap %>% dplyr::filter(.data$locus_name == "DYRK1A_antisense") %>%
        dplyr::arrange(.data$TAR_ID, .data$sample))

cat("\n")
cat("====================================================\n")
cat("TCF4 ANTISENSE\n")
cat("====================================================\n\n")

print(best_overlap %>% dplyr::filter(.data$locus_name == "TCF4_antisense") %>%
        dplyr::arrange(.data$TAR_ID, .data$sample))

cat("\n")
cat("====================================================\n")
cat("SMAD4-NEAR INTERGENIC\n")
cat("====================================================\n\n")

print(best_overlap %>% dplyr::filter(.data$locus_name == "SMAD4_near_intergenic") %>%
        dplyr::arrange(.data$TAR_ID, .data$sample))


# ------------------------------------------------------------
# 46. Output files
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("OUTPUT FILES\n")
cat("====================================================\n\n")

print(list.files(OUTPUT_DIR, recursive = TRUE, full.names = FALSE))


# ------------------------------------------------------------
# 47. Interpretation
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INTERPRETATION\n")
cat("====================================================\n\n")

cat(
  "TAR_reads.bed.gz contains aggregated TAR intervals.\n\n",
  "This analysis evaluates:\n",
  "  - spatial recurrence\n",
  "  - boundary reproducibility\n",
  "  - overlap between samples\n",
  "  - relative interval support\n",
  "  - approximate locus continuity\n\n",
  "This analysis does NOT directly determine:\n",
  "  - splice junctions\n",
  "  - exon-intron structure\n",
  "  - mature isoforms\n",
  "  - full transcript structure\n\n",
  "Reduced segments represent recurrent TAR regions,\n",
  "NOT demonstrated mature transcripts.\n"
)


# ------------------------------------------------------------
# 48. Done
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOCUS REFINEMENT FINISHED\n")
cat("====================================================\n\n")

cat("Results saved to:\n", OUTPUT_DIR, "\n")
