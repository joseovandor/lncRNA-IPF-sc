# ============================================================
# 08_prioritize_candidates_and_extract_sequences.R
#
# lncRNA-IPF-sc | Step 8: candidate prioritization and sequence extraction
# Author: Jose A. Ovando-Ricardez
#
# Prioritizes the unannotated candidates detected in 3/3 IPF samples
# (by control detection, expression level and consistency) and extracts
# their strand-specific genomic sequences (GRCh38) with basic sequence
# metrics, as input for coding-potential analysis.
#
# Note: sequences are the full genomic TAR interval, not necessarily
# the mature transcript.
#
# Requires: uncompressed, indexed GRCh38 FASTA in data/reference/
#           (Homo_sapiens.GRCh38.dna.primary_assembly.fa; .fai is created if missing)
#
# Inputs:  data/TAR_catalog/quantification_3of3_IPF/Unannotated_lncRNA_candidates_3of3_IPF.csv
#          data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# Outputs: data/TAR_catalog/coding_potential/
#          (master table CSV/RDS, per-priority CSVs, BED and genomic FASTA files)
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(Rsamtools)
  library(Biostrings)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

INPUT_FILE <- file.path(
  BASE_DIR, "data/TAR_catalog/quantification_3of3_IPF/Unannotated_lncRNA_candidates_3of3_IPF.csv"
)
OUTPUT_DIR      <- file.path(BASE_DIR, "data/TAR_catalog/coding_potential")
REFERENCE_DIR   <- file.path(BASE_DIR, "data/reference")
REFERENCE_FASTA <- file.path(REFERENCE_DIR, "Homo_sapiens.GRCh38.dna.primary_assembly.fa")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REFERENCE_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Check reference genome
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GENOME REFERENCE\n")
cat("====================================================\n\n")

cat("Expected FASTA:\n", REFERENCE_FASTA, "\n\n")

if (!file.exists(REFERENCE_FASTA)) {
  stop(
    paste0(
      "\nUncompressed FASTA not found:\n",
      REFERENCE_FASTA,
      "\n\n",
      "Before running this script, run in a terminal:\n\n",
      "mkdir -p /home/jovando/IPF_lncRNA_pilot/data/reference\n\n",
      "gzip -dc /app/project/compartido/genomes/",
      "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ",
      "> /home/jovando/IPF_lncRNA_pilot/data/reference/",
      "Homo_sapiens.GRCh38.dna.primary_assembly.fa\n\n",
      "samtools faidx /home/jovando/IPF_lncRNA_pilot/data/reference/",
      "Homo_sapiens.GRCh38.dna.primary_assembly.fa\n"
    )
  )
}


# ------------------------------------------------------------
# 3. Check input
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(paste("Candidate file not found:", INPUT_FILE))
}


# ------------------------------------------------------------
# 4. Load candidates
# ------------------------------------------------------------

candidates <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

cat("\n")
cat("====================================================\n")
cat("INPUT CANDIDATES\n")
cat("====================================================\n\n")

cat("Number of candidates:", nrow(candidates), "\n\n")


# ------------------------------------------------------------
# 5. Required columns
# ------------------------------------------------------------

required_columns <- c(
  "TAR_ID",
  "chr", "start", "end", "strand", "length",
  "n_IPF", "n_Control",
  "IPF_mean_counts", "IPF_min_counts", "IPF_CV_counts",
  "IPF_mean_pct_cells", "IPF_min_pct_cells",
  "counts_VUILD53", "counts_VUILD63", "counts_VUILD64", "counts_VUHD67"
)

missing_columns <- setdiff(required_columns, colnames(candidates))

if (length(missing_columns) > 0) {
  stop(paste("Missing columns:", paste(missing_columns, collapse = ", ")))
}


# ------------------------------------------------------------
# 6. Basic checks
# ------------------------------------------------------------

cat("All n_IPF == 3:", all(candidates$n_IPF == 3), "\n")
cat("All >=200 bp:", all(candidates$length >= 200), "\n")
cat("Duplicated IDs:", sum(duplicated(candidates$TAR_ID)), "\n\n")


# ------------------------------------------------------------
# 7. Control detection
# ------------------------------------------------------------

candidates$pilot_control_status <- ifelse(
  candidates$n_Control == 0, "Not_detected_control", "Detected_control"
)


# ------------------------------------------------------------
# 8. Expression support
# ------------------------------------------------------------

candidates$expression_support <- case_when(
  candidates$IPF_min_counts >= 100 & candidates$IPF_CV_counts < 1 ~ "Very_high_support",
  candidates$IPF_min_counts >= 10 & candidates$IPF_CV_counts < 1 ~ "High_support",
  candidates$IPF_min_counts >= 3 & candidates$IPF_CV_counts < 1 ~ "Moderate_support",
  TRUE ~ "Low_or_heterogeneous"
)


# ------------------------------------------------------------
# 9. Pilot priority
# ------------------------------------------------------------

candidates$priority_group <- case_when(
  candidates$n_Control == 0 & candidates$IPF_min_counts >= 10 & candidates$IPF_CV_counts < 1 ~ "Priority_A",
  candidates$n_Control == 0 ~ "Priority_B",
  candidates$n_Control > 0 ~ "Priority_C",
  TRUE ~ "Unclassified"
)


# ------------------------------------------------------------
# 10. High-confidence pilot
# ------------------------------------------------------------

candidates$high_confidence_pilot <-
  candidates$n_Control == 0 &
  candidates$IPF_min_counts >= 10 &
  candidates$IPF_CV_counts < 1 &
  candidates$IPF_min_pct_cells >= 0.1


# ------------------------------------------------------------
# 11. Priority summary
# ------------------------------------------------------------

cat("====================================================\n")
cat("PRIORITIES\n")
cat("====================================================\n\n")

print(table(candidates$priority_group))
cat("\n")
print(table(candidates$expression_support))
cat("\n")

cat("High-confidence pilot:", sum(candidates$high_confidence_pilot, na.rm = TRUE), "\n\n")


# ------------------------------------------------------------
# 12. GRanges
# ------------------------------------------------------------

candidate_gr <- GRanges(
  seqnames = candidates$chr,
  ranges = IRanges(start = candidates$start, end = candidates$end),
  strand = candidates$strand
)

names(candidate_gr) <- candidates$TAR_ID

mcols(candidate_gr)$TAR_ID         <- candidates$TAR_ID
mcols(candidate_gr)$priority_group <- candidates$priority_group
mcols(candidate_gr)$n_Control      <- candidates$n_Control
mcols(candidate_gr)$IPF_min_counts <- candidates$IPF_min_counts
mcols(candidate_gr)$IPF_CV_counts  <- candidates$IPF_CV_counts


# ------------------------------------------------------------
# 13. Convert chr* to Ensembl style
# ------------------------------------------------------------
# Catalogue uses chr1..chrX/chrY; the Ensembl FASTA uses 1..X/Y.

cat("\n")
cat("====================================================\n")
cat("CHROMOSOME CONVERSION\n")
cat("====================================================\n\n")

old_seqlevels <- GenomeInfoDb::seqlevels(candidate_gr)

new_seqlevels <- sub("^chr", "", old_seqlevels)
new_seqlevels[new_seqlevels == "M"] <- "MT"

seqlevel_conversion <- data.frame(original = old_seqlevels, ensembl = new_seqlevels, stringsAsFactors = FALSE)

cat("Seqlevel conversion:\n")
print(seqlevel_conversion)

rename_map <- setNames(new_seqlevels, old_seqlevels)

candidate_gr <- GenomeInfoDb::renameSeqlevels(candidate_gr, rename_map)

cat("\nSeqlevels after conversion:\n")
print(GenomeInfoDb::seqlevels(candidate_gr))

cat("\nExample seqnames:\n")
print(head(unique(as.character(seqnames(candidate_gr))), 30))


# ------------------------------------------------------------
# 14. Prepare FASTA
# ------------------------------------------------------------

fa <- FaFile(REFERENCE_FASTA)

FAI_FILE <- paste0(REFERENCE_FASTA, ".fai")

if (!file.exists(FAI_FILE)) {
  cat("\nNo .fai index found.\n")
  cat("Creating index with Rsamtools::indexFa()...\n")
  Rsamtools::indexFa(REFERENCE_FASTA)
}


# ------------------------------------------------------------
# 15. Read FASTA index
# ------------------------------------------------------------

fasta_index <- scanFaIndex(fa)

fasta_seqnames <- as.character(seqnames(fasta_index))
tar_seqnames   <- unique(as.character(seqnames(candidate_gr)))

common_chr  <- intersect(tar_seqnames, fasta_seqnames)
missing_chr <- setdiff(tar_seqnames, fasta_seqnames)

cat("\n")
cat("====================================================\n")
cat("CHROMOSOME CHECK\n")
cat("====================================================\n\n")

cat("Candidate chromosomes:", length(tar_seqnames), "\n")
cat("Chromosomes shared with FASTA:", length(common_chr), "\n")

if (length(missing_chr) > 0) {
  cat("\nTAR chromosomes missing from FASTA:\n")
  print(missing_chr)
}


# ------------------------------------------------------------
# 16. Keep regions on available chromosomes
# ------------------------------------------------------------

valid_gr <- candidate_gr[as.character(seqnames(candidate_gr)) %in% fasta_seqnames]

cat("\nCandidates on available chromosomes:", length(valid_gr), "/", length(candidate_gr), "\n")


# ------------------------------------------------------------
# 17. Check chromosome bounds
# ------------------------------------------------------------

fasta_lengths <- width(fasta_index)
names(fasta_lengths) <- as.character(seqnames(fasta_index))

valid_chr_lengths <- fasta_lengths[as.character(seqnames(valid_gr))]

within_bounds <- start(valid_gr) >= 1 & end(valid_gr) <= valid_chr_lengths

cat("Regions within chromosome bounds:", sum(within_bounds), "/", length(within_bounds), "\n")

if (any(!within_bounds)) {
  cat("\nWARNING: out-of-bounds regions:\n")
  print(valid_gr[!within_bounds])
}

valid_gr <- valid_gr[within_bounds]


# ------------------------------------------------------------
# 18. Extract sequences
# ------------------------------------------------------------
# scanFa ignores strand: extract the reference sequence, then
# reverse-complement minus-strand TARs.

query_gr <- valid_gr
strand(query_gr) <- "*"

cat("\nExtracting sequences...\n")

seqs <- scanFa(fa, param = query_gr)
names(seqs) <- names(valid_gr)


# ------------------------------------------------------------
# 19. Orient by strand
# ------------------------------------------------------------

negative_idx <- which(as.character(strand(valid_gr)) == "-")

if (length(negative_idx) > 0) {
  seqs[negative_idx] <- reverseComplement(seqs[negative_idx])
}

cat("Reverse-complemented sequences:", length(negative_idx), "\n")


# ------------------------------------------------------------
# 20. Check lengths
# ------------------------------------------------------------

sequence_length <- width(seqs)
expected_length <- width(valid_gr)
length_match <- sequence_length == expected_length

cat("\n")
cat("====================================================\n")
cat("SEQUENCE CHECK\n")
cat("====================================================\n\n")

cat("Extracted sequences:", length(seqs), "\n")
cat("FASTA length == TAR length:", sum(length_match), "/", length(length_match), "\n")

if (any(!length_match)) {
  cat("Mismatched lengths:", sum(!length_match), "\n")
}


# ------------------------------------------------------------
# 21. Sequence metrics
# ------------------------------------------------------------

sequence_df <- data.frame(TAR_ID = names(seqs), sequence_length = width(seqs), stringsAsFactors = FALSE)

gc_freq <- letterFrequency(seqs, letters = c("G", "C"))

sequence_df$GC_content <- rowSums(gc_freq) / sequence_df$sequence_length * 100

n_count <- letterFrequency(seqs, letters = "N")[, 1]

sequence_df$N_count <- n_count
sequence_df$N_fraction <- n_count / sequence_df$sequence_length


# ------------------------------------------------------------
# 22. Join metrics
# ------------------------------------------------------------

candidates <- candidates %>%
  left_join(sequence_df, by = "TAR_ID")


# ------------------------------------------------------------
# 23. Sort by priority
# ------------------------------------------------------------

candidates <- candidates %>%
  mutate(
    priority_order = case_when(
      priority_group == "Priority_A" ~ 1,
      priority_group == "Priority_B" ~ 2,
      priority_group == "Priority_C" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  arrange(priority_order, desc(IPF_min_counts), IPF_CV_counts, desc(IPF_min_pct_cells))


# ------------------------------------------------------------
# 24. Subsets
# ------------------------------------------------------------

priority_A <- candidates %>% filter(priority_group == "Priority_A")
priority_B <- candidates %>% filter(priority_group == "Priority_B")
priority_C <- candidates %>% filter(priority_group == "Priority_C")
no_control <- candidates %>% filter(n_Control == 0)
high_confidence <- candidates %>% filter(high_confidence_pilot)


# ------------------------------------------------------------
# 25. Master table
# ------------------------------------------------------------

write.csv(
  candidates, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_candidates_for_coding_potential.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Per-priority tables
# ------------------------------------------------------------

write.csv(priority_A, file.path(OUTPUT_DIR, "Priority_A_candidates.csv"), row.names = FALSE)
write.csv(priority_B, file.path(OUTPUT_DIR, "Priority_B_candidates.csv"), row.names = FALSE)
write.csv(priority_C, file.path(OUTPUT_DIR, "Priority_C_candidates.csv"), row.names = FALSE)
write.csv(high_confidence, file.path(OUTPUT_DIR, "High_confidence_pilot_candidates.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 27. BED helper (0-based start)
# ------------------------------------------------------------

make_bed <- function(df) {
  df %>%
    transmute(chr = chr, bed_start = pmax(start - 1, 0), end = end, name = TAR_ID, score = 0, strand = strand)
}


# ------------------------------------------------------------
# 28. BED: all candidates
# ------------------------------------------------------------

write.table(
  make_bed(candidates),
  file = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_candidates.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


# ------------------------------------------------------------
# 29. BED: 0/1 control
# ------------------------------------------------------------

write.table(
  make_bed(no_control),
  file = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_0of1_control.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


# ------------------------------------------------------------
# 30. BED: Priority A
# ------------------------------------------------------------

write.table(
  make_bed(priority_A),
  file = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_Priority_A.bed"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


# ------------------------------------------------------------
# 31. FASTA: all candidates
# ------------------------------------------------------------

writeXStringSet(seqs, filepath = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_candidates_genomic.fa"))


# ------------------------------------------------------------
# 32. FASTA: 0/1 control
# ------------------------------------------------------------

ids_no_control <- no_control$TAR_ID

seqs_no_control <- seqs[names(seqs) %in% ids_no_control]

writeXStringSet(
  seqs_no_control,
  filepath = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_0of1_control_genomic.fa")
)


# ------------------------------------------------------------
# 33. FASTA: Priority A
# ------------------------------------------------------------

ids_A <- priority_A$TAR_ID

seqs_A <- seqs[names(seqs) %in% ids_A]

writeXStringSet(seqs_A, filepath = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_Priority_A_genomic.fa"))


# ------------------------------------------------------------
# 34. FASTA: high-confidence pilot
# ------------------------------------------------------------

ids_high_confidence <- high_confidence$TAR_ID

seqs_high_confidence <- seqs[names(seqs) %in% ids_high_confidence]

writeXStringSet(
  seqs_high_confidence,
  filepath = file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_high_confidence_pilot_genomic.fa")
)


# ------------------------------------------------------------
# 35. Sequence summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SEQUENCES\n")
cat("====================================================\n\n")

cat("Length:\n")
print(summary(candidates$sequence_length))

cat("\nGC (%):\n")
print(summary(candidates$GC_content))

cat("\nSequences with N > 0:", sum(candidates$N_count > 0, na.rm = TRUE), "\n")
cat("Sequences with >1% N:", sum(candidates$N_fraction > 0.01, na.rm = TRUE), "\n")


# ------------------------------------------------------------
# 36. Priority A summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("PRIORITY A\n")
cat("====================================================\n\n")

cat("Number Priority A:", nrow(priority_A), "\n")
cat("Number Priority B:", nrow(priority_B), "\n")
cat("Number Priority C:", nrow(priority_C), "\n")
cat("High-confidence pilot:", nrow(high_confidence), "\n\n")

cat("Min counts IPF - Priority A:\n")
print(summary(priority_A$IPF_min_counts))

cat("\nCV IPF - Priority A:\n")
print(summary(priority_A$IPF_CV_counts))

cat("\nMin % cells IPF - Priority A:\n")
print(summary(priority_A$IPF_min_pct_cells))


# ------------------------------------------------------------
# 37. Top 30 Priority A
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP 30 PRIORITY A\n")
cat("====================================================\n\n")

top_columns <- intersect(
  c(
    "TAR_ID", "chr", "start", "end", "strand", "length",
    "counts_VUILD53", "counts_VUILD63", "counts_VUILD64", "counts_VUHD67",
    "IPF_mean_counts", "IPF_min_counts", "IPF_CV_counts",
    "IPF_mean_pct_cells", "IPF_min_pct_cells",
    "sequence_length", "GC_content", "N_count",
    "priority_group", "high_confidence_pilot"
  ),
  colnames(priority_A)
)

print(head(priority_A[, top_columns, drop = FALSE], 30))


# ------------------------------------------------------------
# 38. Generated files
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GENERATED FILES\n")
cat("====================================================\n\n")

print(list.files(OUTPUT_DIR, full.names = FALSE))


# ------------------------------------------------------------
# 39. RDS
# ------------------------------------------------------------

saveRDS(candidates, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_candidates_for_coding_potential.rds"))


# ------------------------------------------------------------
# 40. Done
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CANDIDATE PRIORITIZATION FINISHED\n")
cat("====================================================\n\n")

cat("Output:\n", OUTPUT_DIR, "\n\n")

cat(
  "IMPORTANT:\n",
  "Sequences correspond to full genomic TAR intervals.\n",
  "They should not yet be interpreted as mature transcript sequences.\n"
)
