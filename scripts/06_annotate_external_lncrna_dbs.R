# ============================================================
# 06_annotate_external_lncrna_dbs.R
#
# lncRNA-IPF-sc | Step 6: external annotation of consensus uTARs
# Author: Jose A. Ovando-Ricardez
#
# Databases: LncBook, NONCODE, LNCipedia, FANTOM-CAT
#
# Compares the uTARs that are unannotated in GENCODE v47 with external
# lncRNA databases (strand-specific when the database carries strand).
#
# Note: "Unannotated_candidate" does NOT yet mean "novel lncRNA".
#
# Inputs:  data/TAR_catalog/annotation/IPF_uTAR_consensus_annotated.csv
#          data/annotations/external/{LncBook.gtf.gz, NONCODE.bed.gz,
#                                     LNCipedia.bed, FANTOM_lncRNA.bed}
# Outputs: data/TAR_catalog/external_annotation/,
#          data/plots/04_lncRNA/03_external_annotation/
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

INPUT_FILE   <- file.path(BASE_DIR, "data/TAR_catalog/annotation/IPF_uTAR_consensus_annotated.csv")
EXTERNAL_DIR <- file.path(BASE_DIR, "data/annotations/external")
OUTPUT_DIR   <- file.path(BASE_DIR, "data/TAR_catalog/external_annotation")
PLOT_DIR     <- file.path(BASE_DIR, "data/plots/04_lncRNA/03_external_annotation")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Check input
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(paste("Catalogue not found:", INPUT_FILE))
}

cat("\n")
cat("====================================================\n")
cat("EXTERNAL ANNOTATION OF CONSENSUS uTARs\n")
cat("====================================================\n\n")
cat("Input catalogue:\n", INPUT_FILE, "\n\n")


# ------------------------------------------------------------
# 3. Read catalogue
# ------------------------------------------------------------

consensus_df <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

cat("Total consensus uTARs:", nrow(consensus_df), "\n")


# ------------------------------------------------------------
# 4. Required columns
# ------------------------------------------------------------

required_columns <- c("TAR_ID", "chr", "start", "end", "strand", "annotation_class",
                      "n_samples", "n_IPF", "n_Control", "presence_class")

missing_columns <- setdiff(required_columns, colnames(consensus_df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}


# ------------------------------------------------------------
# 5. Length
# ------------------------------------------------------------

if (!"length" %in% colnames(consensus_df)) {
  consensus_df$length <- consensus_df$end - consensus_df$start + 1
}


# ------------------------------------------------------------
# 6. Candidates unannotated in GENCODE v47
#    ("Candidate_novel" from the previous step really means
#    "uTAR without GENCODE v47 annotation")
# ------------------------------------------------------------

consensus_df$GENCODE_unannotated <- consensus_df$annotation_class == "Candidate_novel"

cat("Without GENCODE v47 annotation:", sum(consensus_df$GENCODE_unannotated, na.rm = TRUE), "\n\n")


# ------------------------------------------------------------
# 7. Candidate GRanges
# ------------------------------------------------------------

candidate_df <- consensus_df[consensus_df$GENCODE_unannotated, , drop = FALSE]

candidate_gr <- GRanges(
  seqnames = candidate_df$chr,
  ranges   = IRanges(start = candidate_df$start, end = candidate_df$end),
  strand   = candidate_df$strand
)

mcols(candidate_gr)$TAR_ID <- candidate_df$TAR_ID

cat("Candidate GRanges:", length(candidate_gr), "\n\n")


# ------------------------------------------------------------
# 8. External database files
# ------------------------------------------------------------

database_files <- c(
  LncBook   = file.path(EXTERNAL_DIR, "LncBook.gtf.gz"),
  NONCODE   = file.path(EXTERNAL_DIR, "NONCODE.bed.gz"),
  LNCipedia = file.path(EXTERNAL_DIR, "LNCipedia.bed"),
  FANTOM    = file.path(EXTERNAL_DIR, "FANTOM_lncRNA.bed")
)

database_files[!file.exists(database_files)] <- NA_character_

cat("====================================================\n")
cat("EXTERNAL DATABASE FILES\n")
cat("====================================================\n\n")

print(data.frame(database = names(database_files), file = unname(database_files),
                 loaded = !is.na(database_files), row.names = NULL))


# ------------------------------------------------------------
# 9. Loading status
# ------------------------------------------------------------

database_status <- data.frame(
  database             = names(database_files),
  file                 = unname(database_files),
  loaded               = !is.na(database_files),
  n_features           = NA_integer_,
  n_candidate_overlaps = NA_integer_,
  n_unique_TAR         = NA_integer_,
  overlap_mode         = NA_character_,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10. Helper: feature identifier
# ------------------------------------------------------------

get_feature_id <- function(gr) {

  possible_columns <- c("gene_name", "gene_id", "transcript_name", "transcript_id",
                        "Name", "ID", "name")

  available <- intersect(possible_columns, colnames(mcols(gr)))

  if (length(available) == 0) {
    return(paste0("feature_", seq_along(gr)))
  }

  x <- as.character(mcols(gr)[[available[1]]])

  missing_id <- is.na(x) | x == ""
  x[missing_id] <- paste0("feature_", which(missing_id))

  return(x)
}


# ------------------------------------------------------------
# 11. Helper: strand information
# ------------------------------------------------------------

normalize_strand <- function(gr) {

  s <- as.character(strand(gr))
  informative <- s %in% c("+", "-")

  list(
    gr                     = gr,
    has_informative_strand = any(informative),
    informative_fraction   = mean(informative)
  )
}


# ------------------------------------------------------------
# 12. Initialise evidence columns
# ------------------------------------------------------------

for (db in names(database_files)) {
  consensus_df[[db]] <- NA
  consensus_df[[paste0(db, "_ID")]] <- NA_character_
}


# ------------------------------------------------------------
# 13. Import and compare each database
# ------------------------------------------------------------

external_overlap_list <- list()

for (db in names(database_files)) {

  cat("\n")
  cat("====================================================\n")
  cat("DATABASE:", db, "\n")
  cat("====================================================\n")

  db_file <- database_files[[db]]

  # Database not available
  if (is.na(db_file)) {
    cat("File NOT found.\n")
    cat("This database will be marked as NOT EVALUATED.\n")
    next
  }

  cat("Importing:\n", db_file, "\n")

  # Import
  db_gr <- rtracklayer::import(db_file)

  cat("Features imported:", length(db_gr), "\n")

  database_status$n_features[database_status$database == db] <- length(db_gr)

  # Remove invalid coordinates
  valid <- !is.na(start(db_gr)) & !is.na(end(db_gr)) & start(db_gr) >= 1 & end(db_gr) >= start(db_gr)

  db_gr <- db_gr[valid]

  cat("Features with valid coordinates:", length(db_gr), "\n")

  # Identifier
  db_ids <- get_feature_id(db_gr)
  mcols(db_gr)$external_id <- db_ids

  # Shared sequence names
  common_seqlevels <- intersect(GenomeInfoDb::seqlevels(candidate_gr), GenomeInfoDb::seqlevels(db_gr))

  cat("Shared seqlevels:", length(common_seqlevels), "\n")

  if (length(common_seqlevels) == 0) {
    warning(paste(db, "shares no chromosome names with the catalogue.",
                  "Check chr1 vs 1 and the GRCh38 assembly."))
    next
  }

  # Keep shared chromosomes only
  candidate_db_gr <- GenomeInfoDb::keepSeqlevels(candidate_gr, common_seqlevels, pruning.mode = "coarse")
  db_gr           <- GenomeInfoDb::keepSeqlevels(db_gr, common_seqlevels, pruning.mode = "coarse")

  cat("Candidates after seqlevel filtering:", length(candidate_db_gr), "\n")
  cat("External features after seqlevel filtering:", length(db_gr), "\n")

  # Strand
  strand_info <- normalize_strand(db_gr)

  cat("Fraction of features with +/- strand:", round(strand_info$informative_fraction, 4), "\n")

  # Overlap
  if (strand_info$informative_fraction > 0.5) {

    cat("Strand-specific overlap.\n")

    hits <- GenomicRanges::findOverlaps(candidate_db_gr, db_gr, ignore.strand = FALSE)
    overlap_mode <- "strand_specific"

  } else {

    cat("Strand missing/insufficient.\n")
    cat("Overlap ignoring strand.\n")

    hits <- GenomicRanges::findOverlaps(candidate_db_gr, db_gr, ignore.strand = TRUE)
    overlap_mode <- "strand_ignored"
  }

  database_status$overlap_mode[database_status$database == db] <- overlap_mode

  cat("Overlaps found:", length(hits), "\n")

  # No hits
  if (length(hits) == 0) {
    database_status$n_candidate_overlaps[database_status$database == db] <- 0
    database_status$n_unique_TAR[database_status$database == db] <- 0
    consensus_df[[db]][consensus_df$GENCODE_unannotated] <- FALSE
    next
  }

  # Evidence table
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)

  tmp <- data.frame(
    TAR_ID          = mcols(candidate_db_gr)$TAR_ID[q],
    database        = db,
    external_id     = mcols(db_gr)$external_id[s],
    external_chr    = as.character(seqnames(db_gr)[s]),
    external_start  = start(db_gr)[s],
    external_end    = end(db_gr)[s],
    external_strand = as.character(strand(db_gr)[s]),
    overlap_mode    = overlap_mode,
    stringsAsFactors = FALSE
  )

  tmp <- unique(tmp)

  external_overlap_list[[db]] <- tmp

  database_status$n_candidate_overlaps[database_status$database == db] <- nrow(tmp)

  unique_tar <- unique(tmp$TAR_ID)

  database_status$n_unique_TAR[database_status$database == db] <- length(unique_tar)

  cat("Unique uTARs with evidence:", length(unique_tar), "\n")

  # TRUE/FALSE flags
  candidate_ids <- consensus_df$TAR_ID[consensus_df$GENCODE_unannotated]

  consensus_df[[db]][consensus_df$GENCODE_unannotated] <- candidate_ids %in% unique_tar

  # External IDs
  id_by_tar <- tapply(tmp$external_id, tmp$TAR_ID, function(x) {
    paste(unique(x), collapse = ";")
  })

  id_column <- paste0(db, "_ID")

  matching_rows <- match(names(id_by_tar), consensus_df$TAR_ID)

  consensus_df[matching_rows, id_column] <- unname(id_by_tar)
}


# ------------------------------------------------------------
# 14. Combine overlap tables
# ------------------------------------------------------------

if (length(external_overlap_list) > 0) {

  external_overlap_df <- do.call(rbind, external_overlap_list)
  rownames(external_overlap_df) <- NULL

} else {

  external_overlap_df <- data.frame(
    TAR_ID          = character(),
    database        = character(),
    external_id     = character(),
    external_chr    = character(),
    external_start  = integer(),
    external_end    = integer(),
    external_strand = character(),
    overlap_mode    = character(),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 15. Number of external sources
# ------------------------------------------------------------

external_columns <- c("LncBook", "NONCODE", "LNCipedia", "FANTOM")

consensus_df$n_external_lncRNA_sources <- rowSums(
  sapply(external_columns, function(x) {
    consensus_df[[x]] %in% TRUE
  }),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 16. Final classification
# ------------------------------------------------------------

consensus_df$annotation_class_external <- NA_character_

# Already a GENCODE lncRNA
consensus_df$annotation_class_external[
  consensus_df$annotation_class == "GENCODE_lncRNA"
] <- "GENCODE_lncRNA"

# Other GENCODE genes
consensus_df$annotation_class_external[
  consensus_df$annotation_class == "GENCODE_other"
] <- "GENCODE_other"

# External lncRNA evidence
consensus_df$annotation_class_external[
  consensus_df$GENCODE_unannotated & consensus_df$n_external_lncRNA_sources > 0
] <- "External_lncRNA"

# No external evidence
consensus_df$annotation_class_external[
  consensus_df$GENCODE_unannotated & consensus_df$n_external_lncRNA_sources == 0
] <- "Unannotated_candidate"


# ------------------------------------------------------------
# 17. Prioritisation flags
# ------------------------------------------------------------

consensus_df$length_ge200 <- consensus_df$length >= 200

consensus_df$unannotated_ge200 <-
  consensus_df$annotation_class_external == "Unannotated_candidate" &
  consensus_df$length >= 200

consensus_df$unannotated_IPF_recurrent <- consensus_df$unannotated_ge200 & consensus_df$n_IPF >= 2

consensus_df$unannotated_all3_IPF <- consensus_df$unannotated_ge200 & consensus_df$n_IPF == 3

consensus_df$unannotated_all3_IPF_no_control <-
  consensus_df$unannotated_ge200 &
  consensus_df$n_IPF == 3 &
  consensus_df$n_Control == 0


# ------------------------------------------------------------
# 18. Subsets
# ------------------------------------------------------------

external_lncRNA <- consensus_df[consensus_df$annotation_class_external == "External_lncRNA", , drop = FALSE]

unannotated_candidates <- consensus_df[consensus_df$annotation_class_external == "Unannotated_candidate", , drop = FALSE]

unannotated_ge200 <- consensus_df[consensus_df$unannotated_ge200, , drop = FALSE]

unannotated_recurrent <- consensus_df[consensus_df$unannotated_IPF_recurrent, , drop = FALSE]

unannotated_all3_IPF <- consensus_df[consensus_df$unannotated_all3_IPF, , drop = FALSE]

unannotated_all3_IPF_no_control <- consensus_df[consensus_df$unannotated_all3_IPF_no_control, , drop = FALSE]


# ------------------------------------------------------------
# 19. Save tables
# ------------------------------------------------------------

write.csv(consensus_df, file.path(OUTPUT_DIR, "IPF_uTAR_consensus_external_annotated.csv"), row.names = FALSE)
write.csv(external_overlap_df, file.path(OUTPUT_DIR, "external_lncRNA_overlap_details.csv"), row.names = FALSE)
write.csv(database_status, file.path(OUTPUT_DIR, "external_database_status.csv"), row.names = FALSE)
write.csv(external_lncRNA, file.path(OUTPUT_DIR, "IPF_uTAR_external_lncRNA.csv"), row.names = FALSE)
write.csv(unannotated_candidates, file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_candidates.csv"), row.names = FALSE)
write.csv(unannotated_ge200, file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_candidates_ge200bp.csv"), row.names = FALSE)
write.csv(unannotated_recurrent, file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_recurrent_ge2IPF.csv"), row.names = FALSE)
write.csv(unannotated_all3_IPF, file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_all3IPF.csv"), row.names = FALSE)
write.csv(unannotated_all3_IPF_no_control, file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_all3IPF_noControl.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 20. BED of candidates >= 200 bp (0-based start)
# ------------------------------------------------------------

if (nrow(unannotated_ge200) > 0) {

  bed_df <- data.frame(
    chr    = unannotated_ge200$chr,
    start  = pmax(unannotated_ge200$start - 1, 0),
    end    = unannotated_ge200$end,
    name   = unannotated_ge200$TAR_ID,
    score  = 0,
    strand = unannotated_ge200$strand
  )

  write.table(bed_df, file = file.path(OUTPUT_DIR, "IPF_uTAR_unannotated_candidates_ge200bp.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# ------------------------------------------------------------
# 21. Console summary
# ------------------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("MAIN RESULTS\n")
cat("====================================================\n\n")

cat("1. EXTERNAL DATABASE STATUS\n")
cat("----------------------------------------------------\n")
print(database_status)

cat("\n")
cat("2. FINAL CLASSIFICATION\n")
cat("----------------------------------------------------\n")
print(table(consensus_df$annotation_class_external))

cat("\n")
cat("3. EVIDENCE PER EXTERNAL DATABASE\n")
cat("----------------------------------------------------\n")

for (db in external_columns) {
  cat(db, ":", sum(consensus_df[[db]] %in% TRUE, na.rm = TRUE), "\n")
}

cat("\n")
cat("4. NUMBER OF EXTERNAL SOURCES PER uTAR\n")
cat("----------------------------------------------------\n")
print(table(consensus_df$n_external_lncRNA_sources[consensus_df$GENCODE_unannotated]))

cat("\n")
cat("5. UNANNOTATED CANDIDATES\n")
cat("----------------------------------------------------\n")
cat("No external annotation:",
    sum(consensus_df$annotation_class_external == "Unannotated_candidate", na.rm = TRUE), "\n")
cat("Unannotated >=200 bp:",
    sum(consensus_df$unannotated_ge200, na.rm = TRUE), "\n")
cat("Unannotated >=200 bp and present in >=2 IPF:",
    sum(consensus_df$unannotated_IPF_recurrent, na.rm = TRUE), "\n")
cat("Unannotated >=200 bp and present in 3/3 IPF:",
    sum(consensus_df$unannotated_all3_IPF, na.rm = TRUE), "\n")
cat("Unannotated >=200 bp, present in 3/3 IPF and 0/1 control:",
    sum(consensus_df$unannotated_all3_IPF_no_control, na.rm = TRUE), "\n")

cat("\n")
cat("6. IPF RECURRENCE OF CANDIDATES >=200 bp\n")
cat("----------------------------------------------------\n")
print(table(consensus_df$n_IPF[consensus_df$unannotated_ge200]))

cat("\n")
cat("7. IPF / CONTROL PRESENCE\n")
cat("----------------------------------------------------\n")
print(table(consensus_df$presence_class[consensus_df$unannotated_ge200]))


# ------------------------------------------------------------
# 22. Figure 1: final classification
# ------------------------------------------------------------

classification_table <- table(consensus_df$annotation_class_external)

pdf(file.path(PLOT_DIR, "01_external_annotation_classification.pdf"), width = 8, height = 6)
barplot(classification_table, las = 2, ylab = "Number of uTARs",
        main = "Classification after external annotation")
dev.off()

png(file.path(PLOT_DIR, "01_external_annotation_classification.png"), width = 1800, height = 1400, res = 200)
barplot(classification_table, las = 2, ylab = "Number of uTARs",
        main = "Classification after external annotation")
dev.off()


# ------------------------------------------------------------
# 23. Figure 2: evidence per database
# ------------------------------------------------------------

db_counts <- sapply(external_columns, function(db) {
  sum(consensus_df[[db]] %in% TRUE, na.rm = TRUE)
})

pdf(file.path(PLOT_DIR, "02_external_database_hits.pdf"), width = 8, height = 6)
barplot(db_counts, ylab = "Unique uTARs", main = "uTARs with evidence in external databases")
dev.off()

png(file.path(PLOT_DIR, "02_external_database_hits.png"), width = 1800, height = 1400, res = 200)
barplot(db_counts, ylab = "Unique uTARs", main = "uTARs with evidence in external databases")
dev.off()


# ------------------------------------------------------------
# 24. Figure 3: number of sources
# ------------------------------------------------------------

source_table <- table(consensus_df$n_external_lncRNA_sources[consensus_df$GENCODE_unannotated])

pdf(file.path(PLOT_DIR, "03_number_external_sources.pdf"), width = 8, height = 6)
barplot(source_table, xlab = "Number of external databases", ylab = "Number of uTARs",
        main = "Number of external sources per uTAR")
dev.off()

png(file.path(PLOT_DIR, "03_number_external_sources.png"), width = 1800, height = 1400, res = 200)
barplot(source_table, xlab = "Number of external databases", ylab = "Number of uTARs",
        main = "Number of external sources per uTAR")
dev.off()


# ------------------------------------------------------------
# 25. Figure 4: IPF recurrence
# ------------------------------------------------------------

recurrence_table <- table(consensus_df$n_IPF[consensus_df$unannotated_ge200])

pdf(file.path(PLOT_DIR, "04_unannotated_recurrence_IPF.pdf"), width = 8, height = 6)
barplot(recurrence_table, xlab = "Number of IPF samples", ylab = "Number of candidates",
        main = "Recurrence of unannotated candidates >=200 bp")
dev.off()

png(file.path(PLOT_DIR, "04_unannotated_recurrence_IPF.png"), width = 1800, height = 1400, res = 200)
barplot(recurrence_table, xlab = "Number of IPF samples", ylab = "Number of candidates",
        main = "Recurrence of unannotated candidates >=200 bp")
dev.off()


# ------------------------------------------------------------
# 26. Figure 5: IPF / control presence
# ------------------------------------------------------------

presence_table <- table(consensus_df$presence_class[consensus_df$unannotated_ge200])

pdf(file.path(PLOT_DIR, "05_unannotated_presence_IPF_control.pdf"), width = 9, height = 6)
barplot(presence_table, las = 2, ylab = "Number of candidates",
        main = "Presence of unannotated candidates >=200 bp")
dev.off()

png(file.path(PLOT_DIR, "05_unannotated_presence_IPF_control.png"), width = 1800, height = 1400, res = 200)
barplot(presence_table, las = 2, ylab = "Number of candidates",
        main = "Presence of unannotated candidates >=200 bp")
dev.off()


# ------------------------------------------------------------
# 27. Top candidates: 3/3 IPF and 0/1 control
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP: 3/3 IPF AND ABSENT IN CONTROL\n")
cat("====================================================\n\n")

top_columns <- intersect(
  c("TAR_ID", "chr", "start", "end", "strand", "length", "n_samples", "n_IPF",
    "n_Control", "presence_class", "n_external_lncRNA_sources"),
  colnames(unannotated_all3_IPF_no_control)
)

print(head(unannotated_all3_IPF_no_control[, top_columns, drop = FALSE], 30))


# ------------------------------------------------------------
# 28. RDS
# ------------------------------------------------------------

saveRDS(consensus_df, file.path(OUTPUT_DIR, "IPF_uTAR_consensus_external_annotated.rds"))


# ------------------------------------------------------------
# 29. Done
# ------------------------------------------------------------

cat("\n\n")
cat("====================================================\n")
cat("EXTERNAL ANNOTATION FINISHED\n")
cat("====================================================\n\n")

cat("Results:\n", OUTPUT_DIR, "\n\n")
cat("Figures:\n", PLOT_DIR, "\n\n")

cat(
  "IMPORTANT:\n",
  "Unannotated_candidate = not annotated in GENCODE v47 or in the\n",
  "external databases evaluated. It does NOT yet imply a novel lncRNA.\n",
  "Next steps: coding potential, structure,\n",
  "expression and cell specificity.\n"
)
