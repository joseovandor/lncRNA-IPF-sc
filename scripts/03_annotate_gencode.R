# ============================================================
# 03_annotate_gencode.R
#
# lncRNA-IPF-sc | Step 3: annotate consensus uTARs with GENCODE v47
# Author: Jose A. Ovando-Ricardez
#
# Input: data/TAR_catalog/IPF_uTAR_consensus_catalog.csv
#
# 1. Imports the consensus uTAR catalogue.
# 2. Imports GENCODE v47 (GTF/GFF3) from a local file.
# 3. Finds overlaps with GENCODE lncRNA genes (strand-specific).
# 4. Finds overlaps with other GENCODE genes (strand-specific).
# 5. Optionally looks for evidence in LncBook, NONCODE, LNCipedia and
#    FANTOM BED files if present (the full external annotation is done
#    in 06_annotate_external_lncrna_dbs.R).
# 6. Classifies each IPFTAR as GENCODE_lncRNA, External_lncRNA,
#    GENCODE_other or Candidate_novel.
# 7. Keeps the IPF/control recurrence columns from the catalogue.
# 8. Writes annotation tables and figures (PNG + PDF).
#
# Notes:
# - "Candidate_novel" only means that no evidence was found in the
#   loaded annotations; it is NOT a confirmed novel lncRNA.
# - Length, coding potential (CPAT/RNAsamba), TSS/PAS, splicing,
#   conservation and reproducibility must be evaluated afterwards.
# - No candidate is removed here.
#
# Outputs: data/TAR_catalog/annotation/, data/plots/04_lncRNA/01_annotation/
# ============================================================

source("config/config.R")

library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(ggplot2)
library(Matrix)


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

data_dir <- DATA_DIR

catalog_dir           <- file.path(data_dir, "TAR_catalog")
annotation_dir        <- file.path(data_dir, "annotations")
annotation_output_dir <- file.path(catalog_dir, "annotation")
plots_dir             <- file.path(data_dir, "plots")
annotation_plot_dir   <- file.path(plots_dir, "04_lncRNA", "01_annotation")

dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n====================================================\n")
cat("DIRECTORIES\n")
cat("====================================================\n")
cat("Input annotations:\n", annotation_dir, "\n\n")
cat("Annotation results:\n", annotation_output_dir, "\n\n")
cat("Figures:\n", annotation_plot_dir, "\n")


# ------------------------------------------------------------
# 2. Consensus catalogue file
# ------------------------------------------------------------

consensus_file <- file.path(catalog_dir, "IPF_uTAR_consensus_catalog.csv")

if (!file.exists(consensus_file)) {
  stop(paste0("\nFile not found:\n", consensus_file,
              "\n\nRun the consensus catalogue script first."))
}


# ------------------------------------------------------------
# 3. Annotation files
#    GENCODE v47: data/annotations/gencode.v47.annotation.gtf.gz
#    Optional external BEDs (skipped if absent):
#      LncBook.bed(.gz), NONCODE.bed(.gz), LNCipedia.bed(.gz),
#      FANTOM_lncRNA.bed(.gz)
# ------------------------------------------------------------

gencode_file <- file.path(annotation_dir, "gencode.v47.annotation.gtf.gz")

external_files <- list(
  LncBook   = c(file.path(annotation_dir, "LncBook.bed"), file.path(annotation_dir, "LncBook.bed.gz")),
  NONCODE   = c(file.path(annotation_dir, "NONCODE.bed"), file.path(annotation_dir, "NONCODE.bed.gz")),
  LNCipedia = c(file.path(annotation_dir, "LNCipedia.bed"), file.path(annotation_dir, "LNCipedia.bed.gz")),
  FANTOM    = c(file.path(annotation_dir, "FANTOM_lncRNA.bed"), file.path(annotation_dir, "FANTOM_lncRNA.bed.gz"))
)


# ------------------------------------------------------------
# 4. Check GENCODE
# ------------------------------------------------------------

if (!file.exists(gencode_file)) {
  stop(paste0("\nGENCODE v47 file not found:\n\n", gencode_file, "\n\n",
              "Place the GENCODE v47 GTF in the annotations directory ",
              "before continuing.\n"))
}


# ------------------------------------------------------------
# 5. Read consensus catalogue
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("READING CONSENSUS CATALOGUE\n")
cat("====================================================\n")

consensus_df <- read.csv(consensus_file, stringsAsFactors = FALSE)

cat("\nConsensus regions:", nrow(consensus_df), "\n")

required_columns <- c("TAR_ID", "chr", "start", "end", "strand", "length")
missing_columns  <- setdiff(required_columns, colnames(consensus_df))

if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}


# ------------------------------------------------------------
# 6. Catalogue to GRanges
# ------------------------------------------------------------

consensus_gr <- GRanges(
  seqnames = consensus_df$chr,
  ranges   = IRanges(start = consensus_df$start, end = consensus_df$end),
  strand   = consensus_df$strand
)

mcols(consensus_gr)$TAR_ID <- consensus_df$TAR_ID


# ------------------------------------------------------------
# 7. Import GENCODE v47
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("IMPORTING GENCODE v47\n")
cat("====================================================\n")

gencode <- import(gencode_file)

cat("\nGENCODE features imported:", length(gencode), "\n")


# ------------------------------------------------------------
# 8. Keep gene-level features (avoids counting several exons or
#    transcripts of the same gene as separate evidence)
# ------------------------------------------------------------

if ("type" %in% colnames(mcols(gencode))) {
  gencode_genes <- gencode[mcols(gencode)$type == "gene"]
} else {
  gencode_genes <- gencode
}

cat("GENCODE genes:", length(gencode_genes), "\n")


# ------------------------------------------------------------
# 9. Gene type / name / id
# ------------------------------------------------------------

gencode_meta <- as.data.frame(mcols(gencode_genes))

if ("gene_type" %in% colnames(gencode_meta)) {
  gencode_gene_type <- as.character(gencode_meta$gene_type)
} else if ("gene_biotype" %in% colnames(gencode_meta)) {
  gencode_gene_type <- as.character(gencode_meta$gene_biotype)
} else {
  gencode_gene_type <- rep(NA_character_, length(gencode_genes))
}

if ("gene_name" %in% colnames(gencode_meta)) {
  gencode_gene_name <- as.character(gencode_meta$gene_name)
} else {
  gencode_gene_name <- rep(NA_character_, length(gencode_genes))
}

if ("gene_id" %in% colnames(gencode_meta)) {
  gencode_gene_id <- as.character(gencode_meta$gene_id)
} else {
  gencode_gene_id <- rep(NA_character_, length(gencode_genes))
}

mcols(gencode_genes)$gene_type_clean <- gencode_gene_type
mcols(gencode_genes)$gene_name_clean <- gencode_gene_name
mcols(gencode_genes)$gene_id_clean   <- gencode_gene_id


# ------------------------------------------------------------
# 10. lncRNA-compatible biotypes (flexible match because releases
#     use slightly different biotype names)
# ------------------------------------------------------------

lncrna_pattern <- paste(
  c("lncRNA", "lincRNA", "antisense", "processed_transcript", "sense_intronic",
    "sense_overlapping", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA",
    "macro_lncRNA", "non_coding"),
  collapse = "|"
)

is_gencode_lncrna <- grepl(lncrna_pattern, gencode_gene_type, ignore.case = TRUE)
is_gencode_lncrna[is.na(is_gencode_lncrna)] <- FALSE

gencode_lncrna <- gencode_genes[is_gencode_lncrna]
gencode_other  <- gencode_genes[!is_gencode_lncrna]

cat("\nGenes classified as lncRNA/noncoding:", length(gencode_lncrna), "\n")
cat("Other GENCODE genes:", length(gencode_other), "\n")


# ------------------------------------------------------------
# 11. Check shared seqlevels (chr/contig names are not modified
#     automatically because samples include alternative contigs)
# ------------------------------------------------------------

common_seq_gencode <- intersect(seqlevels(consensus_gr), seqlevels(gencode_genes))

cat("\nSeqlevels shared by catalogue and GENCODE:", length(common_seq_gencode), "\n")

if (length(common_seq_gencode) == 0) {
  stop(paste0("\nNo chromosomes/contigs are shared between the catalogue ",
              "and GENCODE.\nCheck whether one uses 'chr1' and the other '1'."))
}


# ------------------------------------------------------------
# 12. Overlap with GENCODE lncRNA (strand-specific, consistent with
#     the TAR logic)
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("OVERLAP WITH GENCODE lncRNA\n")
cat("====================================================\n")

hits_lnc <- findOverlaps(consensus_gr, gencode_lncrna, ignore.strand = FALSE)

cat("\nLocus-GENCODE lncRNA overlaps:", length(hits_lnc), "\n")

gencode_lnc_hits_df <- data.frame(
  TAR_ID    = mcols(consensus_gr)$TAR_ID[queryHits(hits_lnc)],
  gene_id   = mcols(gencode_lncrna)$gene_id_clean[subjectHits(hits_lnc)],
  gene_name = mcols(gencode_lncrna)$gene_name_clean[subjectHits(hits_lnc)],
  gene_type = mcols(gencode_lncrna)$gene_type_clean[subjectHits(hits_lnc)],
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 13. Overlap with other GENCODE genes
# ------------------------------------------------------------

hits_other <- findOverlaps(consensus_gr, gencode_other, ignore.strand = FALSE)

cat("Overlaps with other GENCODE genes:", length(hits_other), "\n")

gencode_other_hits_df <- data.frame(
  TAR_ID    = mcols(consensus_gr)$TAR_ID[queryHits(hits_other)],
  gene_id   = mcols(gencode_other)$gene_id_clean[subjectHits(hits_other)],
  gene_name = mcols(gencode_other)$gene_name_clean[subjectHits(hits_other)],
  gene_type = mcols(gencode_other)$gene_type_clean[subjectHits(hits_other)],
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14. Summarise GENCODE annotation per TAR
# ------------------------------------------------------------

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = ";")
}

if (nrow(gencode_lnc_hits_df) > 0) {
  lnc_ids   <- split(gencode_lnc_hits_df$gene_id, gencode_lnc_hits_df$TAR_ID)
  lnc_names <- split(gencode_lnc_hits_df$gene_name, gencode_lnc_hits_df$TAR_ID)
  lnc_types <- split(gencode_lnc_hits_df$gene_type, gencode_lnc_hits_df$TAR_ID)
} else {
  lnc_ids   <- list()
  lnc_names <- list()
  lnc_types <- list()
}

if (nrow(gencode_other_hits_df) > 0) {
  other_ids   <- split(gencode_other_hits_df$gene_id, gencode_other_hits_df$TAR_ID)
  other_names <- split(gencode_other_hits_df$gene_name, gencode_other_hits_df$TAR_ID)
  other_types <- split(gencode_other_hits_df$gene_type, gencode_other_hits_df$TAR_ID)
} else {
  other_ids   <- list()
  other_names <- list()
  other_types <- list()
}

consensus_df$GENCODE_lncRNA <- consensus_df$TAR_ID %in% unique(gencode_lnc_hits_df$TAR_ID)

consensus_df$GENCODE_lncRNA_gene_id <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(lnc_ids))) {
      return(NA_character_)
    }
    collapse_unique(lnc_ids[[id]])
  },
  character(1)
)

consensus_df$GENCODE_lncRNA_gene_name <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(lnc_names))) {
      return(NA_character_)
    }
    collapse_unique(lnc_names[[id]])
  },
  character(1)
)

consensus_df$GENCODE_lncRNA_gene_type <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(lnc_types))) {
      return(NA_character_)
    }
    collapse_unique(lnc_types[[id]])
  },
  character(1)
)

consensus_df$GENCODE_other <- consensus_df$TAR_ID %in% unique(gencode_other_hits_df$TAR_ID)

consensus_df$GENCODE_other_gene_id <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(other_ids))) {
      return(NA_character_)
    }
    collapse_unique(other_ids[[id]])
  },
  character(1)
)

consensus_df$GENCODE_other_gene_name <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(other_names))) {
      return(NA_character_)
    }
    collapse_unique(other_names[[id]])
  },
  character(1)
)

consensus_df$GENCODE_other_gene_type <- vapply(
  consensus_df$TAR_ID,
  function(id) {
    if (!(id %in% names(other_types))) {
      return(NA_character_)
    }
    collapse_unique(other_types[[id]])
  },
  character(1)
)


# ------------------------------------------------------------
# 15. Helper: first existing file
# ------------------------------------------------------------

find_existing_file <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}


# ------------------------------------------------------------
# 16. Helper: import an external BED and overlap it with the catalogue
#     (strand-specific when the BED has informative strand)
# ------------------------------------------------------------

import_external_bed <- function(file, db_name) {

  cat("\nImporting", db_name, "from:\n", file, "\n")

  gr <- import(file)

  if (length(gr) == 0) {
    cat("  File has no regions.\n")
    return(NULL)
  }

  strand_values <- unique(as.character(strand(gr)))
  strand_values <- strand_values[!is.na(strand_values)]

  has_informative_strand <- any(strand_values %in% c("+", "-"))

  hits <- findOverlaps(consensus_gr, gr, ignore.strand = !has_informative_strand)

  cat("  Overlapping regions:", length(unique(queryHits(hits))), "\n")

  evidence <- data.frame(
    TAR_ID   = mcols(consensus_gr)$TAR_ID[queryHits(hits)],
    database = db_name,
    stringsAsFactors = FALSE
  )

  # Try to recover a name/ID from the BED
  gr_meta <- as.data.frame(mcols(gr))

  candidate_columns <- c("name", "Name", "gene_name", "gene_id", "ID", "id", "transcript_id")
  name_column <- candidate_columns[candidate_columns %in% colnames(gr_meta)]

  if (length(name_column) > 0) {
    evidence$external_id <- as.character(gr_meta[subjectHits(hits), name_column[1]])
  } else {
    evidence$external_id <- NA_character_
  }

  evidence
}


# ------------------------------------------------------------
# 17. Import available external databases
# ------------------------------------------------------------

external_evidence_list <- list()

for (db in names(external_files)) {

  db_file <- find_existing_file(external_files[[db]])

  if (is.na(db_file)) {
    cat("\n", db, ": file not found; skipped.\n", sep = "")
    next
  }

  evidence <- import_external_bed(db_file, db)

  if (!is.null(evidence)) {
    external_evidence_list[[db]] <- evidence
  }
}


# ------------------------------------------------------------
# 18. Combine external evidence
# ------------------------------------------------------------

if (length(external_evidence_list) > 0) {
  external_evidence <- do.call(rbind, external_evidence_list)
  rownames(external_evidence) <- NULL
} else {
  external_evidence <- data.frame(
    TAR_ID      = character(0),
    database    = character(0),
    external_id = character(0),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 19. One flag column per external database
# ------------------------------------------------------------

external_db_names <- names(external_files)

for (db in external_db_names) {
  consensus_df[[db]] <- consensus_df$TAR_ID %in%
    unique(external_evidence$TAR_ID[external_evidence$database == db])
}


# ------------------------------------------------------------
# 20. External IDs per database
# ------------------------------------------------------------

for (db in external_db_names) {

  db_evidence <- external_evidence[external_evidence$database == db, , drop = FALSE]

  if (nrow(db_evidence) == 0) {
    consensus_df[[paste0(db, "_ID")]] <- NA_character_
    next
  }

  db_ids <- split(db_evidence$external_id, db_evidence$TAR_ID)

  consensus_df[[paste0(db, "_ID")]] <- vapply(
    consensus_df$TAR_ID,
    function(id) {
      if (!(id %in% names(db_ids))) {
        return(NA_character_)
      }
      collapse_unique(db_ids[[id]])
    },
    character(1)
  )
}


# ------------------------------------------------------------
# 21. Number of external databases supporting each TAR
# ------------------------------------------------------------

external_flag_matrix <- as.matrix(consensus_df[, external_db_names, drop = FALSE])
storage.mode(external_flag_matrix) <- "numeric"

consensus_df$n_external_databases <- rowSums(external_flag_matrix, na.rm = TRUE)
consensus_df$any_external_lncRNA  <- consensus_df$n_external_databases > 0


# ------------------------------------------------------------
# 22. Final annotation class
#     Priority: GENCODE_lncRNA > External_lncRNA > GENCODE_other > Candidate_novel
# ------------------------------------------------------------

consensus_df$annotation_class <- "Candidate_novel"
consensus_df$annotation_class[consensus_df$GENCODE_other] <- "GENCODE_other"
consensus_df$annotation_class[consensus_df$any_external_lncRNA] <- "External_lncRNA"
consensus_df$annotation_class[consensus_df$GENCODE_lncRNA] <- "GENCODE_lncRNA"


# ------------------------------------------------------------
# 23. Total evidence
# ------------------------------------------------------------

consensus_df$n_lncRNA_sources <- as.integer(consensus_df$GENCODE_lncRNA) +
  consensus_df$n_external_databases


# ------------------------------------------------------------
# 24. Candidate flags
# ------------------------------------------------------------

consensus_df$candidate_novel       <- consensus_df$annotation_class == "Candidate_novel"
consensus_df$candidate_novel_ge200 <- consensus_df$candidate_novel & consensus_df$length >= 200


# ------------------------------------------------------------
# 25. Candidates recurrent in IPF (>=2 of 3 IPF samples).
#     Prioritisation only; no statistical association with IPF.
# ------------------------------------------------------------

if ("n_IPF" %in% colnames(consensus_df)) {
  consensus_df$candidate_IPF_recurrent <- consensus_df$candidate_novel_ge200 &
    consensus_df$n_IPF >= 2
} else {
  consensus_df$candidate_IPF_recurrent <- FALSE
}


# ------------------------------------------------------------
# 26. Save full annotated table
# ------------------------------------------------------------

annotated_file <- file.path(annotation_output_dir, "IPF_uTAR_consensus_annotated.csv")

write.csv(consensus_df, annotated_file, row.names = FALSE)


# ------------------------------------------------------------
# 27. Save subset tables and overlap details
# ------------------------------------------------------------

write.csv(consensus_df[consensus_df$annotation_class == "GENCODE_lncRNA", , drop = FALSE],
          file.path(annotation_output_dir, "IPF_uTAR_GENCODE_lncRNA.csv"), row.names = FALSE)

write.csv(consensus_df[consensus_df$annotation_class == "External_lncRNA", , drop = FALSE],
          file.path(annotation_output_dir, "IPF_uTAR_external_lncRNA.csv"), row.names = FALSE)

write.csv(consensus_df[consensus_df$candidate_novel, , drop = FALSE],
          file.path(annotation_output_dir, "IPF_uTAR_candidate_novel.csv"), row.names = FALSE)

write.csv(consensus_df[consensus_df$candidate_novel_ge200, , drop = FALSE],
          file.path(annotation_output_dir, "IPF_uTAR_candidate_novel_ge200bp.csv"), row.names = FALSE)

write.csv(consensus_df[consensus_df$candidate_IPF_recurrent, , drop = FALSE],
          file.path(annotation_output_dir, "IPF_uTAR_candidate_novel_IPF_recurrent.csv"), row.names = FALSE)

# Detailed overlap tables
write.csv(gencode_lnc_hits_df,
          file.path(annotation_output_dir, "GENCODE_lncRNA_overlap_details.csv"), row.names = FALSE)

write.csv(gencode_other_hits_df,
          file.path(annotation_output_dir, "GENCODE_other_overlap_details.csv"), row.names = FALSE)

write.csv(external_evidence,
          file.path(annotation_output_dir, "external_lncRNA_overlap_details.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 28. BED of unannotated candidates >= 200 bp
# ------------------------------------------------------------

candidate_df <- consensus_df[consensus_df$candidate_novel_ge200, , drop = FALSE]

candidate_bed <- data.frame(
  chr    = candidate_df$chr,
  start  = candidate_df$start - 1,
  end    = candidate_df$end,
  name   = candidate_df$TAR_ID,
  score  = 0,
  strand = candidate_df$strand
)

write.table(candidate_bed,
            file = file.path(annotation_output_dir, "IPF_uTAR_candidate_novel_ge200bp.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# ------------------------------------------------------------
# 29. Helper: save PNG + PDF
# ------------------------------------------------------------

save_annotation_plot <- function(plot_object, filename, width = 8, height = 6) {
  ggsave(filename = file.path(annotation_plot_dir, paste0(filename, ".png")),
         plot = plot_object, width = width, height = height, dpi = 300)
  ggsave(filename = file.path(annotation_plot_dir, paste0(filename, ".pdf")),
         plot = plot_object, width = width, height = height)
}


# ------------------------------------------------------------
# 30. Figure 1: overall annotation class
# ------------------------------------------------------------

annotation_class_df <- as.data.frame(table(consensus_df$annotation_class))
colnames(annotation_class_df) <- c("class", "count")
annotation_class_df <- annotation_class_df[order(annotation_class_df$count, decreasing = TRUE), ]
annotation_class_df$class <- factor(annotation_class_df$class, levels = rev(annotation_class_df$class))

p1 <- ggplot(annotation_class_df, aes(x = class, y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(title = "Annotation class of consensus uTARs", x = "Class", y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p1)
save_annotation_plot(p1, "01_uTAR_annotation_class")


# ------------------------------------------------------------
# 31. Figure 2: evidence per database
# ------------------------------------------------------------

database_count_df <- data.frame(
  database = c("GENCODE_lncRNA", external_db_names),
  count = c(
    sum(consensus_df$GENCODE_lncRNA),
    vapply(external_db_names, function(db) { sum(consensus_df[[db]], na.rm = TRUE) }, numeric(1))
  )
)

database_count_df <- database_count_df[order(database_count_df$count, decreasing = TRUE), ]
database_count_df$database <- factor(database_count_df$database, levels = rev(database_count_df$database))

p2 <- ggplot(database_count_df, aes(x = database, y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(title = "lncRNA evidence per database", x = "Database", y = "uTARs with evidence") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p2)
save_annotation_plot(p2, "02_lncRNA_evidence_by_database")


# ------------------------------------------------------------
# 32. Figure 3: number of evidence sources
# ------------------------------------------------------------

source_count_df <- as.data.frame(table(consensus_df$n_lncRNA_sources))
colnames(source_count_df) <- c("n_sources", "count")
source_count_df$n_sources <- as.numeric(as.character(source_count_df$n_sources))

p3 <- ggplot(source_count_df, aes(x = factor(n_sources), y = count)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Number of sources with lncRNA evidence per uTAR", x = "Number of sources",
       y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)
save_annotation_plot(p3, "03_number_lncRNA_evidence_sources")


# ------------------------------------------------------------
# 33. Figure 4: candidates vs uTARs with prior evidence
# ------------------------------------------------------------

novel_df <- data.frame(
  class = c("With prior evidence", "Candidate_novel"),
  count = c(sum(!consensus_df$candidate_novel), sum(consensus_df$candidate_novel))
)

p4 <- ggplot(novel_df, aes(x = class, y = count)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "uTARs with prior evidence and unannotated candidates", x = NULL, y = "Number of uTARs") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 15, hjust = 1))

print(p4)
save_annotation_plot(p4, "04_known_vs_candidate_novel")


# ------------------------------------------------------------
# 34. Figure 5: length by annotation class
# ------------------------------------------------------------

p5 <- ggplot(consensus_df, aes(x = annotation_class, y = length)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_y_log10() +
  labs(title = "uTAR length by annotation class", x = "Annotation class", y = "Length (bp, log10)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p5)
save_annotation_plot(p5, "05_uTAR_length_by_annotation_class")


# ------------------------------------------------------------
# 35. Figure 6: recurrence of unannotated candidates in IPF
# ------------------------------------------------------------

if ("n_IPF" %in% colnames(consensus_df)) {

  novel_ipf_df <- consensus_df[consensus_df$candidate_novel, , drop = FALSE]

  novel_ipf_recurrence_df <- as.data.frame(table(novel_ipf_df$n_IPF))
  colnames(novel_ipf_recurrence_df) <- c("n_IPF", "count")
  novel_ipf_recurrence_df$n_IPF <- as.numeric(as.character(novel_ipf_recurrence_df$n_IPF))

  p6 <- ggplot(novel_ipf_recurrence_df, aes(x = factor(n_IPF), y = count)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
    labs(title = "Recurrence of unannotated candidates in IPF samples",
         x = "Number of IPF samples", y = "Number of candidates") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))

  print(p6)
  save_annotation_plot(p6, "06_candidate_novel_IPF_recurrence")
}


# ------------------------------------------------------------
# 36. Figure 7: annotation class vs IPF/control presence
# ------------------------------------------------------------

if ("presence_class" %in% colnames(consensus_df)) {

  annotation_presence_df <- as.data.frame(table(consensus_df$annotation_class, consensus_df$presence_class))
  colnames(annotation_presence_df) <- c("annotation_class", "presence_class", "count")

  p7 <- ggplot(annotation_presence_df, aes(x = annotation_class, y = count, fill = presence_class)) +
    geom_col(position = "stack") +
    labs(title = "Annotation and presence pattern in the pilot", x = "Annotation class",
         y = "Number of uTARs", fill = "Presence") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 20, hjust = 1))

  print(p7)
  save_annotation_plot(p7, "07_annotation_class_vs_presence")
}


# ------------------------------------------------------------
# 37. Figure 8: unannotated candidates >= 200 bp
# ------------------------------------------------------------

candidate_length_df <- data.frame(
  class = c("Candidate_novel <200 bp", "Candidate_novel >=200 bp"),
  count = c(sum(consensus_df$candidate_novel & consensus_df$length < 200),
            sum(consensus_df$candidate_novel_ge200))
)

p8 <- ggplot(candidate_length_df, aes(x = class, y = count)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
  labs(title = "Length of unannotated candidates", x = NULL, y = "Number of candidates") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 15, hjust = 1))

print(p8)
save_annotation_plot(p8, "08_candidate_novel_200bp")


# ------------------------------------------------------------
# 38. Figure 9: candidates recurrent in IPF
# ------------------------------------------------------------

if ("n_IPF" %in% colnames(consensus_df)) {

  recurrent_df <- data.frame(
    class = c("Other candidates", "Candidate_novel >=200 bp\nrecurrent in >=2 IPF"),
    count = c(sum(consensus_df$candidate_novel) - sum(consensus_df$candidate_IPF_recurrent),
              sum(consensus_df$candidate_IPF_recurrent))
  )

  p9 <- ggplot(recurrent_df, aes(x = class, y = count)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = format(count, big.mark = ",")), vjust = -0.4, size = 4) +
    labs(title = "Preliminary prioritisation of IPF-recurrent candidates", x = NULL,
         y = "Number of candidates") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 10, hjust = 1))

  print(p9)
  save_annotation_plot(p9, "09_candidate_novel_IPF_recurrent")
}


# ------------------------------------------------------------
# 39. Figure 10: GENCODE lncRNA biotypes
# ------------------------------------------------------------

if (nrow(gencode_lnc_hits_df) > 0) {

  gene_type_df <- as.data.frame(table(gencode_lnc_hits_df$gene_type))
  colnames(gene_type_df) <- c("gene_type", "count")
  gene_type_df <- gene_type_df[order(gene_type_df$count, decreasing = TRUE), ]
  gene_type_df$gene_type <- factor(gene_type_df$gene_type, levels = rev(gene_type_df$gene_type))

  p10 <- ggplot(gene_type_df, aes(x = gene_type, y = count)) +
    geom_col(width = 0.7) +
    coord_flip() +
    labs(title = "GENCODE biotypes overlapping uTARs", x = "Biotype", y = "Number of overlaps") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))

  print(p10)
  save_annotation_plot(p10, "10_GENCODE_lncRNA_biotypes", width = 9, height = 7)
}


# ------------------------------------------------------------
# 40. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("uTAR ANNOTATION FINISHED\n")
cat("====================================================\n")

cat("\nTotal consensus uTARs:", nrow(consensus_df), "\n")

cat("\nClassification:\n")
print(table(consensus_df$annotation_class))

cat("\nGENCODE lncRNA:", sum(consensus_df$GENCODE_lncRNA), "\n")

for (db in external_db_names) {
  cat(db, ":", sum(consensus_df[[db]], na.rm = TRUE), "\n")
}

cat("\nCandidate_novel:", sum(consensus_df$candidate_novel), "\n")
cat("Candidate_novel >=200 bp:", sum(consensus_df$candidate_novel_ge200), "\n")

if ("n_IPF" %in% colnames(consensus_df)) {
  cat("Candidate_novel >=200 bp recurrent in >=2 IPF:",
      sum(consensus_df$candidate_IPF_recurrent), "\n")
}

cat("\nMain annotated table:\n", annotated_file, "\n")

cat("\nFiles written:\n")
print(list.files(annotation_output_dir))

cat("\nFigures:\n")
print(list.files(annotation_plot_dir))

cat("\n====================================================\n")
cat("INTERPRETATION\n")
cat("====================================================\n")
cat(
  "\nGENCODE_lncRNA:\n",
  "  uTAR that now overlaps a GENCODE v47 lncRNA.\n\n",
  "External_lncRNA:\n",
  "  not a GENCODE lncRNA, but supported by at least\n",
  "  one loaded external database.\n\n",
  "GENCODE_other:\n",
  "  overlaps another GENCODE v47 gene type; review it\n",
  "  before treating it as a novel candidate.\n\n",
  "Candidate_novel:\n",
  "  no evidence in the loaded annotations.\n",
  "  It is only an unannotated candidate and needs\n",
  "  further characterisation.\n"
)

cat("\n====================================================\n")
cat("NEXT STEP\n")
cat("====================================================\n")
cat(
  "\nCharacterise the\n",
  "Candidate_novel >=200 bp loci with:\n",
  "  1. FASTA sequence extraction\n",
  "  2. CPAT\n",
  "  3. RNAsamba\n",
  "  4. TSS/PAS/splicing evidence\n",
  "  5. prioritisation by recurrence and cell type\n"
)
