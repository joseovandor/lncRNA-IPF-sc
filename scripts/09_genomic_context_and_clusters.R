# ============================================================
# 09_genomic_context_and_clusters.R
#
# lncRNA-IPF-sc | Step 9: genomic context and same-strand TAR clusters
# Author: Jose A. Ovando-Ricardez
#
# Characterises the genomic context of the unannotated candidates present
# in 3/3 IPF samples relative to GENCODE v47: gene/exon overlap, strand
# relationship (same_strand_exonic, same_strand_intronic, antisense_exonic,
# antisense_gene_body, intergenic), nearest gene, genes within 10/50/100 kb,
# and same-strand TAR clusters separated by <= 5 kb.
#
# Note: the classification describes the context of the TAR interval only;
# it does not demonstrate a mature transcript structure, and clusters are
# not necessarily a single transcript.
#
# Inputs:  data/TAR_catalog/coding_potential/Unannotated_3of3_IPF_candidates_for_coding_potential.csv
#          data/annotations/**/*gencode*v47*.gtf(.gz) (auto-detected; first match is used)
# Outputs: data/TAR_catalog/genomic_context/
# ============================================================

source("config/config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
})


# ------------------------------------------------------------
# 1. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

INPUT_FILE <- file.path(
  BASE_DIR,
  "data/TAR_catalog/coding_potential/Unannotated_3of3_IPF_candidates_for_coding_potential.csv"
)

OUTPUT_DIR <- file.path(BASE_DIR, "data/TAR_catalog/genomic_context")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Locate GENCODE v47
# ------------------------------------------------------------

ANNOTATION_DIR <- file.path(BASE_DIR, "data/annotations")

gtf_candidates <- list.files(ANNOTATION_DIR, pattern = "\\.gtf(\\.gz)?$", full.names = TRUE,
                             recursive = TRUE, ignore.case = TRUE)

gencode_candidates <- gtf_candidates[
  grepl("gencode", basename(gtf_candidates), ignore.case = TRUE) &
    grepl("v47|47", basename(gtf_candidates), ignore.case = TRUE)
]

cat("\n")
cat("====================================================\n")
cat("GENCODE GTF\n")
cat("====================================================\n\n")

if (length(gencode_candidates) == 0) {
  cat("Could not automatically find a GENCODE v47 GTF.\n\n")
  cat("GTFs available under data/annotations:\n")
  print(gtf_candidates)
  stop("\nSet the GENCODE v47 path manually in GENCODE_GTF.")
}

if (length(gencode_candidates) > 1) {
  cat("Found several GENCODE v47 candidates:\n\n")
  print(gencode_candidates)
  cat("\nUsing the first one. Check that it is correct.\n\n")
}

GENCODE_GTF <- gencode_candidates[1]

cat("GENCODE used:\n", GENCODE_GTF, "\n\n")


# ------------------------------------------------------------
# 3. Check input
# ------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(paste("File not found:", INPUT_FILE))
}


# ------------------------------------------------------------
# 4. Load candidates
# ------------------------------------------------------------

candidates <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)

cat("====================================================\n")
cat("CANDIDATES\n")
cat("====================================================\n\n")

cat("Total candidates:", nrow(candidates), "\n")
cat("Priority A:", sum(candidates$priority_group == "Priority_A", na.rm = TRUE), "\n")
cat("High confidence:", sum(candidates$high_confidence_pilot, na.rm = TRUE), "\n\n")


# ------------------------------------------------------------
# 5. Candidate GRanges
# ------------------------------------------------------------

tar_gr <- GRanges(seqnames = candidates$chr,
                  ranges = IRanges(start = candidates$start, end = candidates$end),
                  strand = candidates$strand)

names(tar_gr) <- candidates$TAR_ID

mcols(tar_gr)$TAR_ID <- candidates$TAR_ID


# ------------------------------------------------------------
# 6. Import GENCODE; keep genes and exons
# ------------------------------------------------------------

cat("Importing GENCODE...\n")

gencode <- rtracklayer::import(GENCODE_GTF)

cat("GENCODE features:", length(gencode), "\n\n")

genes <- gencode[gencode$type == "gene"]

exons <- gencode[gencode$type == "exon"]

cat("Genes:", length(genes), "\n")
cat("Exons:", length(exons), "\n\n")


# ------------------------------------------------------------
# 7. Annotation columns
# ------------------------------------------------------------

gene_name_column <- if ("gene_name" %in% colnames(mcols(genes))) {
  "gene_name"
} else {
  NA_character_
}

gene_id_column <- if ("gene_id" %in% colnames(mcols(genes))) {
  "gene_id"
} else {
  NA_character_
}

gene_type_column <- if ("gene_type" %in% colnames(mcols(genes))) {
  "gene_type"
} else if ("gene_biotype" %in% colnames(mcols(genes))) {
  "gene_biotype"
} else {
  NA_character_
}

cat("gene_name column:", gene_name_column, "\n")
cat("gene_id column:", gene_id_column, "\n")
cat("gene_type column:", gene_type_column, "\n\n")


# ------------------------------------------------------------
# 8. Chromosome compatibility (TAR vs GENCODE)
# ------------------------------------------------------------

tar_chr <- unique(as.character(seqnames(tar_gr)))

gene_chr <- unique(as.character(seqnames(genes)))

common_chr <- intersect(tar_chr, gene_chr)

cat("TAR chromosomes:", length(tar_chr), "\n")
cat("GENCODE chromosomes:", length(gene_chr), "\n")
cat("Shared:", length(common_chr), "\n\n")

if (length(common_chr) == 0) {
  stop("No compatible chromosomes between TAR and GENCODE.")
}


# ------------------------------------------------------------
# 9. Keep shared seqlevels
# ------------------------------------------------------------

tar_gr <- GenomeInfoDb::keepSeqlevels(tar_gr, common_chr, pruning.mode = "coarse")

genes <- GenomeInfoDb::keepSeqlevels(genes, common_chr, pruning.mode = "coarse")

exons <- GenomeInfoDb::keepSeqlevels(exons, common_chr, pruning.mode = "coarse")


# ------------------------------------------------------------
# 10. Gene metadata
# ------------------------------------------------------------

gene_metadata <- data.frame(
  gene_index = seq_along(genes),
  gene_id = if (!is.na(gene_id_column)) {
    as.character(mcols(genes)[[gene_id_column]])
  } else {
    NA_character_
  },
  gene_name = if (!is.na(gene_name_column)) {
    as.character(mcols(genes)[[gene_name_column]])
  } else {
    NA_character_
  },
  gene_type = if (!is.na(gene_type_column)) {
    as.character(mcols(genes)[[gene_type_column]])
  } else {
    NA_character_
  },
  gene_chr = as.character(seqnames(genes)),
  gene_start = start(genes),
  gene_end = end(genes),
  gene_strand = as.character(strand(genes)),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 11. Overlaps with genes and exons (same strand / any strand)
# ------------------------------------------------------------

hits_gene_same <- findOverlaps(tar_gr, genes, ignore.strand = FALSE)

hits_gene_any <- findOverlaps(tar_gr, genes, ignore.strand = TRUE)

hits_exon_same <- findOverlaps(tar_gr, exons, ignore.strand = FALSE)

hits_exon_any <- findOverlaps(tar_gr, exons, ignore.strand = TRUE)


# ------------------------------------------------------------
# 12. Logical overlap vectors
# ------------------------------------------------------------

n_tar <- length(tar_gr)

same_gene <- rep(FALSE, n_tar)

same_gene[unique(queryHits(hits_gene_same))] <- TRUE

any_gene <- rep(FALSE, n_tar)

any_gene[unique(queryHits(hits_gene_any))] <- TRUE

same_exon <- rep(FALSE, n_tar)

same_exon[unique(queryHits(hits_exon_same))] <- TRUE

any_exon <- rep(FALSE, n_tar)

any_exon[unique(queryHits(hits_exon_any))] <- TRUE

# Antisense = overlap on the opposite strand only
antisense_gene <- any_gene & !same_gene

antisense_exon <- any_exon & !same_exon


# ------------------------------------------------------------
# 13. Genomic context classification
# ------------------------------------------------------------
# Precedence: same-strand exon > same-strand gene body (no exon) >
# antisense exon > antisense gene body > intergenic

genomic_context <- case_when(
  same_exon ~ "same_strand_exonic",
  same_gene & !same_exon ~ "same_strand_intronic",
  antisense_exon ~ "antisense_exonic",
  antisense_gene ~ "antisense_gene_body",
  TRUE ~ "intergenic"
)


# ------------------------------------------------------------
# 14. Overlapping genes (counts, names, IDs)
# ------------------------------------------------------------

n_gene_overlaps <- tabulate(queryHits(hits_gene_any), nbins = n_tar)

n_gene_overlaps_same_strand <- tabulate(queryHits(hits_gene_same), nbins = n_tar)

overlap_gene_names <- rep(NA_character_, n_tar)

overlap_gene_ids <- rep(NA_character_, n_tar)

if (length(hits_gene_any) > 0) {
  hit_table <- data.frame(TAR_index = queryHits(hits_gene_any),
                          gene_index = subjectHits(hits_gene_any)) %>%
    left_join(gene_metadata, by = "gene_index")

  hit_names <- hit_table %>%
    group_by(TAR_index) %>%
    summarise(gene_names = paste(unique(na.omit(gene_name)), collapse = ";"),
              gene_ids = paste(unique(na.omit(gene_id)), collapse = ";"),
              .groups = "drop")

  overlap_gene_names[hit_names$TAR_index] <- hit_names$gene_names

  overlap_gene_ids[hit_names$TAR_index] <- hit_names$gene_ids
}


# ------------------------------------------------------------
# 15. Nearest gene (strand-agnostic)
# ------------------------------------------------------------

nearest_hits <- distanceToNearest(tar_gr, genes, ignore.strand = TRUE)

nearest_df <- data.frame(
  TAR_index = queryHits(nearest_hits),
  gene_index = subjectHits(nearest_hits),
  nearest_gene_distance = mcols(nearest_hits)$distance,
  stringsAsFactors = FALSE
) %>%
  left_join(gene_metadata, by = "gene_index")

nearest_gene_name <- rep(NA_character_, n_tar)

nearest_gene_id <- rep(NA_character_, n_tar)

nearest_gene_type <- rep(NA_character_, n_tar)

nearest_gene_distance <- rep(NA_integer_, n_tar)

nearest_gene_chr <- rep(NA_character_, n_tar)

nearest_gene_start <- rep(NA_integer_, n_tar)

nearest_gene_end <- rep(NA_integer_, n_tar)

nearest_gene_strand <- rep(NA_character_, n_tar)

nearest_gene_name[nearest_df$TAR_index] <- nearest_df$gene_name

nearest_gene_id[nearest_df$TAR_index] <- nearest_df$gene_id

nearest_gene_type[nearest_df$TAR_index] <- nearest_df$gene_type

nearest_gene_distance[nearest_df$TAR_index] <- nearest_df$nearest_gene_distance

nearest_gene_chr[nearest_df$TAR_index] <- nearest_df$gene_chr

nearest_gene_start[nearest_df$TAR_index] <- nearest_df$gene_start

nearest_gene_end[nearest_df$TAR_index] <- nearest_df$gene_end

nearest_gene_strand[nearest_df$TAR_index] <- nearest_df$gene_strand


# ------------------------------------------------------------
# 16. Genes within 10 / 50 / 100 kb
# ------------------------------------------------------------

tar_10kb <- resize(tar_gr, width = width(tar_gr) + 20000, fix = "center")

hits_10kb <- findOverlaps(tar_10kb, genes, ignore.strand = TRUE)

genes_within_10kb <- tabulate(queryHits(hits_10kb), nbins = n_tar)

tar_50kb <- resize(tar_gr, width = width(tar_gr) + 100000, fix = "center")

hits_50kb <- findOverlaps(tar_50kb, genes, ignore.strand = TRUE)

genes_within_50kb <- tabulate(queryHits(hits_50kb), nbins = n_tar)

tar_100kb <- resize(tar_gr, width = width(tar_gr) + 200000, fix = "center")

hits_100kb <- findOverlaps(tar_100kb, genes, ignore.strand = TRUE)

genes_within_100kb <- tabulate(queryHits(hits_100kb), nbins = n_tar)


# ------------------------------------------------------------
# 17. Same-strand TAR clustering
# ------------------------------------------------------------
# A cluster groups TARs on the same chromosome and strand separated by
# <= 5 kb. This does NOT imply that they form a single transcript.

MAX_CLUSTER_GAP <- 5000

tar_clusters <- reduce(tar_gr, ignore.strand = FALSE, min.gapwidth = MAX_CLUSTER_GAP + 1)

names(tar_clusters) <- paste0("TARCLUSTER_", sprintf("%05d", seq_along(tar_clusters)))


# ------------------------------------------------------------
# 18. Map TARs to clusters
# ------------------------------------------------------------

cluster_hits <- findOverlaps(tar_gr, tar_clusters, ignore.strand = FALSE)

cluster_map <- data.frame(TAR_index = queryHits(cluster_hits),
                          cluster_index = subjectHits(cluster_hits))

# Each TAR must fall in exactly one cluster
if (any(duplicated(cluster_map$TAR_index))) {
  stop("A TAR was assigned to more than one cluster.")
}

cluster_id <- rep(NA_character_, n_tar)

cluster_id[cluster_map$TAR_index] <- names(tar_clusters)[cluster_map$cluster_index]


# ------------------------------------------------------------
# 19. Cluster size, span and members
# ------------------------------------------------------------

cluster_counts <- table(cluster_id)

cluster_n_TAR <- as.integer(cluster_counts[cluster_id])

cluster_span <- width(tar_clusters)

names(cluster_span) <- names(tar_clusters)

cluster_span_bp <- cluster_span[cluster_id]

cluster_member_table <- data.frame(TAR_ID = names(tar_gr), cluster_id = cluster_id,
                                   stringsAsFactors = FALSE)

cluster_members_summary <- cluster_member_table %>%
  group_by(cluster_id) %>%
  summarise(cluster_members = paste(TAR_ID, collapse = ";"), .groups = "drop")

cluster_members <- cluster_members_summary$cluster_members[
  match(cluster_id, cluster_members_summary$cluster_id)
]


# ------------------------------------------------------------
# 20. Number of same-strand neighbours within 5 kb
# ------------------------------------------------------------

pair_hits <- findOverlaps(tar_gr, tar_gr, maxgap = MAX_CLUSTER_GAP, minoverlap = 0,
                          ignore.strand = FALSE)

pair_df <- data.frame(q = queryHits(pair_hits), s = subjectHits(pair_hits))

pair_df <- pair_df[pair_df$q != pair_df$s, , drop = FALSE]

n_nearby_TAR_5kb <- tabulate(pair_df$q, nbins = n_tar)


# ------------------------------------------------------------
# 21. Context table and merge with candidates
# ------------------------------------------------------------

context_df <- data.frame(
  TAR_ID = names(tar_gr),
  genomic_context = genomic_context,
  n_gene_overlaps = n_gene_overlaps,
  n_gene_overlaps_same_strand = n_gene_overlaps_same_strand,
  overlapping_gene_names = overlap_gene_names,
  overlapping_gene_ids = overlap_gene_ids,
  nearest_gene = nearest_gene_name,
  nearest_gene_id = nearest_gene_id,
  nearest_gene_type = nearest_gene_type,
  nearest_gene_distance = nearest_gene_distance,
  nearest_gene_chr = nearest_gene_chr,
  nearest_gene_start = nearest_gene_start,
  nearest_gene_end = nearest_gene_end,
  nearest_gene_strand = nearest_gene_strand,
  genes_within_10kb = genes_within_10kb,
  genes_within_50kb = genes_within_50kb,
  genes_within_100kb = genes_within_100kb,
  cluster_id = cluster_id,
  cluster_n_TAR = cluster_n_TAR,
  cluster_span_bp = as.integer(cluster_span_bp),
  cluster_members = cluster_members,
  nearby_same_strand_TAR_5kb = n_nearby_TAR_5kb,
  stringsAsFactors = FALSE
)

result <- candidates %>%
  left_join(context_df, by = "TAR_ID")


# ------------------------------------------------------------
# 22. Genomic context summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GENOMIC CONTEXT\n")
cat("====================================================\n\n")

print(table(result$genomic_context))

cat("\nPercentage:\n")

print(round(prop.table(table(result$genomic_context)) * 100, 2))


# ------------------------------------------------------------
# 23. Priority A context
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("PRIORITY A - CONTEXT\n")
cat("====================================================\n\n")

priority_A_context <- result %>%
  filter(priority_group == "Priority_A")

print(table(priority_A_context$genomic_context))


# ------------------------------------------------------------
# 24. Cluster summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TAR CLUSTERS <=5 kb\n")
cat("====================================================\n\n")

cat("Total number of clusters:", length(tar_clusters), "\n")
cat("Clusters with >=2 TAR:", sum(cluster_counts >= 2), "\n")
cat("Clusters with >=3 TAR:", sum(cluster_counts >= 3), "\n")
cat("Maximum TAR per cluster:", max(cluster_counts), "\n\n")

print(table(result$cluster_n_TAR))


# ------------------------------------------------------------
# 25. Priority A in clusters
# ------------------------------------------------------------

priority_A_clusters <- priority_A_context %>%
  filter(cluster_n_TAR >= 2) %>%
  arrange(desc(cluster_n_TAR), cluster_id, start)

cat("\n")
cat("====================================================\n")
cat("PRIORITY A IN CLUSTERS\n")
cat("====================================================\n\n")

cat("Priority A in clusters with >=2 TAR:", nrow(priority_A_clusters), "/",
    nrow(priority_A_context), "\n\n")

print(
  priority_A_clusters %>%
    select(TAR_ID, chr, start, end, strand, length,
           genomic_context,
           nearest_gene, nearest_gene_type, nearest_gene_distance,
           cluster_id, cluster_n_TAR, cluster_span_bp, cluster_members,
           IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells)
)


# ------------------------------------------------------------
# 26. Priority A intergenic
# ------------------------------------------------------------

priority_A_intergenic <- priority_A_context %>%
  filter(genomic_context == "intergenic") %>%
  arrange(desc(IPF_min_counts), IPF_CV_counts)

cat("\n")
cat("====================================================\n")
cat("PRIORITY A INTERGENIC\n")
cat("====================================================\n\n")

cat("N:", nrow(priority_A_intergenic), "\n\n")

print(
  priority_A_intergenic %>%
    select(TAR_ID, chr, start, end, strand, length,
           nearest_gene, nearest_gene_type, nearest_gene_distance,
           genes_within_10kb, genes_within_50kb, genes_within_100kb,
           cluster_id, cluster_n_TAR,
           IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells)
)


# ------------------------------------------------------------
# 27. Priority A antisense
# ------------------------------------------------------------

priority_A_antisense <- priority_A_context %>%
  filter(grepl("^antisense", genomic_context)) %>%
  arrange(desc(IPF_min_counts))

cat("\n")
cat("====================================================\n")
cat("PRIORITY A ANTISENSE\n")
cat("====================================================\n\n")

cat("N:", nrow(priority_A_antisense), "\n\n")

print(
  priority_A_antisense %>%
    select(TAR_ID, chr, start, end, strand, length,
           genomic_context, overlapping_gene_names,
           nearest_gene, nearest_gene_type,
           IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells)
)


# ------------------------------------------------------------
# 28. Top Priority A with context
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP 30 PRIORITY A + CONTEXT\n")
cat("====================================================\n\n")

top_priority_context <- priority_A_context %>%
  arrange(desc(IPF_min_counts), IPF_CV_counts) %>%
  select(TAR_ID, chr, start, end, strand, length,
         genomic_context,
         overlapping_gene_names,
         nearest_gene, nearest_gene_type, nearest_gene_distance,
         cluster_id, cluster_n_TAR, cluster_members,
         counts_VUILD53, counts_VUILD63, counts_VUILD64, counts_VUHD67,
         IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells)

print(head(top_priority_context, 30))


# ------------------------------------------------------------
# 29. Save tables
# ------------------------------------------------------------

write.csv(result, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_genomic_context.csv"), row.names = FALSE)

write.csv(priority_A_context, file.path(OUTPUT_DIR, "Priority_A_genomic_context.csv"), row.names = FALSE)

write.csv(priority_A_intergenic, file.path(OUTPUT_DIR, "Priority_A_intergenic.csv"), row.names = FALSE)

write.csv(priority_A_antisense, file.path(OUTPUT_DIR, "Priority_A_antisense.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 30. Cluster table
# ------------------------------------------------------------

cluster_summary <- data.frame(
  cluster_id = names(tar_clusters),
  chr = as.character(seqnames(tar_clusters)),
  start = start(tar_clusters),
  end = end(tar_clusters),
  strand = as.character(strand(tar_clusters)),
  cluster_span_bp = width(tar_clusters),
  stringsAsFactors = FALSE
)

cluster_summary <- cluster_summary %>%
  left_join(cluster_members_summary, by = "cluster_id")

cluster_summary$cluster_n_TAR <- lengths(strsplit(cluster_summary$cluster_members, ";", fixed = TRUE))

cluster_summary <- cluster_summary %>%
  arrange(desc(cluster_n_TAR), chr, start)

write.csv(cluster_summary, file.path(OUTPUT_DIR, "TAR_clusters_5kb.csv"), row.names = FALSE)

write.csv(cluster_summary %>% filter(cluster_n_TAR >= 2),
          file.path(OUTPUT_DIR, "TAR_clusters_5kb_multiTAR.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 31. Save RDS
# ------------------------------------------------------------

saveRDS(tar_clusters, file.path(OUTPUT_DIR, "TAR_clusters_5kb_GRanges.rds"))

saveRDS(result, file.path(OUTPUT_DIR, "Unannotated_3of3_IPF_genomic_context.rds"))


# ------------------------------------------------------------
# 32. Previously flagged loci
# ------------------------------------------------------------

interesting_ids <- c(
  "IPFTAR049318",
  "IPFTAR053240", "IPFTAR053241", "IPFTAR053242",
  "IPFTAR032367", "IPFTAR032368", "IPFTAR032369"
)

cat("\n")
cat("====================================================\n")
cat("PREVIOUSLY FLAGGED LOCI\n")
cat("====================================================\n\n")

print(
  result %>%
    filter(TAR_ID %in% interesting_ids) %>%
    select(TAR_ID, chr, start, end, strand, length,
           genomic_context,
           overlapping_gene_names,
           nearest_gene, nearest_gene_type, nearest_gene_distance,
           cluster_id, cluster_n_TAR, cluster_members,
           counts_VUILD53, counts_VUILD63, counts_VUILD64, counts_VUHD67,
           IPF_min_counts, IPF_CV_counts, IPF_min_pct_cells) %>%
    arrange(chr, start)
)


# ------------------------------------------------------------
# 33. Generated files
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GENERATED FILES\n")
cat("====================================================\n\n")

print(list.files(OUTPUT_DIR, full.names = FALSE))


# ------------------------------------------------------------
# 34. Done
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("GENOMIC CONTEXT ANALYSIS FINISHED\n")
cat("====================================================\n\n")

cat("Output:\n", OUTPUT_DIR, "\n\n")

cat(
  "IMPORTANT:\n",
  "Clusters indicate same-strand genomic proximity.\n",
  "They should not automatically be interpreted as a single transcript.\n"
)
