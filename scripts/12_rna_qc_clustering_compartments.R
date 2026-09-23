# ============================================================
# 12_rna_qc_clustering_compartments.R
#
# lncRNA-IPF-sc | Step 12: RNA QC, clustering and manual compartment annotation
# Author: Jose A. Ovando-Ricardez
#
# Processes the conventional 10x RNA of the pilot (VUILD53, VUILD63,
# VUILD64 = IPF; VUHD67 = Control): QC, merge, SCTransform, PCA,
# clustering (resolution 0.10, 11 clusters), UMAP and FindAllMarkers.
# Clusters are then manually assigned to broad compartments so that TAR
# expression can be summarised per compartment downstream (barcode_10x is
# the key linking RNA cells to TAR barcodes).
#
# Manual cluster -> compartment mapping (pilot annotation):
#   0 Myeloid              4 Airway_Epithelial    8  Myeloid
#   1 Lymphoid             5 Lymphoid             9  Mesenchymal
#   2 Alveolar_Epithelial  6 Endothelial          10 Lymphoid
#   3 Airway_Epithelial    7 Airway_Epithelial
#
# Note: this annotation is PILOT-level. config.R is sourced after
# rm(list = ls()) so that PROJECT_DIR is not removed.
#
# Inputs:  data/<sample>/cellranger/{matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz}
#          data/<sample>/TAR/TAR_feature_bc_matrix/barcodes.tsv.gz
# Outputs: data/single_cell_RNA_Habermann/ (QC summaries, markers, compartment
#          counts, RNA_cell_metadata_for_TAR_mapping.csv,
#          IPF_pilot_RNA_manual_compartments_res010.rds),
#          data/plots/05_single_cell_RNA_Habermann/
# ============================================================

rm(list = ls())
gc()

source("config/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# ------------------------------------------------------------
# 1. future settings (avoids the 500 MiB globals limit in SCTransform)
# ------------------------------------------------------------

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
  options(future.globals.maxSize = 8 * 1024^3)
  cat("\nfuture set to sequential mode.\n")
  cat("future.globals.maxSize = 8 GiB\n\n")
}


# ------------------------------------------------------------
# 2. Directories
# ------------------------------------------------------------

BASE_DIR <- PROJECT_DIR

DATA_DIR   <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(DATA_DIR, "single_cell_RNA_Habermann")
PLOT_DIR   <- file.path(DATA_DIR, "plots", "05_single_cell_RNA_Habermann")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Samples
# ------------------------------------------------------------

samples <- c("VUILD53", "VUILD63", "VUILD64", "VUHD67")

conditions <- c(
  VUILD53 = "IPF",
  VUILD63 = "IPF",
  VUILD64 = "IPF",
  VUHD67  = "Control"
)

matrix_dirs <- setNames(file.path(DATA_DIR, samples, "cellranger"), samples)


# ------------------------------------------------------------
# 4. Check Cell Ranger files
# ------------------------------------------------------------

required_files <- c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")

file_check_list <- list()

for (sample_name in samples) {
  sample_dir <- matrix_dirs[sample_name]
  paths <- file.path(sample_dir, required_files)
  file_check_list[[sample_name]] <- data.frame(
    sample = sample_name, file = required_files, path = paths,
    exists = file.exists(paths), stringsAsFactors = FALSE
  )
}

file_check <- dplyr::bind_rows(file_check_list)

cat("\n")
cat("====================================================\n")
cat("CELL RANGER FILES\n")
cat("====================================================\n\n")

print(file_check)

if (any(!file_check$exists)) {
  stop("Missing Cell Ranger files.")
}


# ------------------------------------------------------------
# 5. Read RNA matrices
# ------------------------------------------------------------

seurat_list <- list()
preQC_list <- list()

for (sample_name in samples) {

  cat("\nProcessing:", sample_name, "\n")

  counts <- Read10X(data.dir = matrix_dirs[sample_name])

  # Keep Gene Expression if Read10X returns several modalities
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  original_barcodes <- colnames(counts)

  # Core 10x barcode: AAACCTGAGGCTAGCA-1 -> AAACCTGAGGCTAGCA
  barcode_10x <- sub("-[0-9]+$", "", original_barcodes)

  obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 0, min.features = 0)

  obj$sample <- sample_name
  obj$condition <- unname(conditions[sample_name])

  # Store barcodes BEFORE RenameCells
  obj$barcode_cellranger <- unname(original_barcodes)
  obj$barcode_10x <- unname(barcode_10x)

  # Prefix cell names with the sample to avoid duplicates
  obj <- RenameCells(object = obj, add.cell.id = sample_name)

  obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^MT-")

  preQC_list[[sample_name]] <- data.frame(
    sample = sample_name,
    condition = unname(conditions[sample_name]),
    cells_before = ncol(obj),
    median_nFeature = median(obj$nFeature_RNA, na.rm = TRUE),
    median_nCount = median(obj$nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(obj$percent.mt, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  seurat_list[[sample_name]] <- obj
}


# ------------------------------------------------------------
# 6. Pre-QC summary
# ------------------------------------------------------------

preQC_summary <- dplyr::bind_rows(preQC_list)

cat("\n")
cat("====================================================\n")
cat("PRE-QC\n")
cat("====================================================\n\n")

print(preQC_summary)

write.csv(preQC_summary, file.path(OUTPUT_DIR, "RNA_preQC_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 7. QC plots
# ------------------------------------------------------------

for (sample_name in samples) {
  obj <- seurat_list[[sample_name]]
  p_qc <- VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0)
  ggsave(filename = file.path(PLOT_DIR, paste0("01_QC_pre_", sample_name, ".png")),
         plot = p_qc, width = 12, height = 4, dpi = 300)
}


# ------------------------------------------------------------
# 8. QC filter (pilot criteria: nFeature_RNA >= 1000, percent.mt <= 25)
# ------------------------------------------------------------

qc_list <- list()
qc_summary_list <- list()

for (sample_name in samples) {

  obj <- seurat_list[[sample_name]]

  before <- ncol(obj)
  obj <- subset(x = obj, subset = nFeature_RNA >= 1000 & percent.mt <= 25)
  after <- ncol(obj)

  qc_summary_list[[sample_name]] <- data.frame(
    sample = sample_name,
    condition = unname(conditions[sample_name]),
    cells_before = before,
    cells_after = after,
    removed = before - after,
    retained_pct = 100 * after / before,
    median_nFeature = median(obj$nFeature_RNA, na.rm = TRUE),
    median_nCount = median(obj$nCount_RNA, na.rm = TRUE),
    median_percent_mt = median(obj$percent.mt, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  qc_list[[sample_name]] <- obj
}

qc_summary <- dplyr::bind_rows(qc_summary_list)

cat("\n")
cat("====================================================\n")
cat("POST-QC\n")
cat("====================================================\n\n")

print(qc_summary)

write.csv(qc_summary, file.path(OUTPUT_DIR, "RNA_postQC_summary.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 9. Merge the 4 samples (no Harmony/CCA integration)
# ------------------------------------------------------------

combined <- merge(x = qc_list[[1]], y = qc_list[2:length(qc_list)])

combined <- JoinLayers(object = combined)

DefaultAssay(combined) <- "RNA"

cat("\nCells after QC:", ncol(combined), "\n")
cat("Genes:", nrow(combined), "\n")


# ------------------------------------------------------------
# 10. SCTransform
# ------------------------------------------------------------

gc()

combined <- SCTransform(object = combined, assay = "RNA", new.assay.name = "SCT",
                        variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = TRUE)

DefaultAssay(combined) <- "SCT"


# ------------------------------------------------------------
# 11. PCA and elbow plot
# ------------------------------------------------------------

combined <- RunPCA(object = combined, assay = "SCT", features = VariableFeatures(combined),
                   npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(object = combined, ndims = 50)

ggsave(filename = file.path(PLOT_DIR, "02_ElbowPlot.png"), plot = p_elbow, width = 7, height = 5, dpi = 300)


# ------------------------------------------------------------
# 12. Neighbors
# ------------------------------------------------------------

N_PCS <- 30

combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:N_PCS, verbose = FALSE)


# ------------------------------------------------------------
# 13. Clustering (final pilot resolution 0.10 -> 11 clusters)
# ------------------------------------------------------------

set.seed(1234)

combined <- FindClusters(object = combined, resolution = 0.10, random.seed = 1234, verbose = FALSE)

combined$pilot_clusters <- unname(as.character(combined$seurat_clusters))

Idents(combined) <- "pilot_clusters"

cat("\n")
cat("====================================================\n")
cat("CLUSTERS RESOLUTION 0.10\n")
cat("====================================================\n\n")

print(table(combined$pilot_clusters))

cat("\nNumber of clusters:", length(unique(combined$pilot_clusters)), "\n")


# ------------------------------------------------------------
# 14. UMAP
# ------------------------------------------------------------

set.seed(1234)

combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:N_PCS, seed.use = 1234, verbose = FALSE)

p_clusters <- DimPlot(object = combined, reduction = "umap", group.by = "pilot_clusters",
                      label = TRUE, repel = TRUE) +
  ggtitle("Pilot clusters - resolution 0.10")

ggsave(filename = file.path(PLOT_DIR, "03_UMAP_clusters_res010.png"),
       plot = p_clusters, width = 10, height = 8, dpi = 300)


# ------------------------------------------------------------
# 15. Cluster markers (to reproduce the annotation decision)
# ------------------------------------------------------------

DefaultAssay(combined) <- "SCT"

Idents(combined) <- "pilot_clusters"

combined <- PrepSCTFindMarkers(combined, verbose = FALSE)

markers_res010 <- FindAllMarkers(object = combined, assay = "SCT", only.pos = TRUE,
                                 min.pct = 0.20, logfc.threshold = 0.25, verbose = FALSE)


# ------------------------------------------------------------
# 16. Top 10 markers per cluster
# ------------------------------------------------------------

fc_column <- if ("avg_log2FC" %in% colnames(markers_res010)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

top10_res010 <- markers_res010 %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = .data[[fc_column]], n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

write.csv(top10_res010, file.path(OUTPUT_DIR, "RNA_res010_top10_markers.csv"), row.names = FALSE)

top10_compact <- top10_res010 %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop")

cat("\n")
cat("====================================================\n")
cat("TOP GENES PER CLUSTER\n")
cat("====================================================\n\n")

print(tibble::as_tibble(top10_compact), n = Inf, width = Inf)

write.csv(top10_compact, file.path(OUTPUT_DIR, "RNA_res010_top10_markers_compact.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 17. Manual annotation of the 11 clusters
# ------------------------------------------------------------
# Based on the markers obtained:
#   0  C1QA/C1QB/C1QC            -> Myeloid
#   1  T/NK markers              -> Lymphoid
#   2  AGER/CLDN18               -> Alveolar epithelial
#   3  ciliated                  -> Airway epithelial
#   4  KRT17/SCGB1A1/BPIFB1      -> Airway epithelial
#   5  GNLY/TRDC/IFNG            -> Lymphoid
#   6  ACKR1/VWF/CLDN5           -> Endothelial
#   7  epithelial / airway-like  -> Airway epithelial
#   8  CHIT1/SPP1/IL1B           -> Myeloid
#   9  LUM/DCN/COL1A2/COL3A1     -> Mesenchymal
#   10 immunoglobulins/JCHAIN    -> Lymphoid

cluster_annotation <- c(
  "0"  = "Myeloid",
  "1"  = "Lymphoid",
  "2"  = "Alveolar_Epithelial",
  "3"  = "Airway_Epithelial",
  "4"  = "Airway_Epithelial",
  "5"  = "Lymphoid",
  "6"  = "Endothelial",
  "7"  = "Airway_Epithelial",
  "8"  = "Myeloid",
  "9"  = "Mesenchymal",
  "10" = "Lymphoid"
)


# ------------------------------------------------------------
# 18. Transfer annotation
# ------------------------------------------------------------

manual_annotation <- unname(cluster_annotation[as.character(combined$pilot_clusters)])

if (any(is.na(manual_annotation))) {
  stop("Some clusters have no manual annotation.")
}

combined$cell_compartment <- manual_annotation

combined$cell_compartment <- factor(
  combined$cell_compartment,
  levels = c("Airway_Epithelial", "Alveolar_Epithelial", "Mesenchymal",
             "Endothelial", "Myeloid", "Lymphoid")
)


# ------------------------------------------------------------
# 19. Check annotation
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CLUSTER -> COMPARTMENT\n")
cat("====================================================\n\n")

print(table(combined$pilot_clusters, combined$cell_compartment))

cat("\n")
cat("====================================================\n")
cat("CELLS PER COMPARTMENT\n")
cat("====================================================\n\n")

print(table(combined$cell_compartment))


# ------------------------------------------------------------
# 20. Compartment UMAPs
# ------------------------------------------------------------

p_compartments <- DimPlot(object = combined, reduction = "umap", group.by = "cell_compartment",
                          label = TRUE, repel = TRUE) +
  ggtitle("Broad cell compartments")

ggsave(filename = file.path(PLOT_DIR, "04_UMAP_manual_compartments.png"),
       plot = p_compartments, width = 9, height = 7, dpi = 300)

p_split <- DimPlot(object = combined, reduction = "umap", group.by = "cell_compartment",
                   split.by = "sample", ncol = 2)

ggsave(filename = file.path(PLOT_DIR, "05_UMAP_manual_compartments_by_sample.png"),
       plot = p_split, width = 14, height = 10, dpi = 300)


# ------------------------------------------------------------
# 21. Validation DotPlot
# ------------------------------------------------------------

validation_markers <- c(
  # Airway
  "KRT5", "KRT17", "TP63", "FOXJ1", "SCGB1A1", "SCGB3A1", "BPIFB1",
  # Alveolar
  "AGER", "CLDN18", "SFTPC", "SFTPA1", "SFTPB", "ABCA3",
  # Mesenchymal
  "LUM", "DCN", "COL1A1", "COL1A2", "COL3A1",
  # Endothelial
  "ACKR1", "PECAM1", "VWF", "CLDN5",
  # Myeloid
  "C1QA", "C1QB", "C1QC", "LST1", "TYROBP", "SPP1",
  # Lymphoid
  "CD3D", "CD3E", "TRBC1", "NKG7", "GNLY", "MS4A1", "JCHAIN"
)

validation_markers <- intersect(validation_markers, rownames(combined))

p_dot <- DotPlot(object = combined, features = validation_markers, group.by = "cell_compartment") +
  RotatedAxis() +
  ggtitle("Manual compartment validation")

ggsave(filename = file.path(PLOT_DIR, "06_DotPlot_manual_compartments.png"),
       plot = p_dot, width = 18, height = 7, dpi = 300)


# ------------------------------------------------------------
# 22. Counts per sample
# ------------------------------------------------------------

compartment_counts <- combined@meta.data %>%
  dplyr::count(sample, condition, cell_compartment, name = "n_cells")

compartment_proportions <- compartment_counts %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(total_cells = sum(n_cells), percentage = 100 * n_cells / total_cells) %>%
  dplyr::ungroup()

write.csv(compartment_counts, file.path(OUTPUT_DIR, "RNA_manual_compartment_counts.csv"), row.names = FALSE)
write.csv(compartment_proportions, file.path(OUTPUT_DIR, "RNA_manual_compartment_proportions.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 23. Check barcode_10x (the key that links RNA to TAR)
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("BARCODE 10x EXAMPLES\n")
cat("====================================================\n\n")

print(
  combined@meta.data %>%
    dplyr::select(sample, barcode_cellranger, barcode_10x) %>%
    head(10)
)


# ------------------------------------------------------------
# 24. Check RNA -> TAR barcode match
# ------------------------------------------------------------

barcode_check_list <- list()

for (sample_name in samples) {

  tar_file <- file.path(DATA_DIR, sample_name, "TAR", "TAR_feature_bc_matrix", "barcodes.tsv.gz")

  if (!file.exists(tar_file)) {
    warning(paste("File not found:", tar_file))
    next
  }

  tar_barcodes <- readLines(gzfile(tar_file))

  # Strip the -1 suffix from TAR barcodes
  tar_barcode_10x <- sub("-[0-9]+$", "", tar_barcodes)

  rna_barcode_10x <- combined@meta.data %>%
    dplyr::filter(sample == sample_name) %>%
    dplyr::pull(barcode_10x)

  matches <- sum(rna_barcode_10x %in% tar_barcode_10x)

  barcode_check_list[[sample_name]] <- data.frame(
    sample = sample_name,
    RNA_cells = length(rna_barcode_10x),
    TAR_barcodes = length(tar_barcode_10x),
    RNA_TAR_matches = matches,
    match_pct = 100 * matches / length(rna_barcode_10x),
    stringsAsFactors = FALSE
  )
}

barcode_check <- dplyr::bind_rows(barcode_check_list)

cat("\n")
cat("====================================================\n")
cat("RNA -> TAR MATCH\n")
cat("====================================================\n\n")

print(barcode_check)

write.csv(barcode_check, file.path(OUTPUT_DIR, "RNA_TAR_barcode_match_manual_annotation.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 25. Export cell metadata for the TAR scripts
# ------------------------------------------------------------

cell_metadata <- combined@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(cell, barcode_cellranger, barcode_10x, sample, condition,
                nCount_RNA, nFeature_RNA, percent.mt, pilot_clusters, cell_compartment)

write.csv(cell_metadata, file.path(OUTPUT_DIR, "RNA_cell_metadata_for_TAR_mapping.csv"), row.names = FALSE)


# ------------------------------------------------------------
# 26. Save cluster -> compartment map
# ------------------------------------------------------------

cluster_annotation_table <- data.frame(
  cluster = names(cluster_annotation),
  cell_compartment = unname(cluster_annotation),
  stringsAsFactors = FALSE
)

write.csv(cluster_annotation_table, file.path(OUTPUT_DIR, "RNA_cluster_manual_annotation_res010.csv"),
          row.names = FALSE)


# ------------------------------------------------------------
# 27. Save final object
# ------------------------------------------------------------

saveRDS(combined, file.path(OUTPUT_DIR, "IPF_pilot_RNA_manual_compartments_res010.rds"))


# ------------------------------------------------------------
# 28. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINAL SUMMARY STEP 12\n")
cat("====================================================\n\n")

cat("Cells:", ncol(combined), "\n")
cat("Genes:", nrow(combined), "\n")
cat("Clusters:", length(unique(combined$pilot_clusters)), "\n\n")

cat("CLUSTERS:\n\n")
print(table(combined$pilot_clusters))

cat("\nCOMPARTMENTS:\n\n")
print(table(combined$cell_compartment))

cat("\nBY SAMPLE:\n\n")
print(table(combined$sample, combined$cell_compartment))

cat("\nRNA -> TAR:\n\n")
print(barcode_check)

cat("\n")
cat("====================================================\n")
cat("STEP 12 FINISHED\n")
cat("====================================================\n\n")

cat("Final object:\n", file.path(OUTPUT_DIR, "IPF_pilot_RNA_manual_compartments_res010.rds"), "\n\n")

cat("Metadata for TAR:\n", file.path(OUTPUT_DIR, "RNA_cell_metadata_for_TAR_mapping.csv"), "\n")
