# lncRNA-IPF-sc

**Discovery of unannotated long non-coding RNA loci in idiopathic pulmonary fibrosis (IPF) at single-cell resolution — pilot analysis.**

Laboratorio de Biología Computacional, Instituto Nacional de Enfermedades Respiratorias (INER)
Author: Jose A. Ovando-Ricardez

---

## Overview

Most single-cell RNA-seq (scRNA-seq) pipelines only quantify reads that fall on annotated genes, so transcription from unannotated regions — including many lncRNAs — is discarded. **TAR-scRNA-seq** recovers these signals by calling *transcriptionally active regions* (TARs) directly from the aligned reads (groHMM, 50-bp bins, TARs merged within 500 bp) and quantifying them per cell, alongside the conventional Cell Ranger gene matrix.

This repository contains the R pipeline for a **pilot** that applies this strategy to lung scRNA-seq data from IPF and control donors, following the workflow suggested by Prakrithi et al. (*Nature Methods*, 2026). The pilot asks whether the approach is technically feasible on public IPF data and whether it yields reproducible, IPF-recurrent, unannotated loci worth pursuing in a full cohort.

The pipeline:

1. imports Cell Ranger and TAR matrices and runs QC;
2. builds a cross-sample **consensus catalogue of unannotated TARs (uTARs)** (`IPFTAR` IDs, same-strand overlap merge);
3. annotates the catalogue against **GENCODE v47** and external lncRNA databases (**LncBook, NONCODE, LNCipedia, FANTOM-CAT**);
4. uses known IPF-related lncRNAs (MEG3, FENDRR, MALAT1, NEAT1, H19, …) as a positive control;
5. quantifies loci present in **3/3 IPF samples**, prioritises unannotated candidates (Priority A/B/C) and extracts their sequences for coding-potential analysis;
6. characterises genomic context (antisense, intronic, intergenic, same-strand clusters) and refines Priority A loci structurally from the per-sample TAR BED files;
7. annotates cells into four broad compartments (**Epithelial, Mesenchymal, Immune, Endothelial**) and measures the compartment specificity of known and unannotated TARs.

> **Terminology.** A uTAR / `Unannotated_candidate` / `Candidate_novel` is a region with no evidence in the annotations loaded. It is **not** a confirmed novel lncRNA: length, coding potential, TSS/PAS, splicing, conservation and replication still have to be evaluated. All IPF vs control comparisons in this pilot are **descriptive** (3 IPF vs 1 control); no formal differential-expression testing is performed.

## Pilot design

| Sample  | Condition | Source |
|---------|-----------|--------|
| VUILD53 | IPF       | Habermann et al., *Sci Adv* 2020 (GEO GSE135893) |
| VUILD63 | IPF       | idem |
| VUILD64 | IPF       | idem |
| VUHD67  | Control   | idem |

Sample IDs are hard-coded in the scripts.

## Repository layout

```
lncRNA-IPF-sc/
├── README.md
├── config/
│   └── config.R          # PROJECT_DIR / DATA_DIR (sourced by every script)
├── scripts/              # analysis steps, run in numeric order
└── data/                 # NOT tracked by git (see data/README.md)
```

## Pipeline

| Step | Script | Original name | Purpose |
|------|--------|---------------|---------|
| 1  | `01_import_and_qc.R` | `00_Data_import_and_QC.R` | Import Cell Ranger + TAR matrices, barcode matching, RNA/TAR QC, aTAR (`_1`) vs uTAR (`_0`) split |
| 2  | `02_build_consensus_utar_catalog.R` | `01_Catalog_generation.R` | Cross-sample consensus uTAR catalogue (`IPFTAR` IDs), sample mapping and presence table, CSV/BED |
| 3  | `03_annotate_gencode.R` | `02_Annotation.R` | Strand-specific annotation against GENCODE v47 (lncRNA / other genes) |
| 4  | `04_known_ipf_lncrna_detection.R` | `03_lncRNA_conocidos_en_IPF.R` | Detect known IPF-related lncRNAs in the RNA matrix and TAR catalogue |
| 5  | `05_known_ipf_lncrna_exploratory.R` | `04_analisis_exploratorio_de_conocidos.R` | Descriptive IPF vs control CPM comparison of known lncRNAs (positive control) |
| 6  | `06_annotate_external_lncrna_dbs.R` | `05_external_lncRNA_annotation.R` | Overlap with LncBook, NONCODE, LNCipedia, FANTOM-CAT |
| 7  | `07_quantify_recurrent_ipf_loci.R` | `06_quantify_3of3_IPF_lncRNA.R` | Quantify loci present in 3/3 IPF (GENCODE / external / unannotated) |
| 8  | `08_prioritize_candidates_and_extract_sequences.R` | `07_prepare_3of3_candidates_for_coding_potential.R` | Priority A/B/C, BED + GRCh38 FASTA for coding-potential tools |
| 9  | `09_genomic_context_and_clusters.R` | `08_genomic_context_and_TAR_clusters.R` | Genomic context vs GENCODE v47, nearest genes, same-strand clusters (≤ 5 kb) |
| 10 | `10_refine_focus_loci_tar_bed.R` | `09_refine_priorityA_loci_from_TAR_bed.R` | Locus refinement of focus loci (ADNP, DYRK1A, TCF4 antisense; SMAD4-near intergenic) from `TAR_reads.bed.gz` |
| 11 | `11_priorityA_structural_refinement.R` | `10_systematic_refinement_PriorityA.R` | Structural reproducibility of all Priority A candidates; IGV tracks |
| 12 | `12_rna_qc_clustering_compartments.R` | `11_RNA_Habermann_manual_compartments_res010.R` | RNA QC, SCTransform, clustering (res 0.10), manual compartment annotation |
| 13 | `13_map_tar_to_single_cells.R` | `12_single_cell_TAR_compartment_mapping.R` | **Missing** — see *Known issues* |
| 14 | `14_compartment_specificity.R` | `13_broad_compartment_TAR_specificity.R` | TAR detection per sample × broad compartment, dominant compartment |
| 15 | `15_known_vs_unannotated_specificity.R` | `14_posthoc_TAR_specificity_known_vs_unannotated.R` | Post-hoc specificity: known vs unannotated TARs, descriptive shortlist |

Each script begins with a header listing its exact inputs and outputs.

## Input data

All paths are relative to `data/` (see `data/README.md` for the full layout):

- `data/<sample>/cellranger/` — Cell Ranger filtered matrix (Cell Ranger 10.1.0, GRCh38-2024-A).
- `data/<sample>/TAR/TAR_feature_bc_matrix/` — TAR-scRNA-seq feature-barcode matrix.
- `data/<sample>/TAR/TAR_reads.bed.gz` — per-sample TAR regions (steps 10–11).
- `data/annotations/` — GENCODE v47 GTF (`gencode.v47.annotation.gtf.gz`).
- `data/annotations/external/` — `LncBook.gtf.gz`, `NONCODE.bed.gz`, `LNCipedia.bed`, `FANTOM_lncRNA.bed` (GRCh38).
- `data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa` — uncompressed Ensembl GRCh38 FASTA (step 8; the `.fai` index is created if missing).

The TAR matrices are produced upstream by the TAR-scRNA-seq pipeline; that step is not part of this repository.

## Usage

Scripts are run from the repository root, in numeric order. The project directory defaults to the working directory and can be overridden with an environment variable:

```bash
export LNCRNA_IPF_PROJECT_DIR=/path/to/lncRNA-IPF-sc   # optional
Rscript scripts/01_import_and_qc.R
Rscript scripts/02_build_consensus_utar_catalog.R
# ...
Rscript scripts/15_known_vs_unannotated_specificity.R
```

Outputs are written under `data/` (tables in `data/TAR_catalog/`, `data/known_IPF_lncRNA/`, `data/single_cell_RNA_Habermann/`, `data/single_cell_TAR/`; figures in `data/plots/`).

## Requirements

R ≥ 4.3 with:

- CRAN: `Seurat` (v5; `JoinLayers` is used), `Matrix`, `dplyr`, `tidyr`, `tibble`, `readr`, `ggplot2`, `patchwork`, optionally `future`
- Bioconductor: `GenomicRanges`, `IRanges`, `S4Vectors`, `GenomeInfoDb`, `rtracklayer`, `Rsamtools`, `Biostrings`

```r
install.packages(c("Seurat", "Matrix", "dplyr", "tidyr", "tibble", "readr", "ggplot2", "patchwork", "future"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "IRanges", "S4Vectors", "GenomeInfoDb",
                       "rtracklayer", "Rsamtools", "Biostrings"))
```

## Known issues and limitations

- **Step 13 is missing.** The script that maps consensus TARs to single cells (`13_map_tar_to_single_cells.R`, originally `12_single_cell_TAR_compartment_mapping.R`) was not included in this version. Step 14 requires its output: `data/single_cell_TAR/IPFTAR_consensus_single_cell_matrices.rds`.
- **Incomplete VUILD53 TAR matrix.** About 94 % of VUILD53 TAR counts fall on chr18/20/21/22, and its TAR BED contains far more regions than its feature matrix. Most Priority A candidates lie on these chromosomes, so "present in 3/3 IPF" and "absent in control" results must be re-checked after regenerating VUILD53.
- **Control absence is mostly structural.** Many loci with 0 counts in the control matrix do have a TAR region in the control BED; "absent in control" should be interpreted with the refinement results (steps 10–11).
- **n = 3 IPF vs 1 control.** All IPF vs control comparisons are descriptive.
- **Pilot-level cell annotation.** The cluster-to-compartment mapping in step 12 is manual and specific to this pilot (resolution 0.10).

## Code-cleaning note

Scripts were renamed and cleaned for this repository: comments and console/plot text were translated to English and formatting was compacted. **The analysis logic is unchanged**; this was verified by comparing the R parse tokens of every original and cleaned script (only the `source("config/config.R")` line and the base-path assignment differ). The only change visible in outputs is in step 5, where the `direction` labels are now `"Higher in IPF"` / `"Lower in IPF"` (previously Spanish), and in plot titles/axis labels, which are now English.
