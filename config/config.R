# ============================================================
# config/config.R
#
# Project paths shared by all scripts. Every script is meant to be
# run from the repository root, e.g.:
#
#   Rscript scripts/01_import_and_qc.R
#
# Set the environment variable LNCRNA_IPF_PROJECT_DIR to the folder
# that contains `data/` (inputs and results). If it is not set, the
# repository root (current working directory) is used.
# ============================================================

PROJECT_DIR <- Sys.getenv("LNCRNA_IPF_PROJECT_DIR", unset = normalizePath("."))
DATA_DIR    <- file.path(PROJECT_DIR, "data")
