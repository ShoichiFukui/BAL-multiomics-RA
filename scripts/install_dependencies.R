# ============================================================================
# install_dependencies.R
# One-shot installation of all R packages used in the Hosogaya BAL multi-omics
# analysis pipeline. Run once before executing any of the analysis scripts.
#
# Usage:  Rscript install_dependencies.R
# ============================================================================

# ---- CRAN packages ----
cran_pkgs <- c(
  "tidyverse", "patchwork", "ggrepel", "pheatmap", "corrplot", "gridExtra",
  "openxlsx", "readxl", "pROC", "randomForest", "caret", "glmnet", "effsize",
  "vegan", "future", "reticulate", "hdf5r", "devtools", "BiocManager",
  "Seurat", "harmony"
)
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# ---- Bioconductor packages ----
bioc_pkgs <- c(
  "DESeq2", "apeglm", "edgeR", "clusterProfiler",
  "AnnotationDbi", "org.Hs.eg.db", "biomaRt",
  "GSVA", "SingleCellExperiment", "MOFA2", "fgsea"
)
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# ---- GitHub-only ----
# BayesPrism — Danko Lab (not on CRAN/Bioconductor as of writing)
if (!requireNamespace("BayesPrism", quietly = TRUE)) {
  devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
}

cat("\nAll dependencies installed.\n")
cat("Reminder: MOFA2 requires a working Python environment with the 'mofapy2'\n")
cat("package. See https://biofam.github.io/MOFA2/.\n")
