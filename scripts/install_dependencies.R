# install_dependencies.R
# One-shot installation of all R packages used in the BAL multi-omics
# analysis pipeline. Run once before executing any of the analysis scripts.
#
# Usage:  Rscript install_dependencies.R

# ---- CRAN packages ----
cran_pkgs <- c(
  "tidyverse", "patchwork", "ggrepel", "msigdbr",
  "readxl", "pROC", "effsize",
  "vegan", "future", "reticulate", "hdf5r", "devtools", "BiocManager",
  "Seurat", "harmony", "survival", "coxphf", "logistf", "ppcor", "DescTools",
  "jsonlite", "png", "magick"
)
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

# ---- Bioconductor packages ----
bioc_pkgs <- c(
  "DESeq2", "edgeR", "clusterProfiler",
  "AnnotationDbi", "org.Hs.eg.db",
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


# ---- Python backend of MOFA2 (scripts 04, 05) ----
# MOFA2 trains through the Python package mofapy2 (version used: 0.7.4). Install it in the Python environment given by MOFA_PYTHON:
#   python3 -m pip install mofapy2==0.7.4
# Versions of the environment in which the scripts were checked against the article: ../docs/package_versions.txt
