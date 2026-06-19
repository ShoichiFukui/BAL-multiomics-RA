# =====================================
# RA-ILD BAL Deconvolution v3
# Author Annotation + BayesPrism Deconvolution
# with FCM validation
# Version 3.0
# =====================================
#
# Pipeline overview:
# - Annotation: uses peer-reviewed cell type annotations from the original
#   publications of each public scRNA-seq dataset; datasets without annotation
#   are classified with BAL-specific canonical markers.
#   Macrophage subtypes: FABP4/MARCO (Alveolar) vs FCN1/CD14 (Monocyte-derived).
# - Deconvolution: BayesPrism (Bayesian model with built-in platform correction).
# - FCM validation: Macrophage/Lymphocyte/Neutrophil correlations computed automatically.
# - QC: canonical marker validation and reference quality checks.
#
# Required packages (install beforehand):
# install.packages("devtools")
# devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
# BiocManager::install("SingleCellExperiment")
# =====================================

# =====================================
# Initial setup
# =====================================

# Temporary directory (workaround when default /var/folders/ runs out of space)
tmp_dir <- file.path(getwd(), "tmp_R")
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = tmp_dir)
options(
  future.globals.maxSize = 16000 * 1024^2,
  bitmapType = "cairo"  # Quartz alternative (avoids Quartz temp file errors)
)

library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(BayesPrism)
library(biomaRt)
library(pheatmap)
library(harmony)
library(patchwork)
library(hdf5r)
library(readxl)
library(future)

plan("sequential")

# Working directory
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)

# Output directory
output_dir <- "./output_v3_BayesPrism"
dir.create(output_dir, showWarnings = FALSE)

# Data directory (reference annotation files, etc.)
data_dir <- "."

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  RA-ILD BAL Deconvolution v3                                ║\n")
cat("║  Author Annotation + BayesPrism Deconvolution               ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# =====================================
# PART 1: Load GSE145926 (COVID-19 BALF) data
# =====================================

cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 1: Load GSE145926 COVID-19 BALF data                  ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

covid_processed_path <- "./output_integrated/GSE145926_COVID_processed.rds"

if (file.exists(covid_processed_path)) {
  cat("Loading pre-processed GSE145926 data...\n")
  covid_seurat <- readRDS(covid_processed_path)
  cat(sprintf("✓ Loaded: %d cells\n", ncol(covid_seurat)))
  
  if (!"sample_id" %in% colnames(covid_seurat@meta.data)) {
    covid_seurat$sample_id <- covid_seurat$sample
  }
  
  # Preserve any pre-existing celltype (for later comparison with author annotation)
  if ("celltype" %in% colnames(covid_seurat@meta.data)) {
    covid_seurat$celltype_v2 <- covid_seurat$celltype
    covid_seurat$celltype <- NULL
  }
  
} else {
  cat("Processing GSE145926 from raw H5 files...\n\n")
  
  sample_info_covid <- data.frame(
    GSM = c("GSM4339769", "GSM4339770", "GSM4339772",
            "GSM4339771", "GSM4339773", "GSM4339774",
            "GSM4475051", "GSM4475052", "GSM4475053",
            "GSM4475048", "GSM4475049", "GSM4475050"),
    SampleName = c("M1", "M2", "M3",
                   "S1", "S2", "S3", "S4", "S5", "S6",
                   "HC1", "HC2", "HC3"),
    Disease = c(rep("Moderate_COVID", 3),
                rep("Severe_COVID", 6),
                rep("Healthy", 3)),
    stringsAsFactors = FALSE
  )
  
  covid_data_path <- "./GSE145926_RAW"
  h5_files <- list.files(covid_data_path, pattern = "\\.h5$", full.names = TRUE)
  cat(sprintf("Found %d H5 files\n\n", length(h5_files)))
  
  seurat_list_covid <- list()
  
  for (h5_file in h5_files) {
    filename <- basename(h5_file)
    gsm <- regmatches(filename, regexpr("GSM[0-9]+", filename))
    idx <- which(sample_info_covid$GSM == gsm)
    if (length(idx) == 0) { cat(sprintf("Unknown GSM: %s, skipping\n", gsm)); next }
    
    sample_name <- sample_info_covid$SampleName[idx]
    disease <- sample_info_covid$Disease[idx]
    cat(sprintf("Loading: %s (%s) - %s\n", sample_name, gsm, disease))
    
    counts <- Read10X_h5(h5_file)
    seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name,
                                     min.cells = 3, min.features = 200)
    seurat_obj$sample <- sample_name
    seurat_obj$sample_id <- sample_name
    seurat_obj$disease_status <- disease
    seurat_obj$condition <- "COVID"
    seurat_obj$dataset <- "GSE145926"
    seurat_list_covid[[sample_name]] <- seurat_obj
    cat(sprintf("  -> %d cells\n", ncol(seurat_obj)))
  }
  
  covid_seurat <- merge(seurat_list_covid[[1]], y = seurat_list_covid[-1],
                        add.cell.ids = names(seurat_list_covid))
  cat(sprintf("\n✓ GSE145926 loaded: %d cells\n", ncol(covid_seurat)))
  
  rm(seurat_list_covid); gc()
  
  # QC
  covid_seurat[["percent.mt"]] <- PercentageFeatureSet(covid_seurat, pattern = "^MT-")
  covid_seurat <- subset(covid_seurat,
                         subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  cat(sprintf("After QC: %d cells\n", ncol(covid_seurat)))
  
  # Preprocessing
  covid_seurat <- NormalizeData(covid_seurat, verbose = FALSE)
  covid_seurat <- FindVariableFeatures(covid_seurat, nfeatures = 3000, verbose = FALSE)
  covid_seurat <- ScaleData(covid_seurat, verbose = FALSE)
  covid_seurat <- RunPCA(covid_seurat, npcs = 50, verbose = FALSE)
  covid_seurat <- RunUMAP(covid_seurat, dims = 1:30, verbose = FALSE)
  covid_seurat <- FindNeighbors(covid_seurat, dims = 1:30, verbose = FALSE)
  covid_seurat <- FindClusters(covid_seurat, resolution = 0.5, verbose = FALSE)
  
  saveRDS(covid_seurat, covid_processed_path)
  cat("✓ GSE145926 saved\n")
}

gc()

# =====================================
# PART 2: GSE193782 (Healthy BAL)
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 2: GSE193782 (Healthy BAL)                            ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# --- Memory-efficient loading ---
# GSE193782 contains large assays (ADT/ITG/SCT). Use a pre-trimmed lite RDS if
# available; otherwise load the original file and trim it immediately.

healthy_lite_path <- "./output_v3_BayesPrism/GSE193782_lite.rds"

if (file.exists(healthy_lite_path)) {
  cat("Loading lightweight GSE193782...\n")
  healthy_seurat <- readRDS(healthy_lite_path)
  cat(sprintf("✓ Loaded lite version: %d cells\n", ncol(healthy_seurat)))
  
} else {
  cat("Loading GSE193782 (full object, may require >60GB)...\n")
  cat("  If this fails with memory error, try one of:\n")
  cat("    Option A: echo 'R_MAX_VSIZE=128Gb' >> ~/.Renviron  then restart R\n")
  cat("    Option B: R --max-vsize=128G -f this_script.R\n")
  cat("    Option C: Manually extract RNA counts from the RDS (see below)\n\n")
  
  # First try a direct read
  healthy_seurat <- tryCatch({
    obj <- readRDS("./analysis/GSE193782_SeuratObject.rds")
    cat(sprintf("Original: %d cells\n", ncol(obj)))

    # Immediately drop unneeded assays to free memory
    DefaultAssay(obj) <- "RNA"
    for (assay_name in c("ADT", "ITG", "SCT")) {
      if (assay_name %in% names(obj@assays)) {
        obj[[assay_name]] <- NULL
        cat(sprintf("  ✓ Removed %s assay\n", assay_name))
        gc()
      }
    }
    for (red_name in c("pca.SCT", "umap.SCT")) {
      if (red_name %in% names(obj@reductions)) {
        obj@reductions[[red_name]] <- NULL
      }
    }
    obj@commands <- list()
    gc()
    obj
    
  }, error = function(e) {
    cat(sprintf("\n⚠ Direct loading failed: %s\n", e$message))
    cat("  Attempting alternative extraction via connection streaming...\n\n")
    
    # Alternative: extract only the RNA counts from the RDS via a connection,
    # allowing GC during read
    con <- gzcon(file("./analysis/GSE193782_SeuratObject.rds", "rb"))
    obj <- tryCatch({
      readRDS(con)
    }, error = function(e2) {
      close(con)
      stop(paste0(
        "Out of memory loading GSE193782. Please do the following:\n",
        "  1. In a terminal: echo 'R_MAX_VSIZE=128Gb' >> ~/.Renviron\n",
        "  2. Restart R/RStudio\n",
        "  3. Re-run the script\n",
        "  Original error: ", e2$message
      ))
    })
    close(con)
    
    DefaultAssay(obj) <- "RNA"
    for (a in c("ADT", "ITG", "SCT")) {
      if (a %in% names(obj@assays)) { obj[[a]] <- NULL; gc() }
    }
    obj@commands <- list()
    gc()
    obj
  })
  
  # Subsample (done first to save memory)
  target_healthy <- 20000
  if (ncol(healthy_seurat) > target_healthy) {
    set.seed(123)
    keep_indices <- sample.int(ncol(healthy_seurat), target_healthy)
    healthy_seurat <- healthy_seurat[, keep_indices]
    cat(sprintf("✓ Subsampled to: %d cells\n", ncol(healthy_seurat)))
  }
  gc()
  
  # Save lite version (enables fast loading next time)
  cat("  Saving lightweight version for future use...\n")
  dir.create(dirname(healthy_lite_path), showWarnings = FALSE)
  saveRDS(healthy_seurat, healthy_lite_path)
  cat("  ✓ Saved to:", healthy_lite_path, "\n")
}

healthy_seurat$condition <- "Healthy"
healthy_seurat$sample_id <- ifelse("orig.ident" %in% colnames(healthy_seurat@meta.data),
                                   as.character(healthy_seurat$orig.ident), "Healthy")
healthy_seurat$dataset <- "GSE193782"
gc()

# =====================================
# PART 3: GSE184735 (Sarcoidosis BAL)
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 3: GSE184735 (Sarcoidosis BAL)                        ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

load_sarcoid <- function(data_path, max_per_sample = 1000) {
  all_files <- list.files(data_path, recursive = TRUE, full.names = TRUE)
  mtx_files <- all_files[grepl("matrix\\.mtx", all_files)]
  mtx_dirs <- unique(dirname(mtx_files))
  seurat_list <- list()
  
  for (i in seq_along(mtx_dirs)) {
    path_parts <- strsplit(mtx_dirs[i], "/")[[1]]
    sample_name <- path_parts[grep("^GSM", path_parts)][1]
    if (is.na(sample_name)) sample_name <- basename(dirname(mtx_dirs[i]))
    cat(sprintf("  [%d/%d] %s... ", i, length(mtx_dirs), sample_name))
    
    tryCatch({
      counts <- Read10X(mtx_dirs[i])
      if (ncol(counts) > max_per_sample) {
        set.seed(100 + i)
        counts <- counts[, sample.int(ncol(counts), max_per_sample)]
      }
      obj <- CreateSeuratObject(counts, project = sample_name, min.cells = 3, min.features = 200)
      obj$sample_id <- sample_name
      obj$condition <- "Sarcoidosis"
      obj$dataset <- "GSE184735"
      seurat_list[[sample_name]] <- obj
      cat(ncol(obj), "cells\n")
      rm(counts); gc(verbose = FALSE)
    }, error = function(e) cat("Error\n"))
  }
  return(seurat_list)
}

sarcoid_list <- load_sarcoid("./analysis/GSE184735/extracted", max_per_sample = 1000)

if (length(sarcoid_list) > 1) {
  sarcoid_merged <- merge(sarcoid_list[[1]], y = sarcoid_list[-1],
                          add.cell.ids = names(sarcoid_list))
} else {
  sarcoid_merged <- sarcoid_list[[1]]
}
rm(sarcoid_list); gc()
cat(sprintf("\n✓ Sarcoidosis: %d cells\n", ncol(sarcoid_merged)))

# =====================================
# PART 4: Merge 3 datasets
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 4: Merging 3 Datasets                                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Healthy (GSE193782):     %d cells\n", ncol(healthy_seurat)))
cat(sprintf("  Sarcoidosis (GSE184735): %d cells\n", ncol(sarcoid_merged)))
cat(sprintf("  COVID-19 (GSE145926):    %d cells\n", ncol(covid_seurat)))

common_genes <- Reduce(intersect, list(
  rownames(healthy_seurat), rownames(sarcoid_merged), rownames(covid_seurat)
))
cat(sprintf("\nCommon genes: %d\n", length(common_genes)))

healthy_seurat <- subset(healthy_seurat, features = common_genes)
sarcoid_merged <- subset(sarcoid_merged, features = common_genes)
covid_seurat   <- subset(covid_seurat, features = common_genes)

reference_obj <- merge(
  x = healthy_seurat,
  y = list(sarcoid_merged, covid_seurat),
  add.cell.ids = c("Healthy", "Sarcoid", "COVID"),
  project = "BAL_Reference"
)

rm(healthy_seurat, sarcoid_merged, covid_seurat); gc()
cat(sprintf("\n✓ Total: %d cells × %d genes\n", ncol(reference_obj), nrow(reference_obj)))
print(table(reference_obj$dataset))

# =====================================
# PART 5: QC filtering
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 5: Quality Control                                    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

reference_obj[["percent.mt"]] <- PercentageFeatureSet(reference_obj, pattern = "^MT-")
reference_obj[["percent.ribo"]] <- PercentageFeatureSet(reference_obj, pattern = "^RP[SL]")

cat(sprintf("Before QC: %d cells\n", ncol(reference_obj)))

reference_obj <- subset(reference_obj,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                          nCount_RNA > 500 & percent.mt < 20)

cat(sprintf("After QC: %d cells\n", ncol(reference_obj)))
gc()

# =====================================
# PART 6: Normalization, dimensionality reduction, Harmony integration
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 6: Normalization, PCA & Harmony Integration           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

reference_obj <- NormalizeData(reference_obj, verbose = FALSE)
cat("✓ Normalization\n")

reference_obj <- FindVariableFeatures(reference_obj, nfeatures = 2000, verbose = FALSE)
cat("✓ Variable features\n")

reference_obj <- ScaleData(reference_obj, features = VariableFeatures(reference_obj), verbose = FALSE)
cat("✓ Scaling\n")

reference_obj <- RunPCA(reference_obj, npcs = 50, verbose = FALSE)
cat("✓ PCA\n")

# Harmony batch correction (for visualization; not used for deconvolution)
reference_obj <- RunHarmony(reference_obj, "dataset", verbose = TRUE)
cat("✓ Harmony completed\n")

reference_obj <- RunUMAP(reference_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
reference_obj <- FindNeighbors(reference_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
reference_obj <- FindClusters(reference_obj, resolution = 0.8, verbose = FALSE)
cat("✓ Clustering & UMAP (Harmony-corrected)\n")
cat(sprintf("  Clusters: %d\n", length(unique(reference_obj$seurat_clusters))))

# Join layers (Seurat v5)
reference_obj <- JoinLayers(reference_obj)
cat("✓ Layers joined\n")
gc()

# =====================================
# PART 7: Author annotation
# =====================================
#
# Strategy: use the peer-reviewed cell type annotations from each public
#   dataset's original publication directly
#   ("Cell type labels from the original publications were used"),
#   avoiding reference-mismatch criticism.
#
# Fallback (only when author annotation is absent from metadata):
#   scoring-based classification using BAL canonical marker genes.
#
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 7: Cell Type Annotation                                ║\n")
cat("║  Author peer-reviewed annotation + canonical marker check    ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# ===== Remove leftover columns from any prior pipeline =====
# Columns carried over from a previous run can be falsely matched during
# annotation lookup, so remove them completely here.
v2_remove_cols <- c("celltype", "celltype_v2", "celltype_fine",
                    "Neutrophil_score1", "Alveolar_Mac_score1",
                    "Monocyte_Mac_score1", "T_cell_score1",
                    "NK_score1", "B_cell_score1", "Epithelial_score1")

cat("Removing leftover pipeline artifact columns...\n")
for (col_name in v2_remove_cols) {
  if (col_name %in% colnames(reference_obj@meta.data)) {
    reference_obj@meta.data[[col_name]] <- NULL
    cat(sprintf("  ✓ Removed: %s\n", col_name))
  }
}

# Initialize fresh annotation columns (start from NA)
reference_obj$celltype <- NA_character_
reference_obj$annotation_source <- NA_character_
reference_obj$original_author_label <- NA_character_
cat("  ✓ Initialized clean annotation columns\n\n")

# --- 7a: Define canonical BAL marker genes & harmonization function ---
# Define functions/constants needed before applying author annotations

cat("=== Preparatory: Define canonical BAL marker genes & harmonization function ===\n\n")

bal_markers <- list(
  Macrophage  = c("CD68", "CD163", "MARCO", "FABP4", "MSR1", "MRC1", "APOE", "PPARG"),
  T_cell      = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "TRBC2"),
  NK          = c("NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "FCGR3A"),
  B_cell      = c("CD79A", "CD79B", "MS4A1", "CD19", "PAX5"),
  Plasma      = c("JCHAIN", "MZB1", "IGHG1", "IGHG2", "SDC1", "XBP1"),
  Neutrophil  = c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9", "MMP9"),
  DC          = c("FCER1A", "CD1C", "CLEC9A", "IRF7", "IRF8", "LILRA4"),
  Epithelial  = c("EPCAM", "KRT18", "KRT19", "SFTPC", "SFTPB", "SCGB1A1"),
  Mast        = c("TPSAB1", "TPSB2", "CPA3", "KIT", "HDC")
)

# Macrophage subtype markers
mac_subtype_markers <- list(
  Alveolar_Mac    = c("FABP4", "MARCO", "MRC1", "PPARG", "MCEMP1"),
  Monocyte_Mac    = c("FCN1", "CD14", "S100A12", "VCAN", "CCR2")
)

cat("Canonical BAL cell type markers defined:\n")
for (ct in names(bal_markers)) {
  avail <- intersect(bal_markers[[ct]], rownames(reference_obj))
  cat(sprintf("  %-15s: %d/%d markers available\n", ct, length(avail), length(bal_markers[[ct]])))
}
cat("\n")

# --- harmonize_celltype: label harmonization function ---
# Map each dataset's author labels to unified labels for BAL deconvolution

harmonize_celltype <- function(label) {
  label <- as.character(label)

  result <- case_when(
    # === GSE193782 Cell.Types specific labels (highest priority) ===
    # AMs = Alveolar Macrophages (major BAL cell population)
    grepl("^AMs$|^AMs\\.", label) ~ "Macrophage",
    # FOLR2.IMs = FOLR2+ Interstitial Macrophages
    grepl("FOLR2\\.IMs|^IMs$", label) ~ "Macrophage",
    # Cyc.Mye = Cycling Myeloid cells (proliferating macrophage lineage)
    grepl("Cyc\\.Mye", label) ~ "Macrophage",
    # Cyc.Lym = Cycling Lymphocytes
    grepl("Cyc\\.Lym", label) ~ "T_cell",
    # Lym = Lymphocytes (T cell predominant in BAL)
    grepl("^Lym$", label) ~ "T_cell",
    # Mig.DCs = Migratory DCs
    grepl("Mig\\.DCs|^DC1$|^DC2$", label) ~ "DC",
    # Epi = Epithelial
    grepl("^Epi$", label) ~ "Epithelial",
    # Mono = Monocytes
    grepl("^Mono$|^Mono\\.", label) ~ "Macrophage",

    # === General labels (GSE145926, GSE184735 marker labels, etc.) ===
    # Macrophage lineage (most important in BAL)
    grepl("Macrophage|macrophage|Macro|MΦ|Alveolar.Mac|Alveolar_Mac|Mono.Mac|Monocyte_Mac|AM$|MDM$|BALF.Mac",
          label, ignore.case = TRUE) ~ "Macrophage",
    grepl("^Monocyte|^monocyte|^CD14.Mono|^CD16.Mono|^cMono|^ncMono|^intMono",
          label, ignore.case = TRUE) ~ "Macrophage",
    
    # DC
    grepl("^DC$|^DC[0-9]|Dendritic|^pDC|^cDC|^mDC|plasmacytoid.DC", label, ignore.case = TRUE) ~ "DC",
    
    # T cells
    grepl("T.cell|T_cell|Tcell|^CD4|^CD8|Treg|^T$|CTL|Th1|Th2|Th17|MAIT|gdT|gamma.delta",
          label, ignore.case = TRUE) ~ "T_cell",
    
    # NK
    grepl("^NK|Natural.killer|natural.killer|NK.cell", label, ignore.case = TRUE) ~ "NK",
    
    # B cells
    grepl("^B.cell|^B_cell|^B$|Bcell|Naive.B|Memory.B|^Bn$|^Bm$",
          label, ignore.case = TRUE) ~ "B_cell",
    
    # Plasma
    grepl("Plasma|^PC$|plasmablast", label, ignore.case = TRUE) ~ "Plasma",
    
    # Neutrophil
    grepl("Neutrophil|neutrophil|^Neut$|PMN", label, ignore.case = TRUE) ~ "Neutrophil",
    
    # Epithelial
    grepl("Epithelial|epithelial|Ciliated|Secretory|Club|AT1|AT2|Basal|^Epi$",
          label, ignore.case = TRUE) ~ "Epithelial",
    
    # Mast
    grepl("Mast|mast", label, ignore.case = TRUE) ~ "Mast",
    
    # Granulocyte (other)
    grepl("Eosinophil|Basophil", label, ignore.case = TRUE) ~ "Granulocyte",
    
    # Unknown / Other
    TRUE ~ "Other"
  )
  
  return(result)
}


# --- 7b: Load and apply original-author annotations ---
cat("=== Step 1: Loading author annotations ===\n\n")

# =====================================================================
# GSE145926: author annotation from Liao et al. Nature Medicine 2020
# Source: https://github.com/zhangzlab/covid_balf/
#   - all.cell.annotation.meta.txt (all cell types)
#   - myeloid.cell.annotation.meta.txt (macrophage subtypes)
# =====================================================================

liao_annotation_path <- file.path(data_dir, "reference", "GSE145926_all_cell_annotation.txt")
liao_myeloid_path <- file.path(data_dir, "reference", "GSE145926_myeloid_annotation.txt")

if (file.exists(liao_annotation_path)) {
  cat("--- GSE145926: Loading Liao et al. author annotations ---\n")
  
  liao_meta <- read.delim(liao_annotation_path, sep = "\t", header = TRUE)
  cat(sprintf("  Loaded %d cells from annotation file\n", nrow(liao_meta)))
  cat("  Author cell types:\n")
  print(table(liao_meta$celltype))
  
  # Build barcode matching key
  # Liao: ID = "AAACCTGAGACACTAA_1", sample_new = "HC1"
  #   -> key = "HC1_AAACCTGAGACACTAA"
  liao_barcode <- gsub("_[0-9]+$", "", liao_meta$ID)  # strip "_1" suffix
  liao_key <- paste0(liao_meta$sample_new, "_", liao_barcode)
  liao_meta$match_key <- liao_key

  # RDS barcode = "COVID_M1_AAACCTGAGATGTCGG-1"
  #   -> key = "M1_AAACCTGAGATGTCGG"
  covid_cells <- which(reference_obj$dataset == "GSE145926")
  covid_barcodes <- colnames(reference_obj)[covid_cells]
  
  # "COVID_{sample}_{barcode}-{N}" -> "{sample}_{barcode}"
  rds_key <- gsub("^COVID_", "", covid_barcodes)           # "M1_AAACCTGAGATGTCGG-1"
  rds_key <- gsub("-[0-9]+$", "", rds_key)                 # "M1_AAACCTGAGATGTCGG"

  # Run matching
  match_idx <- match(rds_key, liao_meta$match_key)
  n_matched <- sum(!is.na(match_idx))
  cat(sprintf("\n  Barcode matching: %d / %d matched (%.1f%%)\n",
              n_matched, length(covid_cells), n_matched / length(covid_cells) * 100))
  
  if (n_matched > 0) {
    # Apply Liao et al. annotation to matched cells
    matched_mask <- !is.na(match_idx)
    matched_labels <- liao_meta$celltype[match_idx[matched_mask]]

    reference_obj$original_author_label[covid_cells[matched_mask]] <- as.character(matched_labels)
    reference_obj$annotation_source[covid_cells[matched_mask]] <- "author_Liao2020"

    # Map Liao et al. labels -> unified BAL labels
    # Liao labels: Macrophages, T, NK, B, Neutrophil, mDC, pDC, Epithelial, Mast, Plasma
    liao_harmonize <- function(label) {
      case_when(
        grepl("Macrophage", label) ~ "Macrophage",
        grepl("^T$", label)       ~ "T_cell",
        grepl("^NK$", label)      ~ "NK",
        grepl("^B$", label)       ~ "B_cell",
        grepl("Neutrophil", label) ~ "Neutrophil",
        grepl("mDC|pDC", label)   ~ "DC",
        grepl("Epithelial", label) ~ "Epithelial",
        grepl("Mast", label)      ~ "Mast",
        grepl("Plasma", label)    ~ "Plasma",
        TRUE                      ~ "Other"
      )
    }
    
    harmonized <- liao_harmonize(matched_labels)
    reference_obj$celltype[covid_cells[matched_mask]] <- harmonized
    
    cat("  Liao et al. → Harmonized mapping:\n")
    for (orig in sort(unique(as.character(matched_labels)))) {
      mapped <- liao_harmonize(orig)
      n <- sum(matched_labels == orig)
      cat(sprintf("    %-20s → %-15s (n=%d)\n", orig, mapped, n))
    }
    
    # Handle unmatched cells
    # Likely QC-filtered by Liao et al., so exclude from the reference
    # (applying the original authors' quality standards, not cherry-picking)
    unmatched_cells <- covid_cells[!matched_mask]
    if (length(unmatched_cells) > 0) {
      cat(sprintf("\n  Excluding %d unmatched GSE145926 cells (not in Liao et al. published annotation)\n",
                  length(unmatched_cells)))
      cat("  Rationale: These cells were likely QC-filtered by the original authors.\n")
      cat("  Applying the authors' quality standards is standard practice.\n")
      
      cat(sprintf("  Before exclusion: %d cells\n", ncol(reference_obj)))
      
      cells_to_keep <- setdiff(colnames(reference_obj), colnames(reference_obj)[unmatched_cells])
      reference_obj <- subset(reference_obj, cells = cells_to_keep)
      gc()
      
      cat(sprintf("  After exclusion: %d cells\n", ncol(reference_obj)))
    }
  }
  
  # Also load myeloid subtypes (for macrophage subclassification)
  if (file.exists(liao_myeloid_path)) {
    liao_myeloid <- read.delim(liao_myeloid_path, sep = "\t", header = TRUE)
    cat(sprintf("\n  Myeloid subtype file: %d cells loaded\n", nrow(liao_myeloid)))
    cat("  Myeloid groups:\n")
    print(table(liao_myeloid$celltype))
  }
  
} else {
  cat("⚠ GSE145926 annotation file not found at:", liao_annotation_path, "\n")
  cat("  Download from: https://github.com/zhangzlab/covid_balf/\n")
  cat("  Will use kNN transfer or canonical markers instead.\n")
}

cat("\n")

# =====================================================================
# GSE193782: "Cell.Types" column in RDS metadata (author annotation)
# GSE184735: no author annotation -> kNN transfer
# =====================================================================

cat("--- GSE193782 & GSE184735: Checking RDS metadata ---\n")

# GSE193782
gse193782_cells <- which(reference_obj$dataset == "GSE193782")
if ("Cell.Types" %in% colnames(reference_obj@meta.data)) {
  ct_labels <- reference_obj@meta.data[gse193782_cells, "Cell.Types"]
  n_valid <- sum(!is.na(ct_labels))
  if (n_valid > 0) {
    cat(sprintf("  GSE193782: Found 'Cell.Types' (%d cells with labels)\n", n_valid))
    cat("  Labels: ", paste(sort(unique(ct_labels[!is.na(ct_labels)])), collapse = ", "), "\n")
    
    reference_obj$original_author_label[gse193782_cells] <- as.character(ct_labels)
    harmonized_193782 <- harmonize_celltype(ct_labels)
    reference_obj$celltype[gse193782_cells] <- harmonized_193782
    reference_obj$annotation_source[gse193782_cells] <- "author_GSE193782"
    
    cat("  GSE193782 → Harmonized mapping:\n")
    for (orig in sort(unique(as.character(ct_labels[!is.na(ct_labels)])))) {
      mapped <- harmonize_celltype(orig)
      n <- sum(ct_labels == orig, na.rm = TRUE)
      cat(sprintf("    %-30s → %-15s (n=%d)\n", orig, mapped, n))
    }
  } else {
    cat("  GSE193782: Cell.Types column exists but all NA\n")
  }
} else {
  cat("  GSE193782: 'Cell.Types' column not found in metadata\n")
}

# GSE184735
gse184735_cells <- which(reference_obj$dataset == "GSE184735")
cat(sprintf("\n  GSE184735: %d cells → will use kNN transfer\n", length(gse184735_cells)))

# Tally unannotated cells
n_annotated <- sum(!is.na(reference_obj$celltype))
n_unannotated <- sum(is.na(reference_obj$celltype))
cat(sprintf("\n  Summary: %d annotated, %d unannotated (to be resolved by kNN)\n\n",
            n_annotated, n_unannotated))


# --- 7c: kNN label transfer for unannotated cells ---
cat("=== Step 2: kNN label transfer for unannotated cells ===\n\n")

# Process cells that did not receive an author annotation
# (all of GSE184735 + unmatched cells from GSE145926/193782)

unannotated_cells <- which(is.na(reference_obj$celltype))
annotated_cells <- which(!is.na(reference_obj$celltype))

cat(sprintf("Annotated cells: %d\n", length(annotated_cells)))
cat(sprintf("Unannotated cells: %d\n", length(unannotated_cells)))

if (length(unannotated_cells) > 0 && length(annotated_cells) > 100) {
  
  if ("harmony" %in% names(reference_obj@reductions)) {
    cat("\nUsing kNN label transfer in Harmony space...\n")
    
    harmony_embeddings <- Embeddings(reference_obj, "harmony")
    n_dims <- min(30, ncol(harmony_embeddings))
    harmony_embeddings <- harmony_embeddings[, 1:n_dims]
    
    k <- 20
    
    ref_embeddings <- harmony_embeddings[annotated_cells, , drop = FALSE]
    ref_labels <- reference_obj$celltype[annotated_cells]
    query_embeddings <- harmony_embeddings[unannotated_cells, , drop = FALSE]
    
    # Batch processing
    batch_size <- 5000
    n_batches <- ceiling(length(unannotated_cells) / batch_size)
    knn_labels <- character(length(unannotated_cells))
    knn_confidence <- numeric(length(unannotated_cells))

    # Precompute on the reference side
    b2 <- rowSums(ref_embeddings^2)
    
    for (b in 1:n_batches) {
      start_idx <- (b - 1) * batch_size + 1
      end_idx <- min(b * batch_size, length(unannotated_cells))
      batch_emb <- query_embeddings[start_idx:end_idx, , drop = FALSE]
      
      a2 <- rowSums(batch_emb^2)
      dist_matrix <- outer(a2, b2, "+") - 2 * tcrossprod(batch_emb, ref_embeddings)
      
      for (i in 1:(end_idx - start_idx + 1)) {
        knn_idx <- order(dist_matrix[i, ])[1:min(k, ncol(dist_matrix))]
        knn_votes <- table(ref_labels[knn_idx])
        winner <- names(knn_votes)[which.max(knn_votes)]
        confidence <- max(knn_votes) / sum(knn_votes)
        
        knn_labels[start_idx + i - 1] <- winner
        knn_confidence[start_idx + i - 1] <- confidence
      }
      
      cat(sprintf("  Batch %d/%d complete\n", b, n_batches))
    }
    
    reference_obj$celltype[unannotated_cells] <- knn_labels
    reference_obj$annotation_source[unannotated_cells] <- "kNN_transfer"
    
    cat(sprintf("\nkNN transfer results (mean confidence: %.2f):\n", mean(knn_confidence)))
    print(table(knn_labels))
    
    low_conf <- sum(knn_confidence < 0.5)
    if (low_conf > 0) {
      cat(sprintf("⚠ Low confidence cells (<50%%): %d (%.1f%%)\n",
                  low_conf, low_conf / length(unannotated_cells) * 100))
    }
    
  } else {
    # Harmony not run -> canonical marker fallback
    cat("Harmony not available → using canonical marker classification...\n")
    
    expr_data <- GetAssayData(reference_obj, layer = "data")
    ds_expr <- expr_data[, unannotated_cells, drop = FALSE]
    
    ct_scores <- matrix(0, nrow = length(unannotated_cells), ncol = length(bal_markers))
    colnames(ct_scores) <- names(bal_markers)
    
    for (ct in names(bal_markers)) {
      avail_markers <- intersect(bal_markers[[ct]], rownames(ds_expr))
      if (length(avail_markers) >= 2) {
        ct_scores[, ct] <- colMeans(ds_expr[avail_markers, , drop = FALSE])
      } else if (length(avail_markers) == 1) {
        ct_scores[, ct] <- ds_expr[avail_markers, ]
      }
    }
    
    max_scores <- apply(ct_scores, 1, max)
    assigned <- colnames(ct_scores)[apply(ct_scores, 1, which.max)]
    assigned[max_scores < 0.1] <- "Other"
    
    reference_obj$celltype[unannotated_cells] <- assigned
    reference_obj$annotation_source[unannotated_cells] <- "marker"
    
    cat("Marker-based classification results:\n")
    print(table(assigned))
  }
  
} else if (length(unannotated_cells) == 0) {
  cat("✓ All cells have author annotations. No kNN transfer needed.\n")
} else {
  cat("⚠ Insufficient annotated cells for kNN transfer.\n")
}

# Final annotation summary
cat("\n=== Final annotation summary ===\n")
cat("Cell type distribution:\n")
print(table(reference_obj$celltype))
cat("\nAnnotation source:\n")
print(table(reference_obj$annotation_source))
cat("\nBy dataset:\n")
for (ds in c("GSE145926", "GSE193782", "GSE184735")) {
  ds_cells <- which(reference_obj$dataset == ds)
  cat(sprintf("  %s: %s\n", ds,
              paste(names(table(reference_obj$annotation_source[ds_cells])),
                    table(reference_obj$annotation_source[ds_cells]),
                    sep = "=", collapse = ", ")))
}

# Handle any remaining NA
n_na <- sum(is.na(reference_obj$celltype))
if (n_na > 0) {
  cat(sprintf("\n⚠ %d cells with NA annotation → assigning 'Other'\n", n_na))
  reference_obj$celltype[is.na(reference_obj$celltype)] <- "Other"
}

# --- 7d: Macrophage subtype classification ---
cat("\n=== Step 3: Macrophage subtype classification ===\n")
cat("  Alveolar Mac (FABP4/MARCO) vs Monocyte-derived Mac (FCN1/CD14)\n\n")

mac_cells <- which(reference_obj$celltype == "Macrophage")
cat(sprintf("  Total Macrophages: %d\n", length(mac_cells)))

# Initialize celltype_fine (same as the coarse celltype)
reference_obj$celltype_fine <- reference_obj$celltype

if (length(mac_cells) > 0) {
  expr_data <- GetAssayData(reference_obj, layer = "data")
  
  avail_alv <- intersect(mac_subtype_markers$Alveolar_Mac, rownames(reference_obj))
  avail_mono <- intersect(mac_subtype_markers$Monocyte_Mac, rownames(reference_obj))
  
  cat(sprintf("  Alveolar markers available: %s\n", paste(avail_alv, collapse = ", ")))
  cat(sprintf("  Monocyte markers available: %s\n", paste(avail_mono, collapse = ", ")))
  
  if (length(avail_alv) >= 2 & length(avail_mono) >= 2) {
    alv_score <- colMeans(expr_data[avail_alv, mac_cells, drop = FALSE])
    mono_score <- colMeans(expr_data[avail_mono, mac_cells, drop = FALSE])
    
    mac_subtype <- ifelse(alv_score > mono_score, "Alveolar_Mac", "Monocyte_Mac")
    
    # GSE193782: cells with Cell.Types == "AMs" are Alveolar_Mac by definition
    # Override using original_author_label
    for (i in seq_along(mac_cells)) {
      ci <- mac_cells[i]
      orig <- reference_obj$original_author_label[ci]
      if (!is.na(orig)) {
        if (grepl("^AMs$|^AMs\\.", orig)) {
          mac_subtype[i] <- "Alveolar_Mac"
        } else if (grepl("^Mono$|^Mono\\.|FOLR2", orig)) {
          mac_subtype[i] <- "Monocyte_Mac"
        }
        # Also use GSE145926 pre-existing labels
        if (grepl("Alveolar_Mac", orig)) {
          mac_subtype[i] <- "Alveolar_Mac"
        } else if (grepl("Monocyte_Mac", orig)) {
          mac_subtype[i] <- "Monocyte_Mac"
        }
      }
    }
    
    reference_obj$celltype_fine[mac_cells] <- mac_subtype
    
    cat(sprintf("  ✓ Alveolar Mac: %d (%.1f%%)\n",
                sum(mac_subtype == "Alveolar_Mac"),
                sum(mac_subtype == "Alveolar_Mac") / length(mac_cells) * 100))
    cat(sprintf("  ✓ Monocyte Mac: %d (%.1f%%)\n",
                sum(mac_subtype == "Monocyte_Mac"),
                sum(mac_subtype == "Monocyte_Mac") / length(mac_cells) * 100))
  } else {
    cat("  ⚠ Insufficient markers for subclassification\n")
  }
}

# --- 7e: Annotation quality validation with canonical markers ---
cat("\n=== Step 4: Annotation quality validation with canonical markers ===\n\n")

# Check representative marker expression for each cell type
expr_data <- GetAssayData(reference_obj, layer = "data")

validation_markers <- list(
  Macrophage  = c("CD68", "MARCO", "FABP4"),
  T_cell      = c("CD3D", "CD3E"),
  NK          = c("NKG7", "GNLY"),
  B_cell      = c("CD79A", "MS4A1"),
  Neutrophil  = c("FCGR3B", "CSF3R"),
  DC          = c("FCER1A", "CD1C"),
  Epithelial  = c("EPCAM", "KRT18"),
  Plasma      = c("JCHAIN", "MZB1")
)

cat("Marker expression validation (mean log-normalized expression):\n")
cat(sprintf("%-15s %-15s %10s %10s %10s\n", "Cell Type", "Marker", "In-type", "Other", "Fold"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (ct in names(validation_markers)) {
  ct_cells <- which(reference_obj$celltype == ct)
  other_cells <- which(reference_obj$celltype != ct)
  
  if (length(ct_cells) == 0) next
  
  for (marker in validation_markers[[ct]]) {
    if (marker %in% rownames(expr_data)) {
      in_type_mean <- mean(expr_data[marker, ct_cells])
      other_mean <- mean(expr_data[marker, other_cells])
      enrichment <- ifelse(other_mean > 0, in_type_mean / other_mean, Inf)
      status <- ifelse(enrichment > 1.5, "✓", ifelse(enrichment > 1.0, "~", "⚠"))
      cat(sprintf("%-15s %-15s %10.2f %10.2f %9.1fx %s\n",
                  ct, marker, in_type_mean, other_mean, enrichment, status))
    }
  }
}

# --- 7f: Annotation summary ---
cat("\n\n=== Cell Type Distribution (Final) ===\n")
print(table(reference_obj$celltype))

cat("\n=== Annotation Source ===\n")
print(table(reference_obj$annotation_source))

cat("\n=== Dataset × Cell Type ===\n")
print(table(reference_obj$dataset, reference_obj$celltype))

cat("\n=== Fine Cell Types ===\n")
print(table(reference_obj$celltype_fine))

# --- 7g: UMAP visualization ---
cat("\nGenerating annotation plots...\n")

p1 <- DimPlot(reference_obj, group.by = "dataset") + ggtitle("Datasets (Harmony)")
p2 <- DimPlot(reference_obj, group.by = "celltype", label = TRUE, repel = TRUE) +
  ggtitle("Cell Types (Author Annotation)")
p3 <- DimPlot(reference_obj, group.by = "celltype_fine", label = TRUE, repel = TRUE) +
  ggtitle("Fine Cell Types")

ggsave(file.path(output_dir, "UMAP_datasets.pdf"), p1, width = 10, height = 8)
ggsave(file.path(output_dir, "UMAP_celltypes_author.pdf"), p2, width = 12, height = 10)
ggsave(file.path(output_dir, "UMAP_celltypes_fine.pdf"), p3, width = 14, height = 10)

# Annotation source plot
p4 <- DimPlot(reference_obj, group.by = "annotation_source") +
  ggtitle("Annotation Source (author vs marker)")
ggsave(file.path(output_dir, "UMAP_annotation_source.pdf"), p4, width = 10, height = 8)

# Canonical marker dot plot (annotation quality evidence)
marker_genes_for_plot <- c(
  "CD68", "MARCO", "FABP4", "FCN1",        # Macrophage
  "CD3D", "CD3E",                            # T cell
  "NKG7", "GNLY", "KLRD1",                  # NK
  "CD79A", "MS4A1",                          # B cell
  "JCHAIN", "MZB1",                          # Plasma
  "FCGR3B", "CSF3R",                         # Neutrophil
  "FCER1A", "CD1C",                          # DC
  "EPCAM", "KRT18"                           # Epithelial
)
marker_genes_for_plot <- intersect(marker_genes_for_plot, rownames(reference_obj))

if (length(marker_genes_for_plot) > 0) {
  p_dot <- DotPlot(reference_obj, features = marker_genes_for_plot, group.by = "celltype") +
    RotatedAxis() +
    ggtitle("Canonical Marker Validation") +
    theme(axis.text.x = element_text(size = 8))
  ggsave(file.path(output_dir, "canonical_marker_dotplot.pdf"), p_dot, width = 16, height = 8)
  cat("✓ Canonical marker dot plot saved\n")
}

# Save reference
saveRDS(reference_obj, file.path(output_dir, "BAL_reference_author_annotated.rds"))
cat("✓ Reference with author annotation saved\n")
gc()

# =====================================
# PART 8: Bulk data preparation
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 8: Bulk Data Preparation                              ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

bulk_path <- "./Hosogaya RNAseq/RawCount_heatmap_table.csv"
bulk_raw <- read.csv(bulk_path, row.names = 1)

if ("Cluster_ID" %in% colnames(bulk_raw)) {
  bulk_raw <- bulk_raw[, -which(colnames(bulk_raw) == "Cluster_ID")]
}

bulk_matrix <- as.matrix(bulk_raw)
cat(sprintf("Bulk data: %d genes × %d samples\n", nrow(bulk_matrix), ncol(bulk_matrix)))

# Gene name conversion (Ensembl -> Gene Symbol)
# Prefer org.Hs.eg.db (local); fall back to biomaRt only on failure
cat("\nConverting Ensembl IDs to Gene Symbols...\n")

gene_mapping <- tryCatch({
  # Method 1: org.Hs.eg.db (offline, fast)
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  ensembl_ids <- gsub("\\..*", "", rownames(bulk_matrix))  # remove version suffix
  mapped <- AnnotationDbi::select(org.Hs.eg.db, keys=ensembl_ids,
                                   columns="SYMBOL", keytype="ENSEMBL")
  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
  mapped <- mapped[!duplicated(mapped$ENSEMBL), ]  # keep first match
  data.frame(ensembl_gene_id=mapped$ENSEMBL, external_gene_name=mapped$SYMBOL)
}, error=function(e) {
  # Method 2: biomaRt (online fallback)
  cat("  org.Hs.eg.db failed, trying biomaRt...\n")
  for(mirror in c("useast","asia","www")) {
    result <- tryCatch({
      ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror=mirror)
      getBM(attributes=c("ensembl_gene_id","external_gene_name"),
            filters="ensembl_gene_id", values=rownames(bulk_matrix), mart=ensembl)
    }, error=function(e2) NULL)
    if(!is.null(result)) return(result)
  }
  stop("All Ensembl mirrors failed and org.Hs.eg.db unavailable")
})

cat(sprintf("  Mapped %d / %d Ensembl IDs to gene symbols\n", nrow(gene_mapping), nrow(bulk_matrix)))

ensembl_ids_clean <- gsub("\\..*", "", rownames(bulk_matrix))
matched <- match(ensembl_ids_clean, gene_mapping$ensembl_gene_id)
new_names <- gene_mapping$external_gene_name[matched]
valid <- !is.na(new_names) & new_names != ""

bulk_matrix <- bulk_matrix[valid, ]
rownames(bulk_matrix) <- new_names[valid]

# Duplicate genes: sum expression counts
if (any(duplicated(rownames(bulk_matrix)))) {
  cat("  Aggregating duplicated gene names...\n")
  bulk_df <- as.data.frame(bulk_matrix)
  bulk_df$gene <- rownames(bulk_matrix)
  bulk_agg <- aggregate(. ~ gene, data = bulk_df, FUN = sum)
  rownames(bulk_agg) <- bulk_agg$gene
  bulk_agg$gene <- NULL
  bulk_matrix <- as.matrix(bulk_agg)
}

cat(sprintf("✓ Bulk matrix: %d genes × %d samples\n", nrow(bulk_matrix), ncol(bulk_matrix)))

# =====================================
# PART 9: BayesPrism Deconvolution
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 9: BayesPrism Deconvolution                           ║\n")
cat("║  - Models platform differences explicitly                   ║\n")
cat("║  - Accounts for uncertainty via Bayesian estimation         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# --- 9a: Prepare BayesPrism input matrices ---
cat("Preparing BayesPrism input matrices...\n")

# scRNA-seq: cell x gene raw count matrix (required by BayesPrism)
sc_counts <- GetAssayData(reference_obj, layer = "counts")

# sparse -> dense conversion & transpose (BayesPrism expects cell x gene)
sc_dat <- t(as.matrix(sc_counts))
cat(sprintf("  sc.dat: %d cells × %d genes\n", nrow(sc_dat), ncol(sc_dat)))

# Bulk: sample x gene count matrix (required by BayesPrism)
bk_dat <- t(bulk_matrix)
cat(sprintf("  bk.dat: %d samples × %d genes\n", nrow(bk_dat), ncol(bk_dat)))

# Cell type labels (BayesPrism supports a 2-layer celltype + cell state structure)
cell_type_labels <- reference_obj$celltype       # coarse type (used for deconvolution)
cell_state_labels <- reference_obj$celltype_fine  # fine type (informational)

cat(sprintf("  Cell types: %s\n", paste(sort(unique(cell_type_labels)), collapse = ", ")))
cat(sprintf("  Cell states: %s\n", paste(sort(unique(cell_state_labels)), collapse = ", ")))

# --- 9b: Remove outlier genes ---
cat("\nFiltering outlier genes...\n")

# BayesPrism recommendation: exclude ribosomal and mitochondrial genes
genes_to_remove <- grep("^RP[SL]|^MT-|^MTRNR", colnames(sc_dat), value = TRUE)
cat(sprintf("  Removing %d ribosomal/mitochondrial genes\n", length(genes_to_remove)))

if (length(genes_to_remove) > 0) {
  keep_genes_sc <- setdiff(colnames(sc_dat), genes_to_remove)
  sc_dat <- sc_dat[, keep_genes_sc]
}

# Remove the same genes from the bulk side
genes_to_remove_bk <- grep("^RP[SL]|^MT-|^MTRNR", colnames(bk_dat), value = TRUE)
if (length(genes_to_remove_bk) > 0) {
  keep_genes_bk <- setdiff(colnames(bk_dat), genes_to_remove_bk)
  bk_dat <- bk_dat[, keep_genes_bk]
}

cat(sprintf("  sc.dat after filter: %d cells × %d genes\n", nrow(sc_dat), ncol(sc_dat)))
cat(sprintf("  bk.dat after filter: %d samples × %d genes\n", nrow(bk_dat), ncol(bk_dat)))

# --- 9c: BayesPrism QC plots ---
cat("\nRunning BayesPrism outlier detection...\n")

sc_stat <- plot.scRNA.outlier(
  input = sc_dat,
  cell.type.labels = cell_type_labels,
  species = "hs",
  return.raw = TRUE
)

bk_stat <- plot.bulk.outlier(
  bulk.input = bk_dat,
  sc.input = sc_dat,
  cell.type.labels = cell_type_labels,
  species = "hs",
  return.raw = TRUE
)

# Restrict genes across both scRNA-seq and bulk
sc_dat_filtered <- cleanup.genes(
  input = sc_dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)

cat(sprintf("  After cleanup: %d genes\n", ncol(sc_dat_filtered)))

# --- 9d: Build Prism object ---
cat("\nConstructing BayesPrism object...\n")

myPrism <- new.prism(
  reference = sc_dat_filtered,
  mixture = bk_dat,
  input.type = "count.matrix",
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_state_labels,
  key = NULL,   # no malignant cell key (non-tumor samples)
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

cat("✓ Prism object created\n")
cat(sprintf("  Reference cell types: %s\n",
            paste(sort(unique(cell_type_labels)), collapse = ", ")))

# --- 9e: Run BayesPrism ---
cat("\nRunning BayesPrism deconvolution (this may take 10-30 minutes)...\n")
cat("  Using Gibbs sampling for posterior estimation...\n")

bp_result <- run.prism(
  prism = myPrism,
  n.cores = parallel::detectCores() - 1,
  update.gibbs = TRUE
)

cat("✓ BayesPrism deconvolution complete!\n")

# --- 9f: Extract results ---
# Cell type fractions (θ)
theta <- get.fraction(
  bp = bp_result,
  which.theta = "final",
  state.or.type = "type"
)

cat("\n=== BayesPrism Cell Type Fractions (%) ===\n")
cat(sprintf("%-20s %10s %10s %10s %10s\n", "Cell Type", "Mean", "SD", "Min", "Max"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (ct in colnames(theta)) {
  ct_vals <- theta[, ct] * 100
  cat(sprintf("%-20s %9.1f%% %9.1f%% %9.1f%% %9.1f%%\n",
              ct, mean(ct_vals), sd(ct_vals), min(ct_vals), max(ct_vals)))
}

# Cell state fractions (fine subtypes)
theta_state <- get.fraction(
  bp = bp_result,
  which.theta = "final",
  state.or.type = "state"
)

cat("\n=== Cell State Fractions (%) ===\n")
cat(sprintf("%-25s %10s\n", "Cell State", "Mean %"))
cat(paste(rep("-", 40), collapse = ""), "\n")
for (cs in colnames(theta_state)) {
  cat(sprintf("%-25s %9.1f%%\n", cs, mean(theta_state[, cs]) * 100))
}

# Save results (RNA-fraction, uncorrected)
write.csv(theta, file.path(output_dir, "celltype_proportions_BayesPrism.csv"))
write.csv(theta_state, file.path(output_dir, "cellstate_proportions_BayesPrism.csv"))

gc()

# =====================================
# PART 9.5: RNA Content correction (RNA fraction -> Cell fraction)
# =====================================
#
# BayesPrism estimates RNA contribution fractions, but FCM measures cell-number
# fractions. Because macrophages contain 5-10x more RNA per cell than lymphocytes,
# converting RNA fraction -> cell fraction is necessary.
#
# Correction factors are estimated directly from the reference scRNA-seq (data-driven):
#   RNA_content_i = median(total UMI per cell) for cell type i
#   cell_fraction_i = (theta_i / RNA_content_i) / sum(theta_j / RNA_content_j)
#
# Reference:
#   Tsoucas et al. (2019) Nature Communications: BayesPrism considers this
#   Newman et al. (2019) Nature Biotechnology: CIBERSORTx S-mode addresses same issue
#

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 9.5: RNA Content Correction                           ║\n")
cat("║  RNA fraction -> Cell fraction (for FCM comparison)         ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Step 1: estimate median RNA content per cell type from the scRNA-seq reference
cat("Estimating RNA content per cell type from scRNA-seq reference...\n")
rna_content <- tapply(reference_obj$nCount_RNA, reference_obj$celltype, median)
cat("\nMedian UMI counts per cell type:\n")
for (ct in sort(names(rna_content))) {
  cat(sprintf("  %-15s: %8.0f UMI\n", ct, rna_content[ct]))
}

# Normalize to a baseline (T cell = 1.0)
if ("T_cell" %in% names(rna_content)) {
  baseline <- rna_content["T_cell"]
} else {
  baseline <- median(rna_content)
}
rna_ratio <- rna_content / as.numeric(baseline)

cat("\nRNA content ratios (relative to T cell):\n")
for (ct in sort(names(rna_ratio))) {
  cat(sprintf("  %-15s: %.2f x\n", ct, rna_ratio[ct]))
}

# Step 2: convert RNA fraction -> Cell fraction
cat("\nConverting RNA fractions to cell fractions...\n")

# Apply correction to each row (sample) of theta
theta_cell <- theta  # copy
for (i in 1:nrow(theta)) {
  # Divide each cell type's RNA fraction by its RNA content ratio
  cell_fraction_raw <- numeric(ncol(theta))
  names(cell_fraction_raw) <- colnames(theta)
  
  for (ct in colnames(theta)) {
    if (ct %in% names(rna_ratio) && rna_ratio[ct] > 0) {
      cell_fraction_raw[ct] <- theta[i, ct] / rna_ratio[ct]
    } else {
      cell_fraction_raw[ct] <- theta[i, ct]
    }
  }
  
  # Normalize so the values sum to 1
  total <- sum(cell_fraction_raw)
  if (total > 0) {
    theta_cell[i, ] <- cell_fraction_raw / total
  }
}

# Compare before vs after correction
cat("\n=== RNA fraction vs Cell fraction (Mean %) ===\n")
cat(sprintf("%-15s %12s %12s %10s\n", "Cell Type", "RNA frac", "Cell frac", "Change"))
cat(paste(rep("-", 55), collapse = ""), "\n")
for (ct in colnames(theta)) {
  rna_pct <- mean(theta[, ct]) * 100
  cell_pct <- mean(theta_cell[, ct]) * 100
  change <- cell_pct - rna_pct
  cat(sprintf("%-15s %11.1f%% %11.1f%% %+9.1f%%\n", ct, rna_pct, cell_pct, change))
}

# Save corrected results
write.csv(theta_cell, file.path(output_dir, "celltype_proportions_cellFraction.csv"))
cat("\n✓ Cell fraction proportions saved\n")

# =====================================
# PART 10: FCM validation (multi-metric)
# =====================================
#
# Validation metrics:
#   1. Spearman rank correlation (robust to nonlinear relationships)
#   2. Pearson correlation (for comparison with conventional methods)
#   3. Group-difference direction concordance (RA-ILD vs Control)
#   4. Bland-Altman analysis (assessment of systematic bias)
#
# Validated on both the RNA-uncorrected (theta) and RNA-corrected (theta_cell) versions

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 10: Flow Cytometry Validation (Multi-metric)          ║\n")
cat("║  Spearman + Pearson + Group concordance + Bland-Altman      ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Load FCM data
fcm_path <- "./FCM_integrated_data_transformed.xlsx"

if (file.exists(fcm_path)) {
  fcm_data <- read_excel(fcm_path, sheet = "BALF_analysis")
  cat(sprintf("FCM data loaded: %d samples\n", nrow(fcm_data)))
  
  # Identify the sample-name column
  id_col <- colnames(fcm_data)[1]
  fcm_data$sample_id <- fcm_data[[id_col]]
  
  # Identify FCM Macrophage/Lymphocyte/Neutrophil columns
  # Prefer _Percent columns (whole-cell fraction %); avoid mistakenly picking subtype marker columns.
  fcm_mac_col <- "BALF_Macrophage_Percent"
  fcm_lym_col <- "BALF_Lymphocyte_Percent"
  fcm_neu_col <- "BALF_Neutrophil_Percent"

  # Fallback when a column is missing
  if (!fcm_mac_col %in% colnames(fcm_data)) {
    fcm_mac_col <- grep("Macrophage_Percent|Macrophage%", colnames(fcm_data), value = TRUE, ignore.case = TRUE)[1]
  }
  if (!fcm_lym_col %in% colnames(fcm_data)) {
    fcm_lym_col <- grep("Lymphocyte_Percent|Lymphocyte%", colnames(fcm_data), value = TRUE, ignore.case = TRUE)[1]
  }
  if (!fcm_neu_col %in% colnames(fcm_data)) {
    fcm_neu_col <- grep("Neutrophil_Percent|Neutrophil%", colnames(fcm_data), value = TRUE, ignore.case = TRUE)[1]
  }

  cat(sprintf("  FCM columns: Mac=%s, Lymph=%s, Neut=%s\n",
              fcm_mac_col, fcm_lym_col, fcm_neu_col))

  # Common samples
  bp_samples <- rownames(theta)
  fcm_samples <- as.character(fcm_data$sample_id)
  common_samples <- intersect(bp_samples, fcm_samples)
  cat(sprintf("  Common samples: %d\n", length(common_samples)))

  if (length(common_samples) > 5) {

    # --- Retrieve FCM values ---
    fcm_idx <- match(common_samples, fcm_data$sample_id)
    fcm_mac_vals  <- as.numeric(fcm_data[[fcm_mac_col]][fcm_idx])
    fcm_lym_vals  <- as.numeric(fcm_data[[fcm_lym_col]][fcm_idx])
    fcm_neut_vals <- as.numeric(fcm_data[[fcm_neu_col]][fcm_idx])
    
    # --- BayesPrism estimates (RNA-uncorrected & corrected) ---
    lymph_cols <- intersect(c("T_cell", "NK", "B_cell", "Plasma"), colnames(theta))
    
    # RNA fraction version
    bp_mac_rna  <- theta[common_samples, "Macrophage"] * 100
    bp_lymph_rna <- rowSums(theta[common_samples, lymph_cols, drop = FALSE]) * 100
    bp_neut_rna <- if ("Neutrophil" %in% colnames(theta)) {
      theta[common_samples, "Neutrophil"] * 100
    } else { rep(0, length(common_samples)) }
    
    # Cell fraction version (RNA content corrected)
    bp_mac_cell  <- theta_cell[common_samples, "Macrophage"] * 100
    bp_lymph_cell <- rowSums(theta_cell[common_samples, lymph_cols, drop = FALSE]) * 100
    bp_neut_cell <- if ("Neutrophil" %in% colnames(theta_cell)) {
      theta_cell[common_samples, "Neutrophil"] * 100
    } else { rep(0, length(common_samples)) }
    
    # --- Correlation analysis (Spearman + Pearson; RNA & Cell versions) ---
    cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
    cat("║  Correlation Analysis: BayesPrism vs FCM                      ║\n")
    cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
    
    run_correlation <- function(bp_vals, fcm_vals, label) {
      valid <- !is.na(fcm_vals) & !is.na(bp_vals)
      n <- sum(valid)
      if (n < 5) return(list(spearman_r = NA, spearman_ci_lo = NA, spearman_ci_hi = NA,
                              pearson_r = NA, n = n))

      sp <- cor.test(bp_vals[valid], fcm_vals[valid], method = "spearman")
      pe <- cor.test(bp_vals[valid], fcm_vals[valid], method = "pearson")

      # Bootstrap 95% CI for Spearman
      set.seed(42)
      boot_rho <- replicate(2000, {
        idx <- sample.int(n, replace = TRUE)
        cor(bp_vals[valid][idx], fcm_vals[valid][idx], method = "spearman")
      })
      sp_ci <- quantile(boot_rho, c(0.025, 0.975), na.rm = TRUE)

      list(
        spearman_r = sp$estimate, spearman_p = sp$p.value,
        spearman_ci_lo = sp_ci[1], spearman_ci_hi = sp_ci[2],
        pearson_r = pe$estimate, pearson_p = pe$p.value,
        n = n, bp = bp_vals[valid], fcm = fcm_vals[valid]
      )
    }
    
    # RNA fraction version
    res_mac_rna <- run_correlation(bp_mac_rna, fcm_mac_vals, "Macrophage_RNA")
    res_lym_rna <- run_correlation(bp_lymph_rna, fcm_lym_vals, "Lymphocyte_RNA")
    res_neu_rna <- run_correlation(bp_neut_rna, fcm_neut_vals, "Neutrophil_RNA")
    
    # Cell fraction version
    res_mac_cell <- run_correlation(bp_mac_cell, fcm_mac_vals, "Macrophage_Cell")
    res_lym_cell <- run_correlation(bp_lymph_cell, fcm_lym_vals, "Lymphocyte_Cell")
    res_neu_cell <- run_correlation(bp_neut_cell, fcm_neut_vals, "Neutrophil_Cell")
    
    cat("--- RNA Fraction (uncorrected) ---\n")
    cat(sprintf("%-15s %12s %12s %12s %12s %5s\n", "Cell Type", "Spearman ρ", "p-value", "Pearson r", "p-value", "n"))
    cat(paste(rep("-", 72), collapse = ""), "\n")
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Macrophage",
                res_mac_rna$spearman_r, res_mac_rna$spearman_p,
                res_mac_rna$pearson_r, res_mac_rna$pearson_p, res_mac_rna$n))
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Lymphocyte",
                res_lym_rna$spearman_r, res_lym_rna$spearman_p,
                res_lym_rna$pearson_r, res_lym_rna$pearson_p, res_lym_rna$n))
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Neutrophil",
                res_neu_rna$spearman_r, res_neu_rna$spearman_p,
                res_neu_rna$pearson_r, res_neu_rna$pearson_p, res_neu_rna$n))
    
    cat("\n--- Cell Fraction (RNA content corrected) ---\n")
    cat(sprintf("%-15s %12s %12s %12s %12s %5s\n", "Cell Type", "Spearman ρ", "p-value", "Pearson r", "p-value", "n"))
    cat(paste(rep("-", 72), collapse = ""), "\n")
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Macrophage",
                res_mac_cell$spearman_r, res_mac_cell$spearman_p,
                res_mac_cell$pearson_r, res_mac_cell$pearson_p, res_mac_cell$n))
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Lymphocyte",
                res_lym_cell$spearman_r, res_lym_cell$spearman_p,
                res_lym_cell$pearson_r, res_lym_cell$pearson_p, res_lym_cell$n))
    cat(sprintf("%-15s %12.3f %12.4f %12.3f %12.4f %5d\n", "Neutrophil",
                res_neu_cell$spearman_r, res_neu_cell$spearman_p,
                res_neu_cell$pearson_r, res_neu_cell$pearson_p, res_neu_cell$n))
    
    # --- Group-difference direction concordance analysis ---
    cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
    cat("║  Group Concordance: RA-ILD vs Control direction agreement     ║\n")
    cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
    
    ky_mask <- grepl("^(KY[0-9]|RA[0-9])", common_samples)
    kyc_mask <- grepl("^(KYC|Sarcoidosis)", common_samples)
    
    if (sum(ky_mask) >= 3 && sum(kyc_mask) >= 3) {
      cat(sprintf("  RA-ILD samples: %d, Control samples: %d\n\n", sum(ky_mask), sum(kyc_mask)))
      
      concordance_table <- data.frame(
        CellType = character(),
        FCM_KY_mean = numeric(), FCM_KYC_mean = numeric(), FCM_direction = character(),
        BP_KY_mean = numeric(), BP_KYC_mean = numeric(), BP_direction = character(),
        Concordant = character(),
        stringsAsFactors = FALSE
      )
      
      # Group comparison using the cell fraction version
      compare_groups <- function(bp_vals, fcm_vals, ky_mask, kyc_mask, cell_name) {
        fcm_ky <- mean(fcm_vals[ky_mask], na.rm = TRUE)
        fcm_kyc <- mean(fcm_vals[kyc_mask], na.rm = TRUE)
        fcm_dir <- ifelse(fcm_ky > fcm_kyc, "KY > KYC", "KY < KYC")
        
        bp_ky <- mean(bp_vals[ky_mask], na.rm = TRUE)
        bp_kyc <- mean(bp_vals[kyc_mask], na.rm = TRUE)
        bp_dir <- ifelse(bp_ky > bp_kyc, "KY > KYC", "KY < KYC")
        
        concordant <- ifelse(fcm_dir == bp_dir, "✓ YES", "✗ NO")
        
        data.frame(
          CellType = cell_name,
          FCM_KY_mean = round(fcm_ky, 1),
          FCM_KYC_mean = round(fcm_kyc, 1),
          FCM_direction = fcm_dir,
          BP_KY_mean = round(bp_ky, 1),
          BP_KYC_mean = round(bp_kyc, 1),
          BP_direction = bp_dir,
          Concordant = concordant,
          stringsAsFactors = FALSE
        )
      }
      
      concordance_table <- rbind(
        compare_groups(bp_mac_cell, fcm_mac_vals, ky_mask, kyc_mask, "Macrophage"),
        compare_groups(bp_lymph_cell, fcm_lym_vals, ky_mask, kyc_mask, "Lymphocyte"),
        compare_groups(bp_neut_cell, fcm_neut_vals, ky_mask, kyc_mask, "Neutrophil")
      )
      
      cat(sprintf("%-12s %10s %10s %12s %10s %10s %12s %10s\n",
                  "CellType", "FCM_KY", "FCM_KYC", "FCM_dir", "BP_KY", "BP_KYC", "BP_dir", "Concordant"))
      cat(paste(rep("-", 95), collapse = ""), "\n")
      for (r in 1:nrow(concordance_table)) {
        cat(sprintf("%-12s %9.1f%% %9.1f%% %12s %9.1f%% %9.1f%% %12s %10s\n",
                    concordance_table$CellType[r],
                    concordance_table$FCM_KY_mean[r], concordance_table$FCM_KYC_mean[r],
                    concordance_table$FCM_direction[r],
                    concordance_table$BP_KY_mean[r], concordance_table$BP_KYC_mean[r],
                    concordance_table$BP_direction[r],
                    concordance_table$Concordant[r]))
      }
      
      n_concordant <- sum(concordance_table$Concordant == "✓ YES")
      cat(sprintf("\nGroup direction concordance: %d/%d (%.0f%%)\n",
                  n_concordant, nrow(concordance_table),
                  n_concordant / nrow(concordance_table) * 100))
      
      write.csv(concordance_table, file.path(output_dir, "FCM_group_concordance.csv"), row.names = FALSE)
    }
    
    # --- Validation plots (6 panels: RNA version x3 + Cell version x3) ---
    cat("\nGenerating validation plots...\n")
    
    tryCatch({
      pdf(file.path(output_dir, "FCM_validation_BayesPrism.pdf"), width = 15, height = 10)
      par(mfrow = c(2, 3), mar = c(5, 5, 4, 2))
      
      plot_validation <- function(bp, fcm, label, color, res) {
        valid <- !is.na(fcm) & !is.na(bp)
        plot(fcm[valid], bp[valid],
             xlab = paste("FCM", label, "(%)"),
             ylab = paste("BayesPrism", label, "(%)"),
             main = sprintf("%s\nSpearman ρ=%.3f (p=%.4f)\nPearson r=%.3f",
                            label, res$spearman_r, res$spearman_p, res$pearson_r),
             pch = 19, col = adjustcolor(color, 0.7), cex = 1.2,
             xlim = c(0, 100), ylim = c(0, 100), cex.main = 0.9)
        abline(0, 1, lty = 2, col = "gray50", lwd = 1.5)
        if (sum(valid) >= 3) {
          abline(lm(bp[valid] ~ fcm[valid]), col = color, lwd = 2)
        }
        # Distinguish KY/KYC by symbol
        ky_pts <- grepl("^(KY[0-9]|RA[0-9])", common_samples) & valid
        kyc_pts <- grepl("^(KYC|Sarcoidosis)", common_samples) & valid
        points(fcm[ky_pts], bp[ky_pts], pch = 19, col = adjustcolor("red", 0.6), cex = 1.2)
        points(fcm[kyc_pts], bp[kyc_pts], pch = 17, col = adjustcolor("blue", 0.6), cex = 1.2)
        legend("topleft", c("RA-ILD (KY)", "Control (KYC)"),
               pch = c(19, 17), col = c("red", "blue"), cex = 0.8, bty = "n")
      }
      
      # Top row: RNA fraction
      plot_validation(bp_mac_rna, fcm_mac_vals, "Macrophage (RNA frac)", "steelblue", res_mac_rna)
      plot_validation(bp_lymph_rna, fcm_lym_vals, "Lymphocyte (RNA frac)", "forestgreen", res_lym_rna)
      plot_validation(bp_neut_rna, fcm_neut_vals, "Neutrophil (RNA frac)", "coral", res_neu_rna)
      
      # Bottom row: Cell fraction (RNA content corrected)
      plot_validation(bp_mac_cell, fcm_mac_vals, "Macrophage (Cell frac)", "steelblue", res_mac_cell)
      plot_validation(bp_lymph_cell, fcm_lym_vals, "Lymphocyte (Cell frac)", "forestgreen", res_lym_cell)
      plot_validation(bp_neut_cell, fcm_neut_vals, "Neutrophil (Cell frac)", "coral", res_neu_cell)
      
      dev.off()
      cat("✓ FCM validation plots saved (2×3 panel)\n")
    }, error = function(e) {
      tryCatch(dev.off(), error = function(x) NULL)
      cat(sprintf("⚠ Plot generation failed: %s\n", e$message))
    })
    
    # --- Validation results CSV (comprehensive) ---
    validation_df <- data.frame(
      Sample = common_samples,
      Group = ifelse(grepl("^(KY[0-9]|RA[0-9])", common_samples), "RA-ILD", "Control"),
      BP_Mac_RNA = bp_mac_rna,
      BP_Lymph_RNA = bp_lymph_rna,
      BP_Neut_RNA = bp_neut_rna,
      BP_Mac_Cell = bp_mac_cell,
      BP_Lymph_Cell = bp_lymph_cell,
      BP_Neut_Cell = bp_neut_cell,
      FCM_Macrophage = fcm_mac_vals,
      FCM_Lymphocyte = fcm_lym_vals,
      FCM_Neutrophil = fcm_neut_vals
    )
    write.csv(validation_df, file.path(output_dir, "FCM_validation_data.csv"), row.names = FALSE)
    
    # Summary CSV (with bootstrap 95% CI)
    summary_corr <- data.frame(
      CellType = rep(c("Macrophage", "Lymphocyte", "Neutrophil"), 2),
      Version = rep(c("RNA_fraction", "Cell_fraction"), each = 3),
      FCM_Column = rep(c(fcm_mac_col, fcm_lym_col, fcm_neu_col), 2),
      Spearman_rho = c(res_mac_rna$spearman_r, res_lym_rna$spearman_r, res_neu_rna$spearman_r,
                       res_mac_cell$spearman_r, res_lym_cell$spearman_r, res_neu_cell$spearman_r),
      Spearman_CI_lo = c(res_mac_rna$spearman_ci_lo, res_lym_rna$spearman_ci_lo, res_neu_rna$spearman_ci_lo,
                          res_mac_cell$spearman_ci_lo, res_lym_cell$spearman_ci_lo, res_neu_cell$spearman_ci_lo),
      Spearman_CI_hi = c(res_mac_rna$spearman_ci_hi, res_lym_rna$spearman_ci_hi, res_neu_rna$spearman_ci_hi,
                          res_mac_cell$spearman_ci_hi, res_lym_cell$spearman_ci_hi, res_neu_cell$spearman_ci_hi),
      Spearman_p = c(res_mac_rna$spearman_p, res_lym_rna$spearman_p, res_neu_rna$spearman_p,
                     res_mac_cell$spearman_p, res_lym_cell$spearman_p, res_neu_cell$spearman_p),
      Pearson_r = c(res_mac_rna$pearson_r, res_lym_rna$pearson_r, res_neu_rna$pearson_r,
                    res_mac_cell$pearson_r, res_lym_cell$pearson_r, res_neu_cell$pearson_r),
      Pearson_p = c(res_mac_rna$pearson_p, res_lym_rna$pearson_p, res_neu_rna$pearson_p,
                    res_mac_cell$pearson_p, res_lym_cell$pearson_p, res_neu_cell$pearson_p),
      N = c(res_mac_rna$n, res_lym_rna$n, res_neu_rna$n,
            res_mac_cell$n, res_lym_cell$n, res_neu_cell$n)
    )
    write.csv(summary_corr, file.path(output_dir, "FCM_correlation_summary.csv"), row.names = FALSE)
    cat("✓ Correlation summary saved\n")
    
  } else {
    cat("  Insufficient common samples for validation\n")
  }
} else {
  cat("  FCM data file not found. Skipping validation.\n")
  cat(sprintf("  Expected path: %s\n", fcm_path))
}

# =====================================
# PART 11: KY vs KYC group comparison
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 11: RA-ILD (KY) vs Control (KYC) Comparison           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

ky_samples  <- grep("^(KY[0-9]|RA[0-9])", rownames(theta), value = TRUE)
kyc_samples <- grep("^(KYC|Sarcoidosis)", rownames(theta), value = TRUE)

cat(sprintf("KY samples (RA-ILD): n=%d\n", length(ky_samples)))
cat(sprintf("KYC samples (Control): n=%d\n", length(kyc_samples)))

# --- 11a: Cell fraction version (RNA content corrected, comparable to FCM) ---
cat("\n--- Cell Fraction (RNA content corrected) ---\n")
cat(sprintf("%-20s %12s %12s %10s %8s\n",
            "Cell Type", "KY (RA-ILD)", "KYC (Ctrl)", "p-value", "Sig"))
cat(paste(rep("-", 65), collapse = ""), "\n")

comparison_results <- data.frame()

for (ct in colnames(theta_cell)) {
  ky_vals  <- theta_cell[ky_samples, ct] * 100
  kyc_vals <- theta_cell[kyc_samples, ct] * 100
  
  if (sd(ky_vals) > 0 & sd(kyc_vals) > 0 & length(ky_vals) >= 3 & length(kyc_vals) >= 3) {
    wtest <- wilcox.test(ky_vals, kyc_vals, exact=TRUE)
    pval <- wtest$p.value
    sig <- ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01, "**",
                         ifelse(pval < 0.05, "*", "ns")))
  } else {
    pval <- NA
    sig <- "NA"
  }
  
  cat(sprintf("%-20s %11.1f%% %11.1f%% %10s %8s\n",
              ct, mean(ky_vals), mean(kyc_vals),
              ifelse(is.na(pval), "NA", sprintf("%.4f", pval)), sig))
  
  comparison_results <- rbind(comparison_results, data.frame(
    CellType = ct,
    Version = "Cell_fraction",
    KY_Mean = mean(ky_vals),
    KYC_Mean = mean(kyc_vals),
    KY_SD = sd(ky_vals),
    KYC_SD = sd(kyc_vals),
    p_value = pval,
    Significance = sig
  ))
}

# --- 11b: RNA fraction version (reference) ---
cat("\n--- RNA Fraction (uncorrected, for reference) ---\n")
cat(sprintf("%-20s %12s %12s %10s %8s\n",
            "Cell Type", "KY (RA-ILD)", "KYC (Ctrl)", "p-value", "Sig"))
cat(paste(rep("-", 65), collapse = ""), "\n")

for (ct in colnames(theta)) {
  ky_vals  <- theta[ky_samples, ct] * 100
  kyc_vals <- theta[kyc_samples, ct] * 100
  
  if (sd(ky_vals) > 0 & sd(kyc_vals) > 0 & length(ky_vals) >= 3 & length(kyc_vals) >= 3) {
    wtest <- wilcox.test(ky_vals, kyc_vals, exact=TRUE)
    pval <- wtest$p.value
    sig <- ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01, "**",
                         ifelse(pval < 0.05, "*", "ns")))
  } else {
    pval <- NA
    sig <- "NA"
  }
  
  cat(sprintf("%-20s %11.1f%% %11.1f%% %10s %8s\n",
              ct, mean(ky_vals), mean(kyc_vals),
              ifelse(is.na(pval), "NA", sprintf("%.4f", pval)), sig))
  
  comparison_results <- rbind(comparison_results, data.frame(
    CellType = ct,
    Version = "RNA_fraction",
    KY_Mean = mean(ky_vals),
    KYC_Mean = mean(kyc_vals),
    KY_SD = sd(ky_vals),
    KYC_SD = sd(kyc_vals),
    p_value = pval,
    Significance = sig
  ))
}

write.csv(comparison_results, file.path(output_dir, "KY_vs_KYC_comparison_BayesPrism.csv"),
          row.names = FALSE)

# Group comparison plot (cell fraction version shown as the main result)
theta_long <- as.data.frame(theta_cell) %>%
  rownames_to_column("Sample") %>%
  mutate(Group = ifelse(grepl("^(KYC|Sarcoidosis)", Sample), "Control (KYC)", "RA-ILD (KY)")) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "CellType", values_to = "Fraction")

p_comparison <- ggplot(theta_long, aes(x = CellType, y = Fraction * 100, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5, size = 1.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(title = "Cell Type Proportions: RA-ILD vs Control (BayesPrism, RNA content corrected)",
       x = "", y = "Proportion (%)") +
  scale_fill_manual(values = c("Control (KYC)" = "steelblue", "RA-ILD (KY)" = "coral"))

ggsave(file.path(output_dir, "KY_vs_KYC_boxplot_BayesPrism.pdf"), p_comparison, width = 14, height = 8)

# =====================================
# PART 12: Visualization and saving results
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║  PART 12: Visualization and Results                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

# Heatmap (pheatmap generates a temporary PNG internally, so set the device explicitly)
tryCatch({
  pdf(file.path(output_dir, "celltype_heatmap_BayesPrism.pdf"), width = 16, height = 10)
  pheatmap(
    t(theta),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
    main = "RA-ILD BAL Cell Type Proportions (BayesPrism)",
    fontsize_row = 10,
    fontsize_col = 8,
    angle_col = 45
  )
  dev.off()
  cat("✓ Heatmap saved\n")
}, error = function(e) {
  tryCatch(dev.off(), error = function(x) NULL)
  cat(sprintf("⚠ Heatmap generation failed: %s\n", e$message))
  cat("  Attempting with cairo_pdf...\n")
  tryCatch({
    cairo_pdf(file.path(output_dir, "celltype_heatmap_BayesPrism.pdf"), width = 16, height = 10)
    pheatmap(
      t(theta),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
      main = "RA-ILD BAL Cell Type Proportions (BayesPrism)",
      fontsize_row = 10,
      fontsize_col = 8,
      angle_col = 45
    )
    dev.off()
    cat("✓ Heatmap saved (cairo_pdf)\n")
  }, error = function(e2) {
    tryCatch(dev.off(), error = function(x) NULL)
    cat(sprintf("⚠ Heatmap skipped: %s\n", e2$message))
  })
})

# Stacked bar plot (cell fraction version)
prop_long <- as.data.frame(theta_cell) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")

p_stacked <- ggplot(prop_long, aes(x = Sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Cell Type Proportions (BayesPrism, RNA content corrected)",
       x = "Sample", y = "Proportion", fill = "Cell Type") +
  scale_fill_brewer(palette = "Set3")

ggsave(file.path(output_dir, "celltype_stacked_barplot_BayesPrism.pdf"), p_stacked, width = 16, height = 8)

# Summary statistics (cell fraction version)
summary_df <- data.frame(
  CellType = colnames(theta_cell),
  Mean = colMeans(theta_cell) * 100,
  SD = apply(theta_cell, 2, sd) * 100,
  Min = apply(theta_cell, 2, min) * 100,
  Max = apply(theta_cell, 2, max) * 100,
  Median = apply(theta_cell, 2, median) * 100
) %>% arrange(desc(Mean))

write.csv(summary_df, file.path(output_dir, "celltype_summary_BayesPrism.csv"), row.names = FALSE)

# Record the annotation strategy
annotation_summary <- data.frame(
  Dataset = c("GSE145926", "GSE193782", "GSE184735"),
  stringsAsFactors = FALSE
)
# Tally annotation source per dataset
for (i in 1:3) {
  ds <- annotation_summary$Dataset[i]
  ds_cells <- which(reference_obj$dataset == ds)
  sources <- table(reference_obj$annotation_source[ds_cells])
  annotation_summary$Annotation_Source[i] <- paste(names(sources), collapse = " + ")
  annotation_summary$N_Cells[i] <- length(ds_cells)
  annotation_summary$Source_Detail[i] <- paste(names(sources), sources, sep = "=", collapse = ", ")
}
write.csv(annotation_summary, file.path(output_dir, "annotation_strategy_summary.csv"),
          row.names = FALSE)

# Save all results
saveRDS(list(
  bp_result = bp_result,
  theta_type = theta,
  theta_state = theta_state,
  summary = summary_df,
  comparison = comparison_results,
  reference_info = list(
    total_cells = ncol(reference_obj),
    datasets = table(reference_obj$dataset),
    celltypes = table(reference_obj$celltype),
    celltypes_fine = table(reference_obj$celltype_fine),
    annotation_method = "Author annotations from original publications + canonical BAL marker fallback",
    annotation_sources = table(reference_obj$annotation_source)
  )
), file.path(output_dir, "full_analysis_results_v3.rds"))

# Session info
sink(file.path(output_dir, "session_info_v3.txt"))
cat("Analysis completed:", as.character(Sys.time()), "\n")
cat("Pipeline: v3 (Author Annotation + BayesPrism)\n\n")
cat("=== Key Parameters ===\n")
cat(sprintf("Reference cells: %d\n", ncol(reference_obj)))
cat(sprintf("Bulk samples: %d\n", ncol(bulk_matrix)))
cat(sprintf("Cell types: %s\n", paste(colnames(theta), collapse = ", ")))
cat("Annotation strategy:\n")
for (ds in c("GSE145926", "GSE193782", "GSE184735")) {
  ds_cells <- which(reference_obj$dataset == ds)
  sources <- unique(reference_obj$annotation_source[ds_cells])
  cat(sprintf("  %s: %s (%d cells)\n", ds, paste(sources, collapse = " + "), length(ds_cells)))
}
cat(sprintf("Deconvolution: BayesPrism (Gibbs sampling)\n\n"))
sessionInfo()
sink()

# Save workspace
save.image(file.path(output_dir, "BAL_analysis_v3.RData"))

# =====================================
# PART 13: Final report
# =====================================

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║                    FINAL REPORT v3                           ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")

cat("\n[ANNOTATION STRATEGY]\n")
for (ds in c("GSE145926", "GSE193782", "GSE184735")) {
  ds_cells <- which(reference_obj$dataset == ds)
  sources <- table(reference_obj$annotation_source[ds_cells])
  cat(sprintf("  %s: %s [%d cells]\n", ds,
              paste(names(sources), sources, sep = "=", collapse = ", "),
              length(ds_cells)))
}

cat("\n[REFERENCE DATA]\n")
cat(sprintf("  Total cells: %d\n", ncol(reference_obj)))
cat(sprintf("  Datasets: GSE193782 (Healthy) + GSE184735 (Sarcoidosis) + GSE145926 (COVID-19)\n"))
cat(sprintf("  All datasets: BAL fluid-derived scRNA-seq\n"))

cat("\n[DECONVOLUTION RESULTS (BayesPrism)]\n")
cat(sprintf("%-20s %10s\n", "Cell Type", "Mean %"))
cat(paste(rep("-", 35), collapse = ""), "\n")
for (i in 1:nrow(summary_df)) {
  cat(sprintf("%-20s %9.1f%%\n", summary_df$CellType[i], summary_df$Mean[i]))
}

# Show FCM validation results if available
if (exists("validation_results") && length(validation_results) > 0) {
  cat("\n[FCM VALIDATION]\n")
  for (ct_name in names(validation_results)) {
    vr <- validation_results[[ct_name]]
    cat(sprintf("  %s: r=%.3f (p=%.4f)\n", ct_name, vr$estimate, vr$p.value))
  }
  
  # Target check
  all_r <- sapply(validation_results, function(x) x$estimate)
  if (all(all_r > 0.6)) {
    cat("\n  ✓ r>0.6 for all cell types -> estimates for FCM-unmeasured cells are also reliable\n")
  } else {
    low_r <- names(all_r)[all_r <= 0.6]
    cat(sprintf("\n  ⚠ Cell types with r<=0.6: %s -> consider further reference improvement\n",
                paste(low_r, collapse = ", ")))
  }
}

cat("\n[METHODS SECTION (for manuscript)]\n")
cat("  Cell type deconvolution was performed using BayesPrism (Chu et al.,\n")
cat("  Nature Cancer, 2022). Three publicly available BAL fluid scRNA-seq\n")
cat("  datasets were used as reference: GSE193782 (healthy controls),\n")
cat("  GSE184735 (sarcoidosis), and GSE145926 (COVID-19, Liao et al.,\n")
cat("  Nature Medicine, 2020). Cell type annotations were adopted from\n")
cat("  the original publications: for GSE145926, annotations were obtained\n")
cat("  from the authors' GitHub repository (github.com/zhangzlab/covid_balf),\n")
cat("  and cells not present in the published annotation were excluded;\n")
cat("  for GSE193782, cell type labels from the published metadata were used.\n")
cat("  For GSE184735, which lacked author-provided annotations, cell type\n")
cat("  labels were transferred from the annotated datasets using k-nearest\n")
cat("  neighbor classification (k=20) in the Harmony-corrected embedding\n")
cat("  space. To convert RNA-fraction estimates to cell-fraction proportions\n")
cat("  comparable with flow cytometry, we applied post-hoc RNA content\n")
cat("  correction using median UMI counts per cell type from the scRNA-seq\n")
cat("  reference (Jew et al., Nature Communications, 2020). Deconvolution\n")
cat("  estimates were validated against paired flow cytometry measurements\n")
cat("  using Spearman rank correlation and concordance of group-level\n")
cat("  differences (RA-ILD vs control).\n")

cat("\n[OUTPUT FILES]\n")
cat(sprintf("Directory: %s\n\n", output_dir))
output_files <- list.files(output_dir)
for (f in output_files) {
  cat(sprintf("  ✓ %s\n", f))
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║              ANALYSIS COMPLETE (v3)                          ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")
