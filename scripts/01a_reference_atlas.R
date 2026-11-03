# 01a_reference_atlas.R -- scRNA-seq BAL reference atlas from public GEO data (GSE145926, GSE193782, GSE184735): merge, QC, Harmony
#   integration and cell-type annotation (author annotations, kNN label transfer, canonical-marker fallback). Figure 1B-E are drawn from it
#   by figures/panels_fig1_atlas.R, and 01b_BayesPrism_n31.R deconvolves the bulk RNA-seq with it.
# Produces: output_BayesPrism/BAL_reference_author_annotated.rds. About 32 GB RAM.
# Initial setup

# Temporary directory for large intermediate files
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
library(harmony)
library(patchwork)
library(hdf5r)
library(future)

plan("sequential")

# Working directory
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
source("scripts/00_cohort_config.R")   # sample_group(), is_sarc(), public_label(): groups from cohort_definition.csv (study IDs)

# Output directory
output_dir <- "./output_BayesPrism"
dir.create(output_dir, showWarnings = FALSE)

# Data directory (reference annotation files, etc.)
data_dir <- "."

# PART 1: Load GSE145926 (COVID-19 BALF) data

covid_processed_path <- "./output_integrated/GSE145926_COVID_processed.rds"

if (file.exists(covid_processed_path)) {
  
  covid_seurat <- readRDS(covid_processed_path)
  
  if (!"sample_id" %in% colnames(covid_seurat@meta.data)) {
    covid_seurat$sample_id <- covid_seurat$sample
  }
  
  # Preserve any pre-existing celltype (for later comparison with author annotation)
  if ("celltype" %in% colnames(covid_seurat@meta.data)) {
    covid_seurat$celltype_alt <- covid_seurat$celltype
    covid_seurat$celltype <- NULL
  }
  
} else {
  
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
  
  seurat_list_covid <- list()
  
  for (h5_file in h5_files) {
    filename <- basename(h5_file)
    gsm <- regmatches(filename, regexpr("GSM[0-9]+", filename))
    idx <- which(sample_info_covid$GSM == gsm)
    if (length(idx) == 0) next
    
    sample_name <- sample_info_covid$SampleName[idx]
    disease <- sample_info_covid$Disease[idx]
    
    counts <- Read10X_h5(h5_file)
    seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name,
                                     min.cells = 3, min.features = 200)
    seurat_obj$sample <- sample_name
    seurat_obj$sample_id <- sample_name
    seurat_obj$disease_status <- disease
    seurat_obj$condition <- "COVID"
    seurat_obj$dataset <- "GSE145926"
    seurat_list_covid[[sample_name]] <- seurat_obj
    
  }
  
  covid_seurat <- merge(seurat_list_covid[[1]], y = seurat_list_covid[-1],
                        add.cell.ids = names(seurat_list_covid))
  
  rm(seurat_list_covid); gc()
  
  # QC
  covid_seurat[["percent.mt"]] <- PercentageFeatureSet(covid_seurat, pattern = "^MT-")
  covid_seurat <- subset(covid_seurat,
                         subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
  
  # Preprocessing
  covid_seurat <- NormalizeData(covid_seurat, verbose = FALSE)
  covid_seurat <- FindVariableFeatures(covid_seurat, nfeatures = 3000, verbose = FALSE)
  covid_seurat <- ScaleData(covid_seurat, verbose = FALSE)
  covid_seurat <- RunPCA(covid_seurat, npcs = 50, verbose = FALSE)
  covid_seurat <- RunUMAP(covid_seurat, dims = 1:30, verbose = FALSE)
  covid_seurat <- FindNeighbors(covid_seurat, dims = 1:30, verbose = FALSE)
  covid_seurat <- FindClusters(covid_seurat, resolution = 0.5, verbose = FALSE)
  
  saveRDS(covid_seurat, covid_processed_path)
  
}

gc()

# PART 2: GSE193782 (Healthy BAL)

# --- Memory-efficient loading ---
# GSE193782 contains large assays (ADT/ITG/SCT). Use a pre-trimmed lite RDS if
# available; otherwise load the original file and trim it immediately.

healthy_lite_path <- "./output_BayesPrism/GSE193782_lite.rds"

if (file.exists(healthy_lite_path)) {
  
  healthy_seurat <- readRDS(healthy_lite_path)
  
} else {
  
  # First try a direct read
  healthy_seurat <- tryCatch({
    obj <- readRDS("./analysis/GSE193782_SeuratObject.rds")
    
    # Immediately drop unneeded assays to free memory
    DefaultAssay(obj) <- "RNA"
    for (assay_name in c("ADT", "ITG", "SCT")) {
      if (assay_name %in% names(obj@assays)) {
        obj[[assay_name]] <- NULL
        
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
    
  }
  gc()
  
  # Save lite version (enables fast loading next time)
  
  dir.create(dirname(healthy_lite_path), showWarnings = FALSE)
  saveRDS(healthy_seurat, healthy_lite_path)
  
}

healthy_seurat$condition <- "Healthy"
healthy_seurat$sample_id <- ifelse("orig.ident" %in% colnames(healthy_seurat@meta.data),
                                   as.character(healthy_seurat$orig.ident), "Healthy")
healthy_seurat$dataset <- "GSE193782"
gc()

# PART 3: GSE184735 (Sarcoidosis BAL)

load_sarcoid <- function(data_path, max_per_sample = 1000) {
  all_files <- list.files(data_path, recursive = TRUE, full.names = TRUE)
  mtx_files <- all_files[grepl("matrix\\.mtx", all_files)]
  mtx_dirs <- unique(dirname(mtx_files))
  seurat_list <- list()
  
  for (i in seq_along(mtx_dirs)) {
    path_parts <- strsplit(mtx_dirs[i], "/")[[1]]
    sample_name <- path_parts[grep("^GSM", path_parts)][1]
    if (is.na(sample_name)) sample_name <- basename(dirname(mtx_dirs[i]))
    
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

# PART 4: Merge 3 datasets

common_genes <- Reduce(intersect, list(
  rownames(healthy_seurat), rownames(sarcoid_merged), rownames(covid_seurat)
))

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

print(table(reference_obj$dataset))

# PART 5: QC filtering

reference_obj[["percent.mt"]] <- PercentageFeatureSet(reference_obj, pattern = "^MT-")
reference_obj[["percent.ribo"]] <- PercentageFeatureSet(reference_obj, pattern = "^RP[SL]")

reference_obj <- subset(reference_obj,
                        subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                          nCount_RNA > 500 & percent.mt < 20)

gc()

# PART 6: Normalization, dimensionality reduction, Harmony integration

reference_obj <- NormalizeData(reference_obj, verbose = FALSE)

reference_obj <- FindVariableFeatures(reference_obj, nfeatures = 2000, verbose = FALSE)

reference_obj <- ScaleData(reference_obj, features = VariableFeatures(reference_obj), verbose = FALSE)

reference_obj <- RunPCA(reference_obj, npcs = 50, verbose = FALSE)

# Harmony batch correction (for visualization; not used for deconvolution)
reference_obj <- RunHarmony(reference_obj, "dataset", verbose = TRUE)

reference_obj <- RunUMAP(reference_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
reference_obj <- FindNeighbors(reference_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
reference_obj <- FindClusters(reference_obj, resolution = 0.8, verbose = FALSE)

# Join layers (Seurat v5)
reference_obj <- JoinLayers(reference_obj)

gc()

# PART 7: Author annotation
#
# Strategy: use the peer-reviewed cell type annotations from each public
#   dataset's original publication directly
#   ("Cell type labels from the original publications were used"),
#   ensuring consistency with the source publications' annotations.
#
# Fallback (only when author annotation is absent from metadata):
#   scoring-based classification using BAL canonical marker genes.
#

# ===== Remove leftover columns from any prior pipeline =====
# Columns carried over from a previous run can be falsely matched during
# annotation lookup, so remove them completely here.
extra_remove_cols <- c("celltype", "celltype_alt", "celltype_fine",
                    "Neutrophil_score1", "Alveolar_Mac_score1",
                    "Monocyte_Mac_score1", "T_cell_score1",
                    "NK_score1", "B_cell_score1", "Epithelial_score1")

for (col_name in extra_remove_cols) {
  if (col_name %in% colnames(reference_obj@meta.data)) {
    reference_obj@meta.data[[col_name]] <- NULL
    
  }
}

# Initialize fresh annotation columns (start from NA)
reference_obj$celltype <- NA_character_
reference_obj$annotation_source <- NA_character_
reference_obj$original_author_label <- NA_character_

# --- 7a: Define canonical BAL marker genes & harmonization function ---
# Define functions/constants needed before applying author annotations

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

# GSE145926: author annotation from Liao et al. Nature Medicine 2020
# Source: https://github.com/zhangzlab/covid_balf/
#   - all.cell.annotation.meta.txt (all cell types)
#   - myeloid.cell.annotation.meta.txt (macrophage subtypes)

liao_annotation_path <- file.path(data_dir, "reference", "GSE145926_all_cell_annotation.txt")
liao_myeloid_path <- file.path(data_dir, "reference", "GSE145926_myeloid_annotation.txt")

for (f in c(liao_annotation_path, liao_myeloid_path)) if (!file.exists(f))
  stop("Liao et al. (2020) annotation file not found: ", f,
       "\n  Download all.cell.annotation.meta.txt and myeloid.cell.annotation.meta.txt from https://github.com/zhangzlab/covid_balf",
       " and save them under reference/ with these names (see data_examples/README.md). The marker-based fallback is not used for GSE145926.")
if (file.exists(liao_annotation_path)) {
  
  liao_meta <- read.delim(liao_annotation_path, sep = "\t", header = TRUE)
  
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

    # Handle unmatched cells
    # Likely QC-filtered by Liao et al., so exclude from the reference
    # (applying the original authors' quality standards, not cherry-picking)
    unmatched_cells <- covid_cells[!matched_mask]
    if (length(unmatched_cells) > 0) {
      
      cells_to_keep <- setdiff(colnames(reference_obj), colnames(reference_obj)[unmatched_cells])
      reference_obj <- subset(reference_obj, cells = cells_to_keep)
      gc()
      
    }
  }
  
  # Also load myeloid subtypes (for macrophage subclassification)
  if (file.exists(liao_myeloid_path)) {
    liao_myeloid <- read.delim(liao_myeloid_path, sep = "\t", header = TRUE)
    
    print(table(liao_myeloid$celltype))
  }
  
}

# GSE193782: "Cell.Types" column in RDS metadata (author annotation)
# GSE184735: no author annotation -> kNN transfer

# GSE193782
gse193782_cells <- which(reference_obj$dataset == "GSE193782")
if ("Cell.Types" %in% colnames(reference_obj@meta.data)) {
  ct_labels <- reference_obj@meta.data[gse193782_cells, "Cell.Types"]
  n_valid <- sum(!is.na(ct_labels))
  if (n_valid > 0) {
    
    reference_obj$original_author_label[gse193782_cells] <- as.character(ct_labels)
    harmonized_193782 <- harmonize_celltype(ct_labels)
    reference_obj$celltype[gse193782_cells] <- harmonized_193782
    reference_obj$annotation_source[gse193782_cells] <- "author_GSE193782"
    
  }
}

# GSE184735
gse184735_cells <- which(reference_obj$dataset == "GSE184735")

# Tally unannotated cells
n_annotated <- sum(!is.na(reference_obj$celltype))
n_unannotated <- sum(is.na(reference_obj$celltype))

# --- 7c: kNN label transfer for unannotated cells ---

# Process cells that did not receive an author annotation
# (all of GSE184735 + unmatched cells from GSE145926/193782)

unannotated_cells <- which(is.na(reference_obj$celltype))
annotated_cells <- which(!is.na(reference_obj$celltype))

if (length(unannotated_cells) > 0 && length(annotated_cells) > 100) {
  
  if ("harmony" %in% names(reference_obj@reductions)) {
    
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
      
    }
    
    reference_obj$celltype[unannotated_cells] <- knn_labels
    reference_obj$annotation_source[unannotated_cells] <- "kNN_transfer"
    
    print(table(knn_labels))
    
    low_conf <- sum(knn_confidence < 0.5)
    
  } else {
    # Harmony not run -> canonical marker fallback
    
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
    
    print(table(assigned))
  }
  
}

# Final annotation summary

print(table(reference_obj$celltype))

print(table(reference_obj$annotation_source))

# Handle any remaining NA
n_na <- sum(is.na(reference_obj$celltype))
if (n_na > 0) {
  
  reference_obj$celltype[is.na(reference_obj$celltype)] <- "Other"
}

# --- 7d: Macrophage subtype classification ---

mac_cells <- which(reference_obj$celltype == "Macrophage")

# Initialize celltype_fine (same as the coarse celltype)
reference_obj$celltype_fine <- reference_obj$celltype

if (length(mac_cells) > 0) {
  expr_data <- GetAssayData(reference_obj, layer = "data")
  
  avail_alv <- intersect(mac_subtype_markers$Alveolar_Mac, rownames(reference_obj))
  avail_mono <- intersect(mac_subtype_markers$Monocyte_Mac, rownames(reference_obj))
  
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
    
  }
}

# --- 7e: Annotation quality validation with canonical markers ---

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

val_tab <- NULL
for (ct in names(validation_markers)) {
  ct_cells <- which(reference_obj$celltype == ct)
  other_cells <- which(reference_obj$celltype != ct)
  
  if (length(ct_cells) == 0) next
  
  for (marker in validation_markers[[ct]]) {
    if (marker %in% rownames(expr_data)) {
      in_type_mean <- mean(expr_data[marker, ct_cells])
      other_mean <- mean(expr_data[marker, other_cells])
      enrichment <- ifelse(other_mean > 0, in_type_mean / other_mean, Inf)
      status <- ifelse(enrichment > 1.5, "enriched", ifelse(enrichment > 1.0, "weak", "not enriched"))
      val_tab <- rbind(val_tab, data.frame(celltype = ct, marker = marker, enrichment = round(enrichment, 2), status = status))
    }
  }
}
print(val_tab, row.names = FALSE)   # marker enrichment in each annotated cell type

# --- 7f: Annotation summary ---

print(table(reference_obj$celltype))

print(table(reference_obj$annotation_source))

print(table(reference_obj$dataset, reference_obj$celltype))

print(table(reference_obj$celltype_fine))

saveRDS(reference_obj, file.path(output_dir, "BAL_reference_author_annotated.rds"))
