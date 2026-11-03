# 01b_BayesPrism_n31.R -- BayesPrism deconvolution of the 31 analysed samples (Fig. 1f-i, 2f; Suppl. Fig. 1; BayesPrism rows of Suppl. Data 2,
#   in its FDR family, and of Suppl. Data 3, outside the 114-marker FDR family).
# Produces: bayesprism_n31/celltype_proportions_BayesPrism_n31.csv (RNA fraction) and celltype_proportions_cellFraction_n31.csv (cell fraction).
# Inputs: output_BayesPrism/BAL_reference_author_annotated.rds (scRNA-seq reference built by scripts/01a_reference_atlas.R)
#   and RNAseq/RawCount_heatmap_table.csv (bulk counts; also in GEO GSE329884).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
suppressPackageStartupMessages({ library(Seurat); library(BayesPrism); library(org.Hs.eg.db) })
source("scripts/00_cohort_config.R")
OD <- file.path(ART_TABLES, "n31", "bayesprism_n31"); dir.create(OD, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
reference_obj <- readRDS("output_BayesPrism/BAL_reference_author_annotated.rds")
# ---- bulk counts: Ensembl -> gene symbol, duplicated symbols summed ----
bulk_raw <- read.csv("RNAseq/RawCount_heatmap_table.csv", row.names = 1)
if ("Cluster_ID" %in% colnames(bulk_raw)) bulk_raw <- bulk_raw[, -which(colnames(bulk_raw) == "Cluster_ID")]
bulk_matrix <- as.matrix(bulk_raw)
ensembl_ids <- gsub("\\..*", "", rownames(bulk_matrix))
mapped <- AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids, columns = "SYMBOL", keytype = "ENSEMBL")
mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]; mapped <- mapped[!duplicated(mapped$ENSEMBL), ]
gene_mapping <- data.frame(ensembl_gene_id = mapped$ENSEMBL, external_gene_name = mapped$SYMBOL)
new_names <- gene_mapping$external_gene_name[match(ensembl_ids, gene_mapping$ensembl_gene_id)]
valid <- !is.na(new_names) & new_names != ""
bulk_matrix <- bulk_matrix[valid, ]; rownames(bulk_matrix) <- new_names[valid]
if (any(duplicated(rownames(bulk_matrix)))) {
  bulk_df <- as.data.frame(bulk_matrix); bulk_df$gene <- rownames(bulk_matrix)
  bulk_agg <- aggregate(. ~ gene, data = bulk_df, FUN = sum); rownames(bulk_agg) <- bulk_agg$gene; bulk_agg$gene <- NULL
  bulk_matrix <- as.matrix(bulk_agg)
}
# ---- the 31 analysed samples ----
keep <- setdiff(colnames(bulk_matrix), c(EXCL_ALL))
stopifnot(length(keep) == 31, all(keep %in% COHORT_TAB$Sample_ID[COHORT_TAB$status == "analysed"]))   # every bulk column is a registered analysed sample
check_n(keep, set_cohort("n31"), "bulk for BayesPrism")
bulk_matrix <- bulk_matrix[, keep]

# ---- BayesPrism ----
sc_dat <- t(as.matrix(GetAssayData(reference_obj, layer = "counts")))
bk_dat <- t(bulk_matrix)
cell_type_labels <- reference_obj$celltype; cell_state_labels <- reference_obj$celltype_fine
g <- grep("^RP[SL]|^MT-|^MTRNR", colnames(sc_dat), value = TRUE); if (length(g)) sc_dat <- sc_dat[, setdiff(colnames(sc_dat), g)]
g <- grep("^RP[SL]|^MT-|^MTRNR", colnames(bk_dat), value = TRUE); if (length(g)) bk_dat <- bk_dat[, setdiff(colnames(bk_dat), g)]
sc_dat_filtered <- cleanup.genes(input = sc_dat, input.type = "count.matrix", species = "hs",
                                 gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"), exp.cells = 5)
rm(sc_dat); invisible(gc())
myPrism <- new.prism(reference = sc_dat_filtered, mixture = bk_dat, input.type = "count.matrix",
                     cell.type.labels = cell_type_labels, cell.state.labels = cell_state_labels, key = NULL,
                     outlier.cut = 0.01, outlier.fraction = 0.1)
bp_result <- run.prism(prism = myPrism, n.cores = parallel::detectCores() - 1, update.gibbs = TRUE)
theta <- get.fraction(bp = bp_result, which.theta = "final", state.or.type = "type")
# ---- RNA fraction -> cell fraction (median RNA content per cell type, relative to T cells) ----
rna_content <- tapply(reference_obj$nCount_RNA, reference_obj$celltype, median)
baseline <- if ("T_cell" %in% names(rna_content)) rna_content["T_cell"] else median(rna_content)
rna_ratio <- rna_content / as.numeric(baseline)
theta_cell <- theta
for (i in 1:nrow(theta)) {
  cf <- sapply(colnames(theta), function(ct) if (ct %in% names(rna_ratio) && rna_ratio[ct] > 0) theta[i, ct] / rna_ratio[ct] else theta[i, ct])
  if (sum(cf) > 0) theta_cell[i, ] <- cf / sum(cf)
}
write.csv(theta, file.path(OD, "celltype_proportions_BayesPrism_n31.csv"))
write.csv(theta_cell, file.path(OD, "celltype_proportions_cellFraction_n31.csv"))

