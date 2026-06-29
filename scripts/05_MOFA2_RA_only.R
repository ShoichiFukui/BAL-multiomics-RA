#!/usr/bin/env Rscript
# 05_MOFA2_RA_only.R — RA-only MOFA2 for infection + CT progression
library(reticulate)
# Set MOFA2_PYTHON to the python interpreter that has 'mofapy2' installed.
mofa_python <- Sys.getenv("MOFA2_PYTHON", unset = "")
if (nzchar(mofa_python)) use_python(mofa_python, required = TRUE)
library(readxl); library(tidyverse); library(ggplot2); library(pROC)
library(MOFA2); library(patchwork)
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
load("results/CT_Multiomics_Results.RData")

# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)

bp_cell <- read.csv("output_v3_BayesPrism/celltype_proportions_cellFraction.csv", row.names=1)
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet="BALF_analysis")
colnames(fcm)[1] <- "Sample_ID"
bp_fcm_common <- intersect(rownames(bp_cell), fcm$Sample_ID)
bp_fcm_common <- bp_fcm_common[!bp_fcm_common %in% EXCL]
ra_ids <- bp_fcm_common[grepl("^RA[0-9]", bp_fcm_common)]
ct_u <- read.csv("results/tables/master_data_with_CT_all.csv")
ct_u <- ct_u[!duplicated(ct_u$Sample_ID),]

cat("RA samples:", length(ra_ids), "\n")

# ── Build views ──
expr_ra <- expr_matrix[, ra_ids]
top_genes <- names(sort(apply(expr_ra, 1, var), decreasing=TRUE))[1:30]
expr_view <- t(expr_ra[top_genes, ])

balf_cyto_cols <- grep("^BALF_", colnames(master_data), value=TRUE)
balf_cyto_cols <- balf_cyto_cols[!grepl("Th17|Treg|Th1|Th2|Macro|Neutro|CD|Lymph|Eosino|Baso|Plasma|Other|Percent", balf_cyto_cols)]
balf_cyto_cols <- balf_cyto_cols[sapply(balf_cyto_cols, function(x) is.numeric(master_data[[x]]))]
balf_view <- as.matrix(master_data[match(ra_ids, master_data$Sample_ID), balf_cyto_cols])
rownames(balf_view) <- ra_ids
balf_view <- balf_view[, colSums(!is.na(balf_view)) >= 18]

serum_cols <- grep("^Serum_", colnames(master_data), value=TRUE)
serum_view <- as.matrix(master_data[match(ra_ids, master_data$Sample_ID), serum_cols])
rownames(serum_view) <- ra_ids
serum_view <- serum_view[, colSums(!is.na(serum_view)) >= 18]

balf_fcm_cols <- intersect(c("BALF_Th17.1_CCR6pos_CXCR3pos","BALF_Th17_CCR6pos_CXCR3neg",
  "BALF_Treg_CD127neg_CD25pos","BALF_Th1_CCR6neg_CXCR3pos","BALF_Th2_CCR6neg_CXCR3neg",
  "BALF_Macrophage_Percent","BALF_Lymphocyte_Percent","BALF_Neutrophil_Percent"), colnames(master_data))
balf_fcm <- as.matrix(master_data[match(ra_ids, master_data$Sample_ID), balf_fcm_cols])
rownames(balf_fcm) <- ra_ids

pb_fcm_cols <- intersect(c("PB_Th17.1_CCR6pos_CXCR3pos","PB_Th17_CCR6pos_CXCR3neg",
  "PB_Treg_CD127neg_CD25pos","PB_Th1_CCR6neg_CXCR3pos","PB_Th2_CCR6neg_CXCR3neg"), colnames(master_data))
pb_fcm <- as.matrix(master_data[match(ra_ids, master_data$Sample_ID), pb_fcm_cols])
rownames(pb_fcm) <- ra_ids

data_list <- list(Expression=t(expr_view), BALF_Cytokine=t(balf_view), Serum_Cytokine=t(serum_view),
                  BALF_FCM=t(balf_fcm), PB_FCM=t(pb_fcm))

cat("Views:", paste(names(data_list), sapply(data_list, ncol), sep="=", collapse=", "), "\n")
cat("Total features:", sum(sapply(data_list, ncol)), "\n")

# Scale
for (v in names(data_list)) {
  m <- data_list[[v]]
  for (j in 1:ncol(m)) {
    s <- sd(m[,j], na.rm=TRUE)
    if (!is.na(s) && s > 0) m[,j] <- scale(m[,j])
  }
  data_list[[v]] <- m
}

# ── Train MOFA2 ──
mofa_ra <- create_mofa(data_list)
mo <- get_default_model_options(mofa_ra); mo$num_factors <- 3
to <- get_default_training_options(mofa_ra)
to$maxiter <- 1000; to$seed <- 42; to$verbose <- FALSE; to$convergence_mode <- "medium"
mofa_ra <- prepare_mofa(mofa_ra, model_options=mo, training_options=to)

cat("\nTraining RA-only MOFA2...\n")
mofa_trained <- run_mofa(mofa_ra, outfile="results/mofa2_RA_only.hdf5", use_basilisk=FALSE)

ra_factors <- get_factors(mofa_trained)[[1]]
ra_r2 <- get_variance_explained(mofa_trained)
ra_weights <- get_weights(mofa_trained)

cat("\nR2 per factor:\n"); print(round(ra_r2$r2_per_factor[[1]], 2))
cat("\nR2 total:\n"); print(round(ra_r2$r2_total[[1]], 2))

# ── Associations ──
inf_st <- master_data$respiratory_infection[match(ra_ids, master_data$Sample_ID)]
dhl <- as.numeric(ct_u$CT_Delta_Healthy_Lung_pct[match(ra_ids, ct_u$Sample_ID)])
prog_st <- ifelse(dhl < -10, 1, 0)

cat("\n=== Factor-outcome associations ===\n")
assoc <- data.frame()
for (f in 1:3) {
  fv <- ra_factors[, f]
  # Infection
  ok_i <- !is.na(inf_st) & !is.na(fv)
  wt_i <- wilcox.test(fv[ok_i & inf_st==1], fv[ok_i & inf_st==0], exact=TRUE)
  r_i <- tryCatch(roc(inf_st[ok_i], fv[ok_i], quiet=TRUE, direction="auto"), error=function(e) NULL)
  # Progression
  ok_p <- !is.na(prog_st) & !is.na(fv)
  wt_p <- wilcox.test(fv[ok_p & prog_st==1], fv[ok_p & prog_st==0], exact=TRUE)
  r_p <- tryCatch(roc(prog_st[ok_p], fv[ok_p], quiet=TRUE, direction="auto"), error=function(e) NULL)
  # CT rho
  sp <- cor.test(fv, dhl, method="spearman", exact=TRUE)

  assoc <- rbind(assoc, data.frame(
    Factor=f,
    Inf_p=wt_i$p.value, Inf_AUC=ifelse(!is.null(r_i), as.numeric(auc(r_i)), NA),
    Prog_p=wt_p$p.value, Prog_AUC=ifelse(!is.null(r_p), as.numeric(auc(r_p)), NA),
    CT_rho=sp$estimate, CT_p=sp$p.value))

  cat(sprintf("  Factor %d: Inf p=%.4f AUC=%.3f | Prog p=%.4f AUC=%.3f | CT rho=%.3f p=%.4f\n",
              f, wt_i$p.value, ifelse(!is.null(r_i), as.numeric(auc(r_i)), NA),
              wt_p$p.value, ifelse(!is.null(r_p), as.numeric(auc(r_p)), NA),
              sp$estimate, sp$p.value))
}

# Top weights
cat("\n=== Top weights ===\n")
for (f in 1:3) {
  fw <- data.frame()
  for (v in names(ra_weights)) {
    wv <- ra_weights[[v]][, f]
    fw <- rbind(fw, data.frame(View=v, Feature=names(wv), Weight=as.numeric(wv)))
  }
  fw <- fw %>% arrange(desc(abs(Weight))) %>% head(8)
  cat(sprintf("\nFactor %d:\n", f))
  for (i in 1:nrow(fw)) cat(sprintf("  %s: %s (%.3f)\n", fw$View[i], fw$Feature[i], fw$Weight[i]))
}

# Full factor weights -> Supplementary Data 9 (RA-specific MOFA2: 5 views, n = 24)
ra_w_long <- data.frame()
for (f in 1:ncol(ra_weights[[1]])) {
  for (v in names(ra_weights)) {
    wv <- ra_weights[[v]][, f]
    ra_w_long <- rbind(ra_w_long, data.frame(
      View = sub("_FACS$", "_FCM", v), Factor = f, Feature = names(wv), Weight = round(as.numeric(wv), 5)))
  }
}
write.csv(ra_w_long, "results/tables/MOFA2_RAonly_weights.csv", row.names = FALSE)
cat(sprintf("\u2713 Wrote MOFA2_RAonly_weights.csv (Suppl. Data 9): %d rows\n", nrow(ra_w_long)))

# Save
save(mofa_trained, ra_factors, ra_r2, ra_weights, assoc,
     file="results/MOFA2_RA_only_Results.RData")
cat("\n✓ Saved.\n")
