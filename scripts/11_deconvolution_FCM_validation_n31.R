# 11_deconvolution_FCM_validation_n31.R -- Figure 1g-i (deconvolution vs the BAL differential counts, columns BALF_*_Percent of the
#   FCM workbook: Spearman, Lin's CCC), Supplementary Fig. 1
#   (Bland-Altman bias and limits of agreement) and Figure 2f (BayesPrism plasma cell and B cell fractions, RA vs sarcoidosis; the
#   plasma cell p = 0.061 and B cell p = 0.30 of the Results), on the n = 31 BayesPrism re-run (01b_BayesPrism_n31.R).
# Produces: <ART_TABLES>/n31/bayesprism_n31/Fig2f_bp_cell_plasma_Bcell.csv and Fig1ghi_SFig1_deconv_vs_differential.csv.
# Inputs: bayesprism_n31/celltype_proportions_cellFraction_n31.csv (01), FCM_integrated_data_transformed.xlsx (controlled access).
# Lymphocytes = T cell + NK + B cell + plasma cell fractions. Wilcoxon and Spearman exact (normal approximation with ties).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(readxl); library(DescTools) })
cfg <- set_cohort("n31"); BPD <- file.path(ART_TABLES, "n31", "bayesprism_n31")

bp_cell <- read.csv(file.path(BPD, "celltype_proportions_cellFraction_n31.csv"), row.names = 1)
bp_cell <- bp_cell[!rownames(bp_cell) %in% cfg$EXCL, ]
check_n(rownames(bp_cell), cfg, "bp_cell")
ga <- sample_group(rownames(bp_cell))

# ---- (3) Fig 2f -------------------------------------------------------------
fig2f <- do.call(rbind, lapply(c("Plasma", "B_cell"), function(ct) {
  v <- bp_cell[[ct]] * 100
  wt <- wilcox.test(v[ga == "RA"], v[ga == "Sarcoidosis"], exact = TRUE)
  data.frame(CellType = ct, n_RA = sum(ga == "RA"), n_Sarc = sum(ga == "Sarcoidosis"),
             median_RA_pct = median(v[ga == "RA"]), median_Sarc_pct = median(v[ga == "Sarcoidosis"]),
             W = unname(wt$statistic), p = wt$p.value)
}))
write.csv(fig2f, file.path(BPD, "Fig2f_bp_cell_plasma_Bcell.csv"), row.names = FALSE)

# ---- (4) Fig 1g-i + SFig 1a-c --------------------------------------------------
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet = "BALF_analysis")
colnames(fcm)[1] <- "Sample_ID"
cf <- intersect(rownames(bp_cell), fcm$Sample_ID)
fi <- match(cf, fcm$Sample_ID)
lc <- intersect(c("T_cell", "NK", "B_cell", "Plasma"), colnames(bp_cell))

sets <- list(
  Macrophage = list(B = bp_cell[cf, "Macrophage"] * 100, F = as.numeric(fcm$BALF_Macrophage_Percent[fi])),
  Lymphocyte = list(B = rowSums(bp_cell[cf, lc]) * 100,  F = as.numeric(fcm$BALF_Lymphocyte_Percent[fi])),
  Neutrophil = list(B = bp_cell[cf, "Neutrophil"] * 100, F = as.numeric(fcm$BALF_Neutrophil_Percent[fi])))

val <- do.call(rbind, lapply(names(sets), function(nm) {
  B <- sets[[nm]]$B; F <- sets[[nm]]$F
  ok <- !is.na(B) & !is.na(F); B <- B[ok]; F <- F[ok]
  sp <- suppressWarnings(cor.test(B, F, method = "spearman", exact = TRUE))
  cc <- DescTools::CCC(B, F)$rho.c
  d <- B - F; bias <- mean(d); s <- sd(d)
  data.frame(CellType = nm, N = length(B), Spearman_rho = unname(sp$estimate),
             Spearman_p = sp$p.value, CCC = cc$est, CCC_lwr = cc$lwr.ci, CCC_upr = cc$upr.ci,
             BA_bias = bias, BA_LoA_lower = bias - 1.96 * s, BA_LoA_upper = bias + 1.96 * s)
}))
write.csv(val, file.path(BPD, "Fig1ghi_SFig1_deconv_vs_differential.csv"), row.names = FALSE)

