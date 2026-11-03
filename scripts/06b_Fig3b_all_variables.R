# 06b_Fig3b_all_variables.R -- Source Data of Figure 3b: every variable compared by subsequent respiratory infection (RA, n = 24), in the
#   order of the Source Data sheet. Figure 3b plots the 12 markers of the 114-marker panel with nominal p < 0.05; the other rows are listed
#   but not plotted: BayesPrism fractions and GSVA scores (derived, outside the FDR family), duplicate clinical-sheet columns (BAL-CD25,
#   CD25, BAL neutrophil CD66b+CD206- %) and constant columns. P: two-sided Wilcoxon rank-sum, exact = TRUE (normal approximation with
#   ties); Cliff's delta (Inf+ vs Inf-); P_adjusted = BH across the 114-marker panel (06), blank outside the panel.
# Produces: common_panel/Fig3b_all_variables.csv.
# Inputs: planA_data_RA24.csv, master_data_with_CT_all.csv, bayesprism_n31 cell fractions (01b), gsva_GO_scores_n31.csv (17),
#   infection_panel.csv (06), results/RA_ILD_Workspace.RData (microbiome diversity).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
R0 <- file.path(ART_TABLES, "n31"); num <- function(x) suppressWarnings(as.numeric(x))
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))
ms <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE)); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
bp <- read.csv(file.path(R0, "bayesprism_n31/celltype_proportions_cellFraction_n31.csv"), row.names = 1)
gv <- read.csv(file.path(R0, "gsva_GO/gsva_GO_scores_n31.csv"))
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e); dv <- e$microbiome_results$alpha_diversity; rm(e)
ip <- read.csv(file.path(R0, "common_panel/infection_panel.csv"))
cd4 <- num(ms[["BAL_lymphocytes_per_mL"]]) * num(ms$CD4) / 100
  VARS <- c("NK" = "Deconvolution", "DC" = "Deconvolution", "Mast" = "Deconvolution", "BALF_MCP3" = "BALF_Cytokine",
  "BALF_Th17.1_CCR6pos_CXCR3pos" = "FCM", "CD3" = "FCM", "Other" = "Deconvolution", "BALF_Th17_CCR6pos_CXCR3neg" = "FCM",
  "BALF_Treg_CD127neg_CD25pos" = "FCM", "BALF_Th2_CCR6neg_CXCR3neg" = "FCM", "Macrophage" = "Deconvolution", "CD56" = "FCM", "BALF_GMCSF" =
  "BALF_Cytokine", "Neutrophil" = "Deconvolution", "BALF_IL2" = "BALF_Cytokine", "CD22" = "FCM", "BALF_Th1_CCR6neg_CXCR3pos" = "FCM", "CD16" =
  "FCM", "BALF_IL15" = "BALF_Cytokine", "Serum_Eotaxin" = "Serum_Cytokine", "Serum_IL6" = "Serum_Cytokine", "BALF_IL1a" = "BALF_Cytokine",
  "PB_Th2_CCR6neg_CXCR3neg" = "FCM", "BALF_IL12p70" = "BALF_Cytokine", "BALF_MDC" = "BALF_Cytokine", "Serum_IP10" = "Serum_Cytokine",
  "BALF_Activated_Macrophage_CD86pos_CD14pos" = "FCM", "BALF_sCD40L" = "BALF_Cytokine", "PB_Neutrophil_CD66bpos_CD206neg" = "FCM",
  "Serum_TNFa" = "Serum_Cytokine", "BALF_Macrophage_Percent" = "FCM", "PB_Th17.1_CCR6pos_CXCR3pos" = "FCM", "PB_Th17_CCR6pos_CXCR3neg" =
  "FCM", "T_cell" = "Deconvolution", "BALF_IL12p40" = "BALF_Cytokine", "BALF_IL17A" = "BALF_Cytokine", "BALF_IL6" = "BALF_Cytokine",
  "Serum_IL10" = "Serum_Cytokine", "PB_Activated_Macrophage_CD86pos_CD14pos" = "FCM", "BALF_FGF2" = "BALF_Cytokine", "BALF_MCP1" =
  "BALF_Cytokine", "BALF_Lymphocyte_Percent" = "FCM", "Serum_MIP1b" = "Serum_Cytokine", "BALF_IL13" = "BALF_Cytokine", "BALF_IL1b" =
  "BALF_Cytokine", "Serum_GCSF" = "Serum_Cytokine", "Serum_sCD40L" = "Serum_Cytokine", "BALF_MIP1a" = "BALF_Cytokine", "BAL-CD25" = "FCM",
  "CD25" = "FCM", "BALF_TNFa" = "BALF_Cytokine", "Serum_TGFa" = "Serum_Cytokine", "BALF_IL8" = "BALF_Cytokine", "Serum_VEGFA" =
  "Serum_Cytokine", "B_cell" = "Deconvolution", "BALF_Fractalkine" = "BALF_Cytokine", "PB_Macrophage_CD14pos_CD16neg" = "FCM", "BALF_TGFa" =
  "BALF_Cytokine", "BALF_IL4" = "BALF_Cytokine", "Serum_GMCSF" = "Serum_Cytokine", "Serum_IL5" = "Serum_Cytokine", "Serum_IL1a" =
  "Serum_Cytokine", "Serum_IL1b" = "Serum_Cytokine", "Serum_IL1RA" = "Serum_Cytokine", "Serum_IL2" = "Serum_Cytokine", "Serum_IL17F" =
  "Serum_Cytokine", "BALF_GCSF" = "BALF_Cytokine", "Serum_FGF2" = "Serum_Cytokine", "Serum_IL13" = "Serum_Cytokine", "Serum_IL17A" =
  "Serum_Cytokine", "Serum_IL7" = "Serum_Cytokine", "Plasma" = "Deconvolution", "BALF_IL5" = "BALF_Cytokine", "BALF_IP10" = "BALF_Cytokine",
  "BALF_Neutrophil_CD66b_pos_CD206_neg(%)" = "FCM", "BALF_Eotaxin" = "BALF_Cytokine", "Serum_IL12p70" = "Serum_Cytokine", "Serum_IFNa2" =
  "Serum_Cytokine", "BALF_IFNg" = "BALF_Cytokine", "Serum_IFNg" = "Serum_Cytokine", "Serum_IL12p40" = "Serum_Cytokine", "Serum_GROa" =
  "Serum_Cytokine", "Serum_IL15" = "Serum_Cytokine", "Serum_MDC" = "Serum_Cytokine", "Epithelial" = "Deconvolution", "Serum_IL22" =
  "Serum_Cytokine", "BALF_Macrophage_CD14pos_CD16neg" = "FCM", "Shannon" = "Microbiome", "BALF_IL10" = "BALF_Cytokine", "Serum_MIP1a" =
  "Serum_Cytokine", "BALF_IL22" = "BALF_Cytokine", "CD4" = "FCM", "BALF_PDGFAA" = "BALF_Cytokine", "Serum_FLT3L" = "Serum_Cytokine",
  "Serum_TNFb" = "Serum_Cytokine", "Serum_Fractalkine" = "Serum_Cytokine", "Serum_IL4" = "Serum_Cytokine", "BALF_IL7" = "BALF_Cytokine",
  "Serum_MCP3" = "Serum_Cytokine", "CD8" = "FCM", "BALF_Neutrophil_CD66bpos_CD206neg" = "FCM", "Serum_MCP1" = "Serum_Cytokine", "BALF_IL27" =
  "BALF_Cytokine", "BALF_TNFb" = "BALF_Cytokine", "PB_Th1_CCR6neg_CXCR3pos" = "FCM", "Serum_IL27" = "Serum_Cytokine",
  "BALF_Neutrophil_Percent" = "FCM", "BALF_IL18" = "BALF_Cytokine", "BALF_GROa" = "BALF_Cytokine", "Serum_EGF" = "Serum_Cytokine",
  "Serum_IL18" = "Serum_Cytokine", "PB_Treg_CD127neg_CD25pos" = "FCM", "Simpson" = "Microbiome", "BALF_IFNa2" = "BALF_Cytokine",
  "BALF_Eosinophil_Percent" = "FCM", "BALF_IL1RA" = "BALF_Cytokine", "BALF_VEGFA" = "BALF_Cytokine", "BALF_EGF" = "BALF_Cytokine", "Serum_IL8"
  = "Serum_Cytokine", "Serum_PDGFAA" = "Serum_Cytokine", "BALF_MIP1b" = "BALF_Cytokine", "CD4_CD8_Ratio" = "FCM", "BALF_FLT3L" =
  "BALF_Cytokine", "BALF_IL17F" = "BALF_Cytokine", "BALF_Basophil_Percent" = "FCM", "BALF_Plasma_Cell_Percent" = "FCM",
  "BALF_Other_Cell_1_Percent" = "FCM", "BALF_Other_Cell_2_Percent" = "FCM", "BAL_Th17.1_count" = "BAL_CD4_subset_count", "BAL_Th17_count" =
  "BAL_CD4_subset_count", "BAL_Th1_count" = "BAL_CD4_subset_count", "BAL_Th2_count" = "BAL_CD4_subset_count", "BAL_Treg_count" =
  "BAL_CD4_subset_count", "GSVA_GO_OXPHOS" = "GSVA_gene_set_score", "GSVA_GO_Surfactant" = "GSVA_gene_set_score")
stopifnot(length(VARS) == 135)
value_of <- function(v, k) {
  if (k == "Deconvolution") return(100 * bp[ra$Sample_ID, v])
  if (k == "GSVA_gene_set_score") return(gv[[sub("GSVA_GO_", "", v)]][match(ra$Sample_ID, gv$Sample_ID)])
  if (k == "Microbiome") return(dv[[v]][match(ra$Sample_ID, dv$Sample_ID)])
  if (k == "BAL_CD4_subset_count") { s <- sub("^BAL_(.*)_count$", "\\1", v); return(cd4 * num(ms[[grep(paste0("^BALF_", gsub(".", "\\.", s, fixed = TRUE), "_"), names(ms), value = TRUE)[1]]]) / 100) }
  col <- c(`BAL-CD25` = "BAL.CD25", `BALF_Neutrophil_CD66b_pos_CD206_neg(%)` = "BAL_neutrophil_CD66b_CD206neg_pct_sheet")[v]
  num(ms[[if (is.na(col)) v else col]])
}
panel_id <- function(v, k) if (k == "Microbiome") paste0("Microbiome_", v) else if (v %in% c("CD3", "CD4", "CD8", "CD22", "CD16", "CD56", "CD4_CD8_Ratio")) paste0("BAL_", v) else v
mag <- function(d) { a <- abs(d); if (a < 0.147) "negligible" else if (a < 0.33) "small" else if (a < 0.474) "medium" else "large" }
inf <- ra$respiratory_infection; stopifnot(sum(inf) == 6)
out <- do.call(rbind, lapply(names(VARS), function(v) { k <- VARS[[v]]; x <- value_of(v, k); stopifnot(length(x) == 24)
  a <- x[is.finite(x) & inf == 1]; b <- x[is.finite(x) & inf == 0]; cd <- mean(outer(a, b, function(p, q) sign(p - q)))
  p <- if (length(unique(c(a, b))) > 1) suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value else NA
  data.frame(Variable = v, Category = k, Median_InfPos = median(a), Median_InfNeg = median(b), N_pos = length(a), N_neg = length(b),
             P_value = p, Cliff_delta = cd, Magnitude = mag(cd), P_adjusted = ip$BH[match(panel_id(v, k), ip$Variable)]) }))
out$Category[grepl("^BALF_.*_Percent$", out$Variable)] <- "BAL_differential_count"   # clinical-laboratory differential counts, as in Supplementary Data 2
stopifnot(sum(!is.na(out$P_adjusted)) == 114, sum(out$P_value < 0.05 & !is.na(out$P_adjusted), na.rm = TRUE) == 12)   # the 12 markers of Figure 3b
write.csv(out, file.path(R0, "common_panel/Fig3b_all_variables.csv"), row.names = FALSE)
