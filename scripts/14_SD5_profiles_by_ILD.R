# 14_SD5_profiles_by_ILD.R -- Supplementary Data 5 (BAL cellular, BAL and serum cytokine and PB flow cytometry profiles by baseline ILD
#   status, RA: ILD+ n = 11, ILD- n = 13): median by ILD status, exact Wilcoxon P (normal approximation with ties), BH adjustment over the
#   100 rows that remain after removal of the two CD66b+CD16+ neutrophil variables (measured in 11 of 31 patients).
# Produces: <ART_TABLES>/n31/supp_tables/SD5_profiles_by_ILD.csv (100 rows; P and BH of Supplementary Data 5).
# Inputs: planA_data_RA24.csv and results/tables/master_data_with_CT_all.csv (controlled access).
# The 102 variables are the rows of the table before the removal (the variables not detected or constant in all samples were not listed).
# The variable list is written out here (computation as run for the manuscript).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
SUPP <- file.path(ART_TABLES, "n31", "supp_tables"); dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
DROP <- c("BALF_Neutrophil_CD66bpos_CD16pos", "PB_Neutrophil_CD66bpos_CD16pos")
VARS <- c("BALF_Treg_CD127neg_CD25pos", "BALF_Macrophage_CD14pos_CD16neg", "BALF_IL1RA", "Serum_FGF2", "Serum_GMCSF", "Serum_IL1b", "BALF_IL18",
  "Serum_IL17A", "Serum_PDGFAA", "Serum_IL17F", "Serum_TGFa", "Serum_IL2", "Serum_Fractalkine", "Serum_IL12p70", "Serum_IL1a",
  "PB_Neutrophil_CD66bpos_CD206neg", "Serum_IL13", "Serum_IL12p40", "Serum_IP10", "BALF_GROa", "BALF_MCP3", "Serum_IL15", "Serum_IL7",
  "BALF_GCSF", "BALF_MDC", "Serum_FLT3L", "Serum_MCP3", "BALF_IFNa2", "Serum_EGF", "BALF_IL12p40", "Serum_GROa", "BALF_Neutrophil_Percent",
  "BALF_VEGFA", "Serum_IL22", "BALF_IL1a", "BALF_PDGFAA", "Serum_IL4", "Serum_IL1RA", "Serum_TNFb", "BALF_Th17_CCR6pos_CXCR3neg",
  "BALF_MIP1b", "BALF_IL4", "PB_Th1_CCR6neg_CXCR3pos", "Serum_IL5", "BALF_IFNg", "PB_Macrophage_CD14pos_CD16neg", "Serum_IFNa2", "BALF_IL6",
  "BALF_Macrophage_Percent", "BALF_Activated_Macrophage_CD86pos_CD14pos", "PB_Th17.1_CCR6pos_CXCR3pos", "Serum_Eotaxin", "Serum_MCP1",
  "BALF_Neutrophil_CD66bpos_CD206neg", "Serum_GCSF", "BALF_TNFb", "Serum_TNFa", "Serum_sCD40L", "BALF_Lymphocyte_Percent", "Serum_MIP1a",
  "BALF_IL10", "Serum_IL8", "BALF_GMCSF", "BALF_Eosinophil_Percent", "BALF_EGF", "PB_Treg_CD127neg_CD25pos", "Serum_MDC", "BALF_IL1b",
  "PB_Th17_CCR6pos_CXCR3neg", "Serum_IL27", "BALF_IL15", "BALF_IL7", "BALF_MCP1", "BALF_IL8", "BALF_Th1_CCR6neg_CXCR3pos", "Serum_VEGFA",
  "BALF_IL12p70", "BALF_MIP1a", "BALF_TNFa", "BALF_sCD40L", "BALF_Fractalkine", "BALF_IL27", "Serum_IL6", "Serum_IFNg", "Serum_IL18",
  "BALF_IL13", "BALF_IL22", "BALF_IP10", "BALF_TGFa", "BALF_Neutrophil_CD66bpos_CD16pos", "BALF_IL17A", "BALF_Th17.1_CCR6pos_CXCR3pos",
  "PB_Th2_CCR6neg_CXCR3neg", "BALF_Eotaxin", "Serum_MIP1b", "PB_Activated_Macrophage_CD86pos_CD14pos", "Serum_IL10",
  "PB_Neutrophil_CD66bpos_CD16pos", "BALF_Th2_CCR6neg_CXCR3neg", "BALF_FGF2", "BALF_IL2", "BALF_IL5")
stopifnot(length(VARS) == 102, all(DROP %in% VARS))
ra <- read.csv(file.path(ART_TABLES, "n31", "planA/planA_data_RA24.csv")); ms <- read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
ild <- as.integer(ra$ILD_Group == "ILD_Positive"); stopifnot(sum(ild == 1) == 11, sum(ild == 0) == 13)
rec <- do.call(rbind, lapply(VARS, function(v) { x <- suppressWarnings(as.numeric(ms[[v]])); a <- x[ild == 1 & !is.na(x)]; b <- x[ild == 0 & !is.na(x)]
  data.frame(Variable = v, median_pos = median(a), median_neg = median(b), n_pos = length(a), n_neg = length(b), P = suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value) }))
new <- rec[!rec$Variable %in% DROP, ]; new$BH <- p.adjust(new$P, "BH"); stopifnot(nrow(new) == 100)

write.csv(new, file.path(SUPP, "SD5_profiles_by_ILD.csv"), row.names = FALSE)
