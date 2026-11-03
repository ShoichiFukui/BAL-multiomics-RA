# 12_MOFA_AUC_DeLong.R -- Results paragraph "Multi-omics integration" and Figure 5c/5i: serum IL-6 AUC with bootstrap 95% CI (RA vs
#   sarcoidosis, n = 31), disease Factor 1 AUC with CI and exact Wilcoxon P, paired DeLong test Factor 1 vs serum IL-6 (Figure 5c),
#   methotrexate use by infection status, and paired DeLong test BAL Th17.1 count vs RA-specific Factor 1 for a relative %FVC decline
#   of 10% or more (INBUILD criterion; Figure 5i). The DeLong test for infection (Figure 5f) is computed in 05_RA_only_MOFA2.R.
# Produces: <ART_TABLES>/n31/MOFA_AUC_DeLong.csv.
# Inputs: diseaseMOFA2_factor_scores.csv (04_disease_MOFA2.R), panels/panel_data_pft.csv (07_pulmonary_function_panels.R),
#   planA_data_RA24.csv and results/tables/master_data_with_CT_all.csv (controlled access).
# pROC; 2,000 bootstrap resamples with set.seed(42) before each ci.auc. For the paired DeLong tests both ROC curves are oriented in the
#   same direction (a lower BAL Th17.1 count scored as higher risk).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages(library(pROC))
R0 <- file.path(ART_TABLES, "n31")
fs <- read.csv(file.path(R0, "diseaseMOFA2_factor_scores.csv"), check.names = FALSE)
ms <- read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE)
id <- fs[[1]]
y <- ifelse(is_sarc(id), 0, 1); il6 <- as.numeric(ms$Serum_IL6[match(id, ms$Sample_ID)])
r <- roc(y, il6, quiet = TRUE, direction = "auto"); set.seed(42); ci <- ci.auc(r, method = "bootstrap", boot.n = 2000)

f1 <- fs[[grep("Factor1|^F1$|factor1", names(fs), ignore.case = TRUE)[1]]]
r2 <- roc(y, f1, quiet = TRUE, direction = "auto"); set.seed(42); c2 <- ci.auc(r2, method = "bootstrap", boot.n = 2000)

ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); mtx <- med_at_BAL(ra$Sample_ID, "MTX"); stopifnot(sum(mtx) == 8, sum(mtx[ra$respiratory_infection == 1]) == 0)   # MTX at BAL
print(table(infection = ra$respiratory_infection, MTX = mtx, useNA = "ifany"))

d <- roc.test(r2, r, method = "delong")
pf <- read.csv(file.path(ART_FIGS, "n31", "panels", "panel_data_pft.csv"))
ok2 <- is.finite(pf$inbuild) & is.finite(pf$t171c) & is.finite(pf$F1)
i1 <- roc(pf$inbuild[ok2], -pf$t171c[ok2], quiet = TRUE, direction = "<"); i2 <- roc(pf$inbuild[ok2], pf$F1[ok2], quiet = TRUE, direction = "<")
i2a <- roc(pf$inbuild[ok2], pf$F1[ok2], quiet = TRUE, direction = "auto")
di <- roc.test(i1, i2, method = "delong", paired = TRUE)
res <- data.frame(
  Statistic = c("Serum IL-6 AUC (RA vs sarcoidosis)", "Serum IL-6 AUC 95% CI low", "Serum IL-6 AUC 95% CI high",
                "Disease Factor 1 AUC", "Disease Factor 1 AUC 95% CI low", "Disease Factor 1 AUC 95% CI high",
                "Disease Factor 1 exact Wilcoxon P", "DeLong P, Factor 1 vs serum IL-6 (Figure 5C)",
                "AUC BAL Th17.1 count, INBUILD decline", "AUC RA Factor 1, INBUILD decline (direction auto)",
                "DeLong P, BAL Th17.1 count vs RA Factor 1 (Figure 5I)"),
  Value = c(auc(r), ci[1], ci[3], auc(r2), c2[1], c2[3], wilcox.test(f1[y == 1], f1[y == 0], exact = TRUE)$p.value, d$p.value,
            auc(i1), auc(i2a), di$p.value))
write.csv(res, file.path(R0, "MOFA_AUC_DeLong.csv"), row.names = FALSE)
