# 10_supplementary_table_values.R -- values of Supplementary Table 4 (ST4_BAL_CD4_counts.csv: BAL CD4+ subset counts by group and correlation with the
#   change in %FVC), Supplementary Table 6 (ST6_clinical_vs_multiomics.csv: clinical vs multi-omics markers) and Supplementary Table 7 (ST7_partial_correlations.csv:
#   partial Spearman correlations of the BAL Th17.1 count, one covariate at a time: ILD score, baseline %FVC, interval, infection, MTX,
#   glucocorticoid, BAL lymphocyte count; medication at BAL), and the per-patient Source Data of Figures 4 and 5 (SourceData_Fig4_Fig5_pft.csv).
# Produces: supp_tables/ST4_*, ST6_*, ST7_*, SourceData_Fig4_Fig5_pft.csv.
# Inputs: planA_data_RA24.csv, master_data_with_CT_all.csv, sarc_BAL_WBC.csv, pft_series_by_sample.tsv, 05 output (RA Factor 1),
#   17 GSVA scores (the OXPHOS row of Supplementary Table 6 is the GSVA GO OXPHOS score).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages(library(pROC)); set.seed(42)
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "supp_tables"); dir.create(OD, recursive = TRUE, showWarnings = FALSE); num <- function(x) suppressWarnings(as.numeric(x))
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); ms <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE))
ms <- ms[!ms$Sample_ID %in% EXCL_ALL, ]; g <- ifelse(is_sarc(ms$Sample_ID), "Sarc", "RA")
sw <- read.csv(file.path(R0, "sarc_BAL_WBC.csv"))
w <- num(ms$BAL.WBC); w[g == "Sarc"] <- num(sw$BALWBC[match(ms$Sample_ID[g == "Sarc"], sw$Sample_ID)])
lym <- w * 1e6 * num(ms$BALF_Lymphocyte_Percent) / 100; cd4 <- lym * num(ms$CD4) / 100
S <- c(Th17.1 = "BALF_Th17.1_CCR6pos_CXCR3pos", Th17 = "BALF_Th17_CCR6pos_CXCR3neg", Th1 = "BALF_Th1_CCR6neg_CXCR3pos", Th2 = "BALF_Th2_CCR6neg_CXCR3neg", Treg = "BALF_Treg_CD127neg_CD25pos")
C <- c(list(`BAL lymphocytes` = lym, `BAL CD4+ T cells` = cd4), setNames(lapply(S, function(s) cd4 * num(ms[[s]]) / 100), paste("BAL", names(S))))
i <- match(ra$Sample_ID, ms$Sample_ID); y <- ra$d_pFVC; inf <- ra$respiratory_infection
sp <- function(a, b) { ok <- is.finite(a) & is.finite(b); r <- suppressWarnings(cor.test(a[ok], b[ok], method = "spearman", exact = TRUE)); c(rho = unname(r$estimate), p = r$p.value, n = sum(ok)) }
wt <- function(a, b) suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value
fmt <- function(v) sprintf("%s [%s–%s]", format(signif(median(v), 3), big.mark = ","), format(signif(quantile(v, .25), 3), big.mark = ","), format(signif(quantile(v, .75), 3), big.mark = ","))
t9 <- do.call(rbind, lapply(names(C), function(n) { x <- C[[n]]; a <- x[g == "RA"]; b <- x[g == "Sarc"]; xr <- x[i]; s <- sp(xr, y)
  data.frame(Variable = n, RA = fmt(a), Sarcoidosis = fmt(b), P_RA_vs_Sarc = wt(a, b), Inf_pos = fmt(xr[inf == 1]), Inf_neg = fmt(xr[inf == 0]), P_inf = wt(xr[inf == 1], xr[inf == 0]), N_dFVC = s[["n"]], rho_dFVC = s[["rho"]], P_dFVC = s[["p"]]) }))
t9$BH_RA_vs_Sarc <- p.adjust(t9$P_RA_vs_Sarc, "BH"); t9$BH_inf <- p.adjust(t9$P_inf, "BH"); t9$BH_dFVC <- p.adjust(t9$P_dFVC, "BH"); write.csv(t9, file.path(OD, "ST4_BAL_CD4_counts.csv"), row.names = FALSE)
rk <- function(v) rank(v, na.last = "keep"); pc <- function(x, y, z) { ok <- is.finite(x) & is.finite(y) & is.finite(z); r <- cor(resid(lm(rk(x[ok]) ~ rk(z[ok]))), resid(lm(rk(y[ok]) ~ rk(z[ok])))); n <- sum(ok); c(rho = r, p = 2 * pt(-abs(r * sqrt((n - 3) / (1 - r^2))), n - 3), n = n) }
x <- C[["BAL Th17.1"]][i]; msr <- ms[i, ]
pft <- read.delim(file.path(R0, "pft_series_by_sample.tsv"), stringsAsFactors = FALSE); pft$days_from_BAL <- as.numeric(pft$days_from_BAL)
iv <- sapply(seq_len(nrow(ra)), function(k) { t <- pft[pft$Sample_ID == ra$Sample_ID[k], ]; if (nrow(t) < 2) return(NA); b <- t[which.min(abs(t$days_from_BAL)), ]; f <- t[which.min(abs(t$days_from_BAL - as.numeric(ra$CT2d[k]))), ]; if (b$days_from_BAL == f$days_from_BAL) NA else as.numeric(f$days_from_BAL - b$days_from_BAL) / 365.25 })   # CT2d = days from BAL to the follow-up CT
cov <- list(`ILD score` = num(ra$ILD_Score), `Baseline %FVC` = ra$PFT1_pFVC, `Interval between tests (years)` = iv, `Infection during follow-up` = num(inf), `Methotrexate use` = med_at_BAL(ra$Sample_ID, "MTX"), `Glucocorticoid use` = med_at_BAL(ra$Sample_ID, "GC"), `BAL lymphocyte count` = lym[i])   # no row adjusted for the BAL Th17.1 % (a count is not adjusted for one of its own factors)
pt <- do.call(rbind, lapply(names(cov), function(n) { r <- pc(x, y, cov[[n]]); data.frame(Adjusted_for = n, N = r[["n"]], partial_rho = r[["rho"]], P = r[["p"]]) })); write.csv(pt, file.path(OD, "ST7_partial_correlations.csv"), row.names = FALSE)
V <- list(Age = num(msr$Age), Sex = num(msr$`gender_Male.1._0.`), CDAI = num(msr$CDAI), `DAS28-CRP` = num(msr$DAS28.CRP), CRP = num(msr$CRP), `RA duration (years)` = num(msr$Years_since_diagnosis_of_RA), `ILD score` = num(ra$ILD_Score), `MTX use` = med_at_BAL(ra$Sample_ID, "MTX"), `Glucocorticoid use` = med_at_BAL(ra$Sample_ID, "GC"), `NSAIDs use` = num(msr$NSAIDS),
  `BAL Th17.1` = num(msr$BALF_Th17.1_CCR6pos_CXCR3pos), `BAL MCP-3` = num(msr$BALF_MCP3), `OXPHOS score (GSVA, GO)` = { gv <- read.csv(file.path(R0, "gsva_GO/gsva_GO_scores_n31.csv")); gv$OXPHOS[match(ra$Sample_ID, gv$Sample_ID)] }, `PB Th17` = num(msr$PB_Th17_CCR6pos_CXCR3neg), `Serum IL-27` = num(msr$Serum_IL27), `Serum IL-8` = num(msr$Serum_IL8), `Serum MDC` = num(msr$Serum_MDC), `BAL Th17.1 count` = x, `BAL Th1 count` = C[["BAL Th1"]][i])
# infection: direction-free AUC (>= 0.5, CI mirrored), Wilcoxon for continuous and Fisher's exact test for binary variables
t6 <- do.call(rbind, lapply(names(V), function(n) { v <- V[[n]]; s <- sp(v, y); ok <- is.finite(v); r <- roc(inf[ok], v[ok], levels = c(0, 1), direction = "auto", quiet = TRUE); set.seed(42); ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, quiet = TRUE)
  a <- as.numeric(auc(r)); lo <- ci[1]; hi <- ci[3]; if (a < 0.5) { a <- 1 - a; lh <- c(1 - hi, 1 - lo); lo <- lh[1]; hi <- lh[2] }
  p <- if (all(v[ok] %in% 0:1)) fisher.test(table(factor(v[ok], 0:1), factor(inf[ok], 0:1)))$p.value else wt(v[ok & inf == 1], v[ok & inf == 0])
  data.frame(Variable = n, N_dFVC = s[["n"]], rho_dFVC = s[["rho"]], P_dFVC = s[["p"]], Inf_AUC = a, Inf_lo = lo, Inf_hi = hi, Inf_P = p) })); write.csv(t6, file.path(OD, "ST6_clinical_vs_multiomics.csv"), row.names = FALSE)
stopifnot(abs(t6$rho_dFVC[t6$Variable == "MTX use"] - 0.277) < 5e-4, abs(t6$P_dFVC[t6$Variable == "MTX use"] - 0.212) < 5e-4, abs(t6$Inf_P[t6$Variable == "MTX use"] - 0.0664) < 5e-4)   # medication at BAL
f1 <- read.csv(file.path(R0, "mofa_RA/SourceData_Fig5_n31.csv"))
sd <- data.frame(Sample_ID = ra$Sample_ID, Delta_pFVC = y, Delta_pDLCO = ra$d_pDLCO, BAL_Th171_count_per_mL = x, BAL_Th171_pct = V[["BAL Th17.1"]], BAL_Th1_count_per_mL = V[["BAL Th1 count"]], PB_Th17 = V[["PB Th17"]], MOFA_RA_Factor1 = f1$RA_Factor1[match(ra$Sample_ID, f1$Sample_ID)], Decline_pFVC_5 = ifelse(is.finite(y), as.integer(y <= -5), NA))
write.csv(sd, file.path(OD, "SourceData_Fig4_Fig5_pft.csv"), row.names = FALSE)
print(t9[, c(1, 4, 7, 9, 10)]); print(pt); print(t6[, c(1, 3, 4, 5, 8)])
