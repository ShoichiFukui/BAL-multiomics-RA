# 06_common_panel_screens.R -- common panel of 114 baseline markers versus infection, change in %FVC and change in total opacity
#   (Suppl. Data 3, 6, 7; Fig. 3B, 4E). Derived variables (BayesPrism fractions, BAL/PB ratios, GSVA GO gene-set scores) are listed outside the
#   FDR family.
# Produces: common_panel/infection_panel.csv (Supplementary Data 3; Figure 3b), infection_extra.csv (derived rows of Supplementary Data 3),
#   dFVC_panel.csv (Supplementary Data 6; Figure 4e), CT_panel.csv (Supplementary Data 7), panel_definition.csv (114 markers).
# Inputs: planA_data_RA24.csv, master_data_with_CT_all.csv, 17 GSVA scores (gsva_GO/gsva_GO_scores_n31.csv), 01 fractions,
#   results/RA_ILD_Workspace.RData.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages(library(pROC)); set.seed(42)
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "common_panel"); dir.create(OD, recursive = TRUE, showWarnings = FALSE); num <- function(x) suppressWarnings(as.numeric(x))
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); ms <- read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
ms <- check_english_columns(ms)
DROP <- c("BALF_Neutrophil_CD66bpos_CD16pos", "PB_Neutrophil_CD66bpos_CD16pos")
M <- list(); C <- c(); add <- function(n, v, k) { if (n %in% DROP) return(invisible()); v <- num(v); if (sum(is.finite(v)) >= 12 && length(unique(v[is.finite(v)])) > 1) { M[[n]] <<- v; C[n] <<- k } }
for (c0 in grep("^Serum_", names(ms), value = TRUE)) add(c0, ms[[c0]], "Serum cytokine")
for (c0 in grep("^BALF_(?!.*(Th|Treg|CCR6|CD25|CD86|CD66b|Percent))", names(ms), perl = TRUE, value = TRUE)) add(c0, ms[[c0]], "BAL cytokine")
for (c0 in c("BALF_Macrophage_Percent", "BALF_Lymphocyte_Percent", "BALF_Neutrophil_Percent", "BALF_Eosinophil_Percent")) add(c0, ms[[c0]], "BAL differential count")
for (c0 in grep("^BALF_.*(Th|Treg|CD86|CD66b|CD14)", names(ms), value = TRUE)) add(c0, ms[[c0]], "BAL FCM")
for (c0 in c("CD3", "CD4", "CD8", "CD22", "CD16", "CD56", "CD4_CD8_Ratio")) add(paste0("BAL_", c0), ms[[c0]], "BAL lymphocyte subsets (clinical)")
for (c0 in grep("^PB_", names(ms), value = TRUE)) add(c0, ms[[c0]], "PB FCM")
cd4 <- num(ms[["BAL_lymphocytes_per_mL"]]) * num(ms$CD4) / 100   # BAL CD4+ T cells per mL
for (s in c("Th17.1", "Th17", "Th1", "Th2", "Treg")) add(paste0("BAL_", s, "_count"), cd4 * num(ms[[grep(paste0("^BALF_", gsub(".", "\\.", s, fixed = TRUE), "_"), names(ms), value = TRUE)[1]]]) / 100, "BAL CD4+ subset count")
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e); dv <- e$microbiome_results$alpha_diversity
for (x in c("Shannon", "Simpson")) add(paste0("Microbiome_", x), dv[[x]][match(ra$Sample_ID, dv$Sample_ID)], "Microbiome diversity")
stopifnot(length(M) == 114); print(table(C)); write.csv(data.frame(Marker = names(M), Category = C), file.path(OD, "panel_definition.csv"), row.names = FALSE)
X <- list(); XC <- c()
bp <- read.csv(file.path(R0, "bayesprism_n31/celltype_proportions_cellFraction_n31.csv"), row.names = 1); for (ct in colnames(bp)) { X[[paste0("Deconv_", ct)]] <- 100 * bp[ra$Sample_ID, ct]; XC[paste0("Deconv_", ct)] <- "Deconvolution (derived; not in FDR family)" }
X$Ratio_Th17.1_BAL_PB <- num(ms$`BALF_Th17.1_CCR6pos_CXCR3pos`) / (num(ms$`PB_Th17.1_CCR6pos_CXCR3pos`) + 0.01); XC["Ratio_Th17.1_BAL_PB"] <- "BAL/PB ratio (derived; not in FDR family)"
X$Ratio_Treg_BAL_PB <- num(ms$BALF_Treg_CD127neg_CD25pos) / (num(ms$PB_Treg_CD127neg_CD25pos) + 0.01); XC["Ratio_Treg_BAL_PB"] <- "BAL/PB ratio (derived; not in FDR family)"
gv <- read.csv(file.path(R0, "gsva_GO/gsva_GO_scores_n31.csv"))   # 17_GSVA_GO_scores.R
for (g in c("OXPHOS", "Surfactant")) { X[[paste0("GSVA_GO_", g)]] <- gv[[g]][match(ra$Sample_ID, gv$Sample_ID)]; XC[paste0("GSVA_GO_", g)] <- "Gene-set score (GSVA; derived; not in FDR family)" }
inf <- ra$respiratory_infection; stopifnot(sum(inf) == 6, length(inf) == 24); cliff <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))
infrow <- function(n, x, k) { ok <- is.finite(x); a <- x[ok & inf == 1]; b <- x[ok & inf == 0]; r <- roc(inf[ok], x[ok], levels = c(0, 1), direction = "auto", quiet = TRUE); set.seed(42)
  ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, quiet = TRUE); cd <- cliff(a, b)
  au <- as.numeric(auc(r)); lo <- ci[1]; hi <- ci[3]; if (au < 0.5) { au <- 1 - au; lh <- c(1 - hi, 1 - lo); lo <- lh[1]; hi <- lh[2] }   # direction-free AUC (>= 0.5), as in Supplementary Data 3
  data.frame(Variable = n, Category = k, Median_pos = median(a), Median_neg = median(b), N_pos = length(a), N_neg = length(b), AUC = au, CI_lo = lo, CI_hi = hi, Cliff = cd, P = suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value, Direction = ifelse(cd > 0, "Up in Inf+", "Down in Inf+")) }
I <- do.call(rbind, lapply(names(M), function(n) infrow(n, M[[n]], C[[n]]))); I$BH <- p.adjust(I$P, "BH"); I <- I[order(-I$AUC, I$P), ]   # order of Supplementary Data 3
IX <- do.call(rbind, lapply(names(X), function(n) infrow(n, X[[n]], XC[[n]]))); IX$BH <- NA
write.csv(I, file.path(OD, "infection_panel.csv"), row.names = FALSE); write.csv(IX, file.path(OD, "infection_extra.csv"), row.names = FALSE)
sp <- function(x, y) { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 8 || length(unique(x[ok])) < 2) return(c(NA, NA, sum(ok))); r <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = TRUE)); c(r$estimate, r$p.value, sum(ok)) }
scr <- function(y, lab) { s <- do.call(rbind, lapply(names(M), function(n) { r <- sp(M[[n]], y); data.frame(Outcome = lab, Marker = n, Category = C[[n]], N = r[3], rho = r[1], P = r[2]) }))
  t <- is.finite(s$P); s$BH <- NA; s$BH[t] <- p.adjust(s$P[t], "BH"); s[order(s$P), ] }
F <- scr(ra$d_pFVC, "Change in %FVC"); O <- scr(round(-ra$d_aer, 2), "Change in total opacity")   # the Healthy Lung screen is not reported
write.csv(F, file.path(OD, "dFVC_panel.csv"), row.names = FALSE); write.csv(O, file.path(OD, "CT_panel.csv"), row.names = FALSE)
k <- function(D, m, col) { id <- if ("Marker" %in% names(D)) D$Marker else D$Variable; D[id == m, col][1] }

stopifnot(nrow(I) == 114, sum(I$P < 0.05) == 12, abs(min(I$BH) - 0.353) < 5e-4, abs(k(IX, "GSVA_GO_OXPHOS", "P") - 0.00942) < 5e-5)   # values of the article
for (D in list(F, O)) { t <- is.finite(D$P) }

stopifnot(sum(is.finite(F$P)) == 113, abs(min(F$BH, na.rm = TRUE) - 0.200) < 5e-4, sum(is.finite(O$P)) == 114, sum(O$P < 0.05, na.rm = TRUE) == 0)   # values of the article
