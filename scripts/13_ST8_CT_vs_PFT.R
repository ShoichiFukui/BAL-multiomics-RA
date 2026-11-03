# 13_ST8_CT_vs_PFT.R -- Supplementary Table 8 (automated CT measures and pulmonary function, RA): Spearman correlations of total lung
#   volume and total opacity (GMM GGO + Consolidation + Dense Tissue = 100 - Air - Healthy Lung) with FVC, %FVC and %DLCO at baseline
#   (test nearest to BAL) and at follow-up (test nearest to the follow-up CT), and of the changes between the two time points,
#   including the change in GGO, Consolidation and Dense Tissue (Figure 4i; Results).
# Produces: <ART_TABLES>/n31/ST8/ST8_CT_vs_PFT.csv.
# Inputs: planA_data_RA24.csv (controlled access).
# Changes are computed at the recorded precision (CT percentages and volume to 2 decimals, %FVC and %DLCO to 1 decimal) so that real ties
#   are kept; Spearman exact (normal approximation with ties). BH adjustment as stated in the table legend: within each time point for
#   total lung volume and total opacity, and across the nine change rows of GGO, Consolidation and Dense Tissue as a separate family.
# The BH step is written out here from the legend definition. Rows not reported (Healthy Lung, Air) are not computed.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "ST8"); dir.create(OD, recursive = TRUE, showWarnings = FALSE)
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); stopifnot(nrow(ra) == 24)
sp <- function(o, y) { ok <- is.finite(o) & is.finite(y); r <- suppressWarnings(cor.test(o[ok], y[ok], method = "spearman", exact = TRUE)); c(n = sum(ok), rho = unname(r$estimate), p = r$p.value) }
row <- function(t, k, j, o, y, fam) { r <- sp(o, y); data.frame(Time = t, CT_measure = k, Pulmonary_function = j, N = r[["n"]], rho = r[["rho"]], P = r[["p"]], Family = fam) }
# ---- baseline and follow-up (total opacity = 100 - aerated lung) ----
lv <- list(Baseline = list(V = ra$CT_T1_total_lung_cm3, O = 100 - ra$T1_aer, P = c(FVC = "PFT1_FVC", `%FVC` = "PFT1_pFVC", `%DLCO` = "PFT1_pDLCO")),
           `Follow-up` = list(V = ra$CT_T2_total_lung_cm3, O = 100 - ra$T2_aer, P = c(FVC = "PFT2_FVC", `%FVC` = "PFT2_pFVC", `%DLCO` = "PFT2_pDLCO")))
out <- do.call(rbind, lapply(names(lv), function(t) { x <- lv[[t]]; rbind(
  do.call(rbind, lapply(names(x$P), function(j) row(t, "Total lung volume (cm3)", j, x$V, ra[[x$P[[j]]]], t))),
  do.call(rbind, lapply(names(x$P), function(j) row(t, "Total opacity (%)", j, x$O, ra[[x$P[[j]]]], t)))) }))
# ---- change rows ----
CT <- list(`Total lung volume (cm3)` = round(ra$CT_T2_total_lung_cm3 - ra$CT_T1_total_lung_cm3, 2), `Total opacity (%)` = round(-ra$d_aer, 2),
           `GGO (%)` = round(ra$CT_Delta_GGO_pct, 2), `Consolidation (%)` = round(ra$CT_Delta_Consolidation_pct, 2), `Dense Tissue (%)` = round(ra$CT_Delta_Dense_Tissue_pct, 2))
PF <- list(`Relative change in FVC (%)` = ra$dFVC_rel, `Change in %FVC` = round(ra$d_pFVC, 1), `Change in %DLCO` = round(ra$d_pDLCO, 1))
ch <- do.call(rbind, lapply(names(CT), function(k) do.call(rbind, lapply(names(PF), function(j)
  row("Change", k, j, CT[[k]], PF[[j]], if (k %in% c("Total lung volume (cm3)", "Total opacity (%)")) "Change" else "Change (abnormal classes)")))))
out <- rbind(out, ch)
out$BH <- ave(out$P, out$Family, FUN = function(p) p.adjust(p, "BH"))
stopifnot(table(out$Family)[["Change (abnormal classes)"]] == 9, all(table(out$Family)[c("Baseline", "Follow-up", "Change")] == 6))
write.csv(out, file.path(OD, "ST8_CT_vs_PFT.csv"), row.names = FALSE)
x <- out; x$rho <- round(x$rho, 2); x$P <- signif(x$P, 3); x$BH <- signif(x$BH, 3); print(x, row.names = FALSE)
