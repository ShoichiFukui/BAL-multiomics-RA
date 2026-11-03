# 15_Th_subsets_BAL_PB.R -- Figure 2h (BAL CD4+ T cell subsets by FCM, RA vs sarcoidosis, exact Wilcoxon) and Figure 2i/2j (BAL vs PB
#   Th17.1 and Treg proportions, Spearman; all 31 patients, and within RA and within sarcoidosis as quoted in the Results).
# Produces: <ART_TABLES>/n31/Fig2h_BAL_Th_subsets.csv, Fig2ij_BAL_vs_PB.csv (value-check inputs of figures/panels_screen_MOFA.R) and
#   Fig2ij_BAL_vs_PB_by_group.csv.
# Inputs: results/RA_ILD_Workspace.RData (master_data; controlled access).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
cfg <- set_cohort("n31")

e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
md <- e$master_data; rm(e); invisible(gc())
md <- md[!md$Sample_ID %in% cfg$EXCL, ]
check_n(md$Sample_ID, cfg, "master_data")
g <- ifelse(md$Sample_Group == "Control", "Sarcoidosis", "RA")

# ---- Fig 2h (corrected order: tc[i] <-> label[i]) ---------------------------
tc <- c(Th1 = "BALF_Th1_CCR6neg_CXCR3pos", `Th2-enriched` = "BALF_Th2_CCR6neg_CXCR3neg",
        Th17 = "BALF_Th17_CCR6pos_CXCR3neg", Th17.1 = "BALF_Th17.1_CCR6pos_CXCR3pos",
        Treg = "BALF_Treg_CD127neg_CD25pos")
f2i <- do.call(rbind, lapply(names(tc), function(lab) {
  v <- as.numeric(md[[tc[[lab]]]]); ok <- !is.na(v)
  a <- v[ok & g == "RA"]; b <- v[ok & g == "Sarcoidosis"]
  wt <- suppressWarnings(wilcox.test(a, b, exact = TRUE))
  data.frame(Subset = lab, Variable = tc[[lab]], n_RA = length(a), n_Sarc = length(b),
             median_RA = median(a), IQR_RA = sprintf("%.1f-%.1f", quantile(a, .25), quantile(a, .75)),
             median_Sarc = median(b), IQR_Sarc = sprintf("%.1f-%.1f", quantile(b, .25), quantile(b, .75)),
             p = wt$p.value, ties = any(duplicated(c(a, b))))
}))
f2i$BH_5subsets <- p.adjust(f2i$p, "BH")
write.csv(f2i, file.path(cfg$tables, "Fig2h_BAL_Th_subsets.csv"), row.names = FALSE)

# ---- Fig 2i/2j ---------------------------------------------------------------
pairs <- list(Th17.1 = c("BALF_Th17.1_CCR6pos_CXCR3pos", "PB_Th17.1_CCR6pos_CXCR3pos"),
              Treg   = c("BALF_Treg_CD127neg_CD25pos", "PB_Treg_CD127neg_CD25pos"))
jk <- do.call(rbind, lapply(names(pairs), function(nm) {
  ba <- as.numeric(md[[pairs[[nm]][1]]]); pb <- as.numeric(md[[pairs[[nm]][2]]])
  ok <- !is.na(ba) & !is.na(pb)
  sp <- suppressWarnings(cor.test(ba[ok], pb[ok], method = "spearman", exact = TRUE))
  data.frame(Subset = nm, N = sum(ok), N_RA = sum(ok & g == "RA"), N_Sarc = sum(ok & g == "Sarcoidosis"),
             Spearman_rho = unname(sp$estimate), p = sp$p.value)
}))
write.csv(jk, file.path(cfg$tables, "Fig2ij_BAL_vs_PB.csv"), row.names = FALSE)

# ---- within-group correlations (Results: within RA) ----
g <- ifelse(md$Sample_Group == "Control", "Sarcoidosis", "RA")
out <- data.frame()
for (s in c("Th17.1_CCR6pos_CXCR3pos", "Treg_CD127neg_CD25pos")) for (grp in c("All", "RA", "Sarcoidosis")) {
  b <- as.numeric(md[[paste0("BALF_", s)]]); p <- as.numeric(md[[paste0("PB_", s)]]); k <- !is.na(b) & !is.na(p) & (grp == "All" | g == grp)
  r <- suppressWarnings(cor.test(b[k], p[k], method = "spearman", exact = TRUE))
  out <- rbind(out, data.frame(Subset = sub("_.*", "", s), Group = grp, n = sum(k), rho = round(unname(r$estimate), 3), p = signif(r$p.value, 3)))
}
print(out, row.names = FALSE)
write.csv(out, file.path(cfg$tables, "Fig2ij_BAL_vs_PB_by_group.csv"), row.names = FALSE)
