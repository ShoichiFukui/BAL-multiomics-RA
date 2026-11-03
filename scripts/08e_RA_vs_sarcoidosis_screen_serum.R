# 08e_RA_vs_sarcoidosis_screen_serum.R -- RA-versus-sarcoidosis screen of all variables (input of 08c; Fig. 2e) and serum cytokines versus the healthy reference (Suppl. Table 1; Suppl. Fig. 2).
# Produces: SD2_effect_sizes_RA_vs_Sarc.csv, ST1_serum_vs_healthy.csv.
# Inputs: results/RA_ILD_Workspace.RData, bayesprism_n31/celltype_proportions_cellFraction_n31.csv (01b_BayesPrism_n31.R),
#   <ART_DIR>/41plex_healthy_control.xlsx.
# Output file ST1_serum_vs_healthy.csv holds the values of Supplementary Table 1 (all 41 serum cytokines: RA, sarcoidosis, healthy
#   reference).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(effsize); library(readxl) })
cfg <- set_cohort("n31")   # n = 31 only (reads the n = 31 BayesPrism run)

# ---- data -------------------------------------------------------------------
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
md <- e$master_data; rm(e); invisible(gc())
md <- md[!md$Sample_ID %in% cfg$EXCL, ]
check_n(md$Sample_ID, cfg, "master_data")

bp_cell <- read.csv(file.path(ART_TABLES, "n31", "bayesprism_n31", "celltype_proportions_cellFraction_n31.csv"), row.names = 1)   # n = 31 BayesPrism re-run (01b_BayesPrism_n31.R)
bp_cell <- bp_cell[!rownames(bp_cell) %in% cfg$EXCL, ]
check_n(rownames(bp_cell), cfg, "bp_cell")

# ---- (1) screen ---------------------------------------------------------------
compute_es <- function(v, g) {
  v1 <- as.numeric(v[g == "RA"]); v2 <- as.numeric(v[g == "Control"])
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  if (length(v1) < 3 || length(v2) < 3)
    return(data.frame(n1 = length(v1), n2 = length(v2), median1 = NA, median2 = NA,
                      p_value = NA, cliff_delta = NA, magnitude = NA))
  wt <- suppressWarnings(wilcox.test(v1, v2, exact = TRUE))
  cd <- cliff.delta(v1, v2)
  data.frame(n1 = length(v1), n2 = length(v2), median1 = median(v1), median2 = median(v2),
             p_value = wt$p.value, cliff_delta = unname(cd$estimate),
             magnitude = as.character(cd$magnitude))
}

cyto <- grep("^BALF_|^Serum_", names(md), value = TRUE)
cyto <- cyto[!grepl("Treg|Th|Macro|Neutro|Lympho|CD", cyto)]
balf_fcm <- grep("^BALF_.*(Treg|Th|Macro|Neutro)", names(md), value = TRUE)
pb_fcm   <- grep("^PB_", names(md), value = TRUE)

rows <- list()
for (v in cyto) {
  x <- suppressWarnings(as.numeric(md[[v]]))
  if (all(is.na(x))) next
  rows[[v]] <- cbind(Variable = v, compute_es(x, md$Sample_Group),
                     Category = ifelse(grepl("^BALF_", v), "BALF_Cytokine", "Serum_Cytokine"))
}
for (v in c(balf_fcm, pb_fcm)) {
  rows[[v]] <- cbind(Variable = v, compute_es(as.numeric(md[[v]]), md$Sample_Group),
                     Category = ifelse(grepl("^BALF_", v), "BALF_FCM", "PB_FCM"))
}
g_bp <- ifelse(is_sarc(rownames(bp_cell)), "Control", "RA")
for (ct in colnames(bp_cell)) {
  rows[[paste0("bp_", ct)]] <- cbind(Variable = ct, compute_es(bp_cell[[ct]] * 100, g_bp),
                                     Category = "Deconvolution")
}
es <- do.call(rbind, rows); rownames(es) <- NULL
es$Category[grepl("^BALF_.*_Percent$", es$Variable)] <- "BAL_differential_count"   # BAL differential counts of the clinical laboratory (Supplementary Data 2, Figure 2e)
es$p_adjusted <- p.adjust(es$p_value, method = "BH")
es <- es[order(es$p_value), ]
p_plasma <- es$p_value[es$Category == "Deconvolution" & es$Variable == "Plasma"]
stopifnot(abs(p_plasma - 0.06095) < 1e-4)   # Figure 2f / Supplementary Data 2: plasma cell p = 0.061 on the n = 31 BayesPrism re-run
write.csv(es, file.path(cfg$tables, "SD2_effect_sizes_RA_vs_Sarc.csv"), row.names = FALSE)

sig <- es[!is.na(es$p_adjusted) & es$p_adjusted < 0.05, ]

print(data.frame(Variable = sig$Variable, Category = sig$Category,
                 delta = round(sig$cliff_delta, 3), p = signif(sig$p_value, 3),
                 BH = signif(sig$p_adjusted, 3)), row.names = FALSE)

# ---- (2) top-10 serum cytokines vs healthy reference -------------------------
hc <- read_excel(file.path(ART_DIR, "41plex_healthy_control.xlsx"),
                 sheet = "FIXED")
cyto_map <- c(GCSF = "G-CSF", IL6 = "IL-6", IL13 = "IL-13", IL17A = "IL-17A", IFNa2 = "IFNa2",
  IL22 = "IL-22", IL4 = "IL-4", IL2 = "IL-2", GROa = "GROa", FGF2 = "FGF-2", TNFb = "TNFb",
  MCP3 = "MCP-3", IL10 = "IL-10", IL12p70 = "IL-12p70", IL17F = "IL-17F", IL1RA = "IL-1RA",
  TGFa = "TGFa", IL5 = "IL-5", IL1b = "IL-1b", IL7 = "IL-7", IL1a = "IL-1a",
  Fractalkine = "Fractalkine", IL15 = "IL-15", IFNg = "IFNg", MIP1a = "MIP-1a", IL8 = "IL-8",
  VEGFA = "VEGF-A", FLT3L = "FLT-3L", IL12p40 = "IL-12p40", EGF = "EGF", Eotaxin = "Eotaxin",
  GMCSF = "GM-CSF", MDC = "MDC", MIP1b = "MIP-1b", PDGFAA = "PDGF-AA", TNFa = "TNFa",
  sCD40L = "sCD40L", IP10 = "IP-10", MCP1 = "MCP-1", IL18 = "IL-18", IL27 = "IL-27")

ser <- es[es$Category == "Serum_Cytokine" & !is.na(es$cliff_delta), ]
ser <- ser[order(-abs(ser$cliff_delta), ser$p_value), ]
top10 <- head(ser, 10)
top10_p <- head(ser[order(ser$p_value), ], 10)$Variable

# Supplementary Table 1: all 41 serum cytokines of the panel (sCD40L included; it is not in the screen above)
hres <- do.call(rbind, lapply(paste0("Serum_", names(cyto_map)), function(v) {
  s <- sub("^Serum_", "", v); hv <- as.numeric(hc[[cyto_map[[s]]]]); hv <- hv[!is.na(hv)]
  x <- as.numeric(md[[v]]); g <- md$Sample_Group
  ra <- x[g == "RA" & !is.na(x)]; sa <- x[g == "Control" & !is.na(x)]
  q <- unname(quantile(hv, c(0.025, 0.975)))
  data.frame(Cytokine = s, Healthy_n = length(hv), Healthy_median = median(hv),
             Healthy_2.5 = q[1], Healthy_97.5 = q[2],
             RA_n = length(ra), RA_median = median(ra),
             RA_Q1 = unname(quantile(ra, .25)), RA_Q3 = unname(quantile(ra, .75)),
             Sarc_n = length(sa), Sarc_median = median(sa),
             Sarc_Q1 = unname(quantile(sa, .25)), Sarc_Q3 = unname(quantile(sa, .75)),
             RA_median_above_upper95 = median(ra) > q[2],
             Sarc_median_above_upper95 = median(sa) > q[2],
             RA_vs_Sarc_p = suppressWarnings(wilcox.test(ra, sa, exact = TRUE))$p.value,
             cliff_delta = mean(sign(outer(ra, sa, "-"))),
             in_top10_this_cohort = v %in% top10$Variable)
}))
hres$RA_vs_Sarc_BH <- p.adjust(hres$RA_vs_Sarc_p, "BH")   # across the 41 serum cytokines (Supplementary Table 1)
stopifnot(nrow(hres) == 41)
write.csv(hres, file.path(cfg$tables, "ST1_serum_vs_healthy.csv"), row.names = FALSE)
