# ============================================================================
# 20260605_recolor_Fig5_ROC_CUD.R
# Recolor the three Figure 5 ROC panels (5c disease, 5f infection, 5i CT
# progression) with Color-Universal-Design (CUD) safe colors so the two curves
# in each panel are easy to tell apart for P/D-type colour vision.
#   Tokyo CUD guideline recommended accent set (p.7):
#     single biomarker curve   -> orange  #FF9900  (warm)
#     integrative MOFA factor  -> blue    #0041FF  (cool)
#     5c nested-LOOCV reference -> CUD green #35A16B (bluish green, dashed)
# ONLY the colour values change; all data/logic/AUC/endpoints/titles/positions
# are reproduced verbatim from the source scripts:
#   5c  20260604_fix_Fig5c_ROC_endpoints.R
#   5f  20260601_make_Fig5_infection_ROC_RAfactor.R
#   5i  06_Comprehensive_Figures.R  (Fig5j_RA_progression_ROC block)
# Overwrites the three panel PNG/PDFs in results/panels/.
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(pROC) })

CUD_ORANGE <- "#FF9900"   # single biomarker
CUD_BLUE   <- "#0041FF"   # integrative MOFA factor
CUD_GREEN  <- "#35A16B"   # nested-LOOCV reference (5c)
panel_dir  <- "results/panels"

load("results/RA_ILD_Workspace.RData")
master_all <- master_data
EXCL_v <- if (exists("EXCL")) EXCL else character(0)

tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0,
                                    margin = ggplot2::margin(0,0,3,0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(), plot.margin = unit(c(2,3,2,3),"mm"))
}
add_ep <- function(spec, sens) {
  ord <- order(spec, sens); spec <- spec[ord]; sens <- sens[ord]
  if (spec[1] != 0 || sens[1] != 0) { spec <- c(0, spec); sens <- c(0, sens) }
  if (tail(spec, 1) != 1 || tail(sens, 1) != 1) { spec <- c(spec, 1); sens <- c(sens, 1) }
  list(spec = spec, sens = sens)
}
sv <- function(p, name, w = 80, h = 75) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width=w, height=h, units="mm", dpi=300, device=cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p, width=w, height=h, units="mm", dpi=300, bg="white")
  cat(sprintf("  Saved: %s\n", name))
}

# ---------------------------------------------------------------------------
# 5c — RA vs Sarcoidosis disease ROC  (Factor 1 vs best serum biomarker)
# ---------------------------------------------------------------------------
e6 <- new.env(); load("results/MOFA2_6views_Results.RData", envir = e6)
factors6 <- get("factors", e6)
group_m <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors6)), "Sarcoidosis", "RA")
y_disease <- ifelse(group_m == "RA", 1, 0)
roc_f1 <- roc(y_disease, factors6[, 1], quiet = TRUE, direction = "auto")
serum <- master_all[match(rownames(factors6), master_all$Sample_ID), ]
best <- NULL; ba <- 0
for (cc in grep("^Serum_", colnames(serum), value = TRUE)) {
  v <- as.numeric(serum[[cc]])
  if (sum(!is.na(v)) >= 30) { r <- tryCatch(roc(y_disease, v, quiet=TRUE, direction="auto"), error=function(x) NULL)
    if (!is.null(r) && as.numeric(auc(r)) > ba) { ba <- as.numeric(auc(r)); best <- cc } }
}
roc_best <- roc(y_disease, as.numeric(serum[[best]]), quiet = TRUE, direction = "auto")
cf <- gsub("^Serum_", "", best); cf <- gsub("^GCSF$", "G-CSF", cf); cf <- gsub("^IL(\\d)", "IL-\\1", cf)
ep1 <- add_ep(1 - roc_f1$specificities,  roc_f1$sensitivities)
ep2 <- add_ep(1 - roc_best$specificities, roc_best$sensitivities)
roc_df <- rbind(
  data.frame(Spec = ep1$spec, Sens = ep1$sens, Model = sprintf("Factor 1 (AUC = %.3f)", auc(roc_f1))),
  data.frame(Spec = ep2$spec, Sens = ep2$sens, Model = sprintf("%s (AUC = %.3f)", cf, auc(roc_best))))
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))
cat(sprintf("5c: Factor 1 AUC=%.3f, %s AUC=%.3f\n", auc(roc_f1), cf, auc(roc_best)))
p5c <- ggplot(roc_df, aes(Spec, Sens, color = Model)) +
  geom_step(linewidth = 0.6, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  annotate("segment", x = 0, xend = 1, y = 1, yend = 1, linetype = "dashed", color = CUD_GREEN, linewidth = 0.4) +
  annotate("text", x = 0.98, y = 0.07, size = 2.3, hjust = 1, color = CUD_GREEN,
           label = "Nested LOOCV (AUC = 0.962)", fontface = "italic") +
  scale_color_manual(values = c(CUD_BLUE, CUD_ORANGE)) +    # Factor 1 = blue, serum = orange
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "on") +
  labs(x = "1 – Specificity", y = "Sensitivity", color = "") + ggtitle("RA vs Sarcoidosis") +
  tn(8) + theme(legend.position = c(0.62, 0.22),
                legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm"))
sv(p5c, "Fig5c_diseaseROC_noacronym", 80, 75)

# ---------------------------------------------------------------------------
# 5f / 5i shared: RA-only MOFA factors + EXCL-filtered master + CT
# ---------------------------------------------------------------------------
er <- new.env(); load("results/MOFA2_RA_only_Results.RData", envir = er)
ra_factors <- get("ra_factors", er)
master_data <- master_all[!master_all$Sample_ID %in% EXCL_v, ]

# ---- 5f — infection prediction (BAL Th17.1 vs RA Factor 1) ----
ids <- rownames(ra_factors)
inf <- master_data$respiratory_infection[match(ids, master_data$Sample_ID)]
th171 <- as.numeric(master_data$BALF_Th17.1_CCR6pos_CXCR3pos[match(ids, master_data$Sample_ID)])
fac1 <- as.numeric(ra_factors[, 1])
ok_b <- !is.na(inf) & !is.na(th171); roc_th171 <- roc(inf[ok_b], th171[ok_b], quiet=TRUE, direction="auto")
ok_f <- !is.na(inf) & !is.na(fac1);  roc_fac1  <- roc(inf[ok_f], fac1[ok_f],  quiet=TRUE, direction="auto")
cat(sprintf("5f: BAL Th17.1 AUC=%.3f, RA Factor 1 AUC=%.3f\n", auc(roc_th171), auc(roc_fac1)))
ep1 <- add_ep(1 - roc_th171$specificities, roc_th171$sensitivities)
ep2 <- add_ep(1 - roc_fac1$specificities,  roc_fac1$sensitivities)
roc_5f <- rbind(
  data.frame(Sens = ep1$sens, Spec = ep1$spec, Model = sprintf("BAL Th17.1 (AUC = %.3f)", auc(roc_th171))),
  data.frame(Sens = ep2$sens, Spec = ep2$spec, Model = sprintf("RA Factor 1 (AUC = %.3f)", auc(roc_fac1))))
roc_5f$Model <- factor(roc_5f$Model, levels = unique(roc_5f$Model))
p5f <- ggplot(roc_5f, aes(Spec, Sens, color = Model)) +
  geom_step(linewidth = 0.7, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c(CUD_ORANGE, CUD_BLUE)) +   # Th17.1 = orange, Factor 1 = blue
  labs(x = "1 – Specificity", y = "Sensitivity", color = "") +
  ggtitle("Infection prediction (RA)") +
  tn(8) + theme(legend.position = c(0.62, 0.22),
                legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm"))
sv(p5f, "Fig5f_infection_ROC_RAfactor", 80, 75)

# ---- 5i — CT progression ROC (PB Th17 vs RA Factor 1) ----
load("results/CT_Multiomics_Results.RData")
ct_data <- read.csv("results/tables/master_data_with_CT_all.csv")
ct_unique <- ct_data[!duplicated(ct_data$Sample_ID), ]
ra_ids_mofa <- rownames(ra_factors)
dhl_mofa <- as.numeric(ct_unique$CT_Delta_Healthy_Lung_pct[match(ra_ids_mofa, ct_unique$Sample_ID)])
prog_mofa <- ifelse(dhl_mofa < -10, 1, 0)
pb_th17_vals <- as.numeric(master_data[["PB_Th17_CCR6pos_CXCR3neg"]][match(ra_ids_mofa, master_data$Sample_ID)])
fv_best <- ra_factors[, 1]
ok_t <- !is.na(prog_mofa) & !is.na(pb_th17_vals); roc_th17_5j <- roc(prog_mofa[ok_t], pb_th17_vals[ok_t], quiet=TRUE, direction="auto")
ok_fb <- !is.na(prog_mofa) & !is.na(fv_best);     roc_fac_5j  <- roc(prog_mofa[ok_fb], fv_best[ok_fb],     quiet=TRUE, direction="auto")
cat(sprintf("5i: PB Th17 AUC=%.3f, RA Factor 1 AUC=%.3f\n", auc(roc_th17_5j), auc(roc_fac_5j)))
ep1 <- add_ep(1 - roc_th17_5j$specificities, roc_th17_5j$sensitivities)
ep2 <- add_ep(1 - roc_fac_5j$specificities,  roc_fac_5j$sensitivities)
roc_5j <- rbind(
  data.frame(Sens = ep1$sens, Spec = ep1$spec, Model = sprintf("PB Th17 (AUC = %.3f)", auc(roc_th17_5j))),
  data.frame(Sens = ep2$sens, Spec = ep2$spec, Model = sprintf("RA Factor 1 (AUC = %.3f)", auc(roc_fac_5j))))
roc_5j$Model <- factor(roc_5j$Model, levels = unique(roc_5j$Model))
p5i <- ggplot(roc_5j, aes(Spec, Sens, color = Model)) +
  geom_step(linewidth = 0.7, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c(CUD_ORANGE, CUD_BLUE)) +   # PB Th17 = orange, Factor 1 = blue
  labs(x = "1 – Specificity", y = "Sensitivity", color = "") +
  ggtitle("CT Progression (RA)") +
  tn(8) + theme(legend.position = c(0.65, 0.22),
                legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm"))
sv(p5i, "Fig5j_RA_progression_ROC", 80, 75)

cat("\nDone: 3 ROC panels recolored (CUD orange/blue).\n")
