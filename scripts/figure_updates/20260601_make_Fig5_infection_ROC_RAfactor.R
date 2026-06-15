# ============================================================================
# 20260601_make_Fig5_infection_ROC_RAfactor.R
# Regenerate the Figure 5 "Infection prediction (RA)" ROC panel so the
# integrated factor curve is the RA-ONLY MOFA Factor 1 instead of the
# full-cohort 6-view Factor 2 (IIV).
#
# Rationale: respiratory infection is an outcome defined only WITHIN RA, so the
# integrated axis should come from the RA-only MOFA (ra_factors[,1]), not from
# the full-cohort model whose Factor 2 (IIV) is optimised across RA+sarcoidosis.
# Kept symmetric to the CT-progression ROC (paper Fig 5j / Fig5j block of
# 06_Comprehensive_Figures.R): a single best biomarker (red) vs the
# RA-only Factor 1 (purple), same theme/legend/size.
#
# Old panel (replaced): BAL Th17.1 (0.870) + IIV (0.750)  [Fig6f_ROC_infection]
# New panel:            BAL Th17.1 (0.870) + RA Factor 1   [Fig5f_infection_ROC_RAfactor]
#
# ROC method = roc(outcome, predictor, direction="auto") verbatim from the
# canonical pipeline (median-based; AUC can be <0.5). No hardcoded AUCs.
# Output: results/panels/Fig5f_infection_ROC_RAfactor.{png,pdf} (80x75mm)
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2); library(pROC) })
EXCL <- character(0)

stopifnot(file.exists("results/RA_ILD_Workspace.RData"))
stopifnot(file.exists("results/MOFA2_RA_only_Results.RData"))
load("results/RA_ILD_Workspace.RData")
e <- new.env(); load("results/MOFA2_RA_only_Results.RData", envir = e)
ra_factors <- get("ra_factors", e)            # 24 x 3 (RA only)
stopifnot(is.matrix(ra_factors), ncol(ra_factors) >= 1)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
panel_dir <- "results/panels"

tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.title = element_text(size = rel(0.95), face = "bold"),
      plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0,
                                margin = ggplot2::margin(0, 0, 3, 0)),
      legend.text = element_text(size = rel(0.8)),
      legend.key.size = unit(3, "mm"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
      axis.line = element_blank(),
      plot.margin = unit(c(2, 3, 2, 3), "mm")
    )
}
sv <- function(p, name, w = 80, h = 75) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p,
         width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s (%dx%d mm)\n", name, w, h))
}
add_ep <- function(spec, sens) {
  ord <- order(spec, sens); spec <- spec[ord]; sens <- sens[ord]
  if (spec[1] != 0 || sens[1] != 0) { spec <- c(0, spec); sens <- c(0, sens) }
  if (tail(spec, 1) != 1 || tail(sens, 1) != 1) { spec <- c(spec, 1); sens <- c(sens, 1) }
  list(spec = spec, sens = sens)
}

# --- align outcomes/predictors to the RA-only MOFA sample order ---
ids <- rownames(ra_factors)
inf <- master_data$respiratory_infection[match(ids, master_data$Sample_ID)]
th171 <- as.numeric(master_data$BALF_Th17.1_CCR6pos_CXCR3pos[
  match(ids, master_data$Sample_ID)])
fac1 <- as.numeric(ra_factors[, 1])

# single best biomarker: BAL Th17.1
ok_b <- !is.na(inf) & !is.na(th171)
roc_th171 <- roc(inf[ok_b], th171[ok_b], quiet = TRUE, direction = "auto")
# RA-only integrated factor: Factor 1
ok_f <- !is.na(inf) & !is.na(fac1)
roc_fac1 <- roc(inf[ok_f], fac1[ok_f], quiet = TRUE, direction = "auto")

cat(sprintf("RA n=%d, infection events=%d\n", sum(!is.na(inf)), sum(inf == 1, na.rm = TRUE)))
cat(sprintf("  BAL Th17.1   AUC = %.3f (n=%d)\n", auc(roc_th171), sum(ok_b)))
cat(sprintf("  RA Factor 1  AUC = %.3f (n=%d)\n", auc(roc_fac1), sum(ok_f)))

ep1 <- add_ep(1 - roc_th171$specificities, roc_th171$sensitivities)
ep2 <- add_ep(1 - roc_fac1$specificities,  roc_fac1$sensitivities)
roc_df <- rbind(
  data.frame(Sens = ep1$sens, Spec = ep1$spec,
             Model = sprintf("BAL Th17.1 (AUC = %.3f)", auc(roc_th171))),
  data.frame(Sens = ep2$sens, Spec = ep2$spec,
             Model = sprintf("RA Factor 1 (AUC = %.3f)", auc(roc_fac1)))
)
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))

sv(ggplot(roc_df, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.7, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                 color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#C0392B", "#9B59B6")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     ggtitle("Infection prediction (RA)") +
     tn(8) + theme(legend.position = c(0.62, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig5f_infection_ROC_RAfactor", 80, 75)
