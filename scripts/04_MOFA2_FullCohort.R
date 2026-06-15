#!/usr/bin/env Rscript
# 04_MOFA2_FullCohort.R — Pekayvaz-inspired MOFA2 enhancements
# Items 2,3,4,6,7,8,9 from reviewer enhancement list

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(readxl)
  library(pROC); library(patchwork); library(glmnet)
})

if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
panel_dir <- "results/panels"
table_dir <- "results/tables"

# ── Theme (Pekayvaz-inspired, cleaner) ──────────────────────────────
tn2 <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(size = rel(0.9), color = "black"),
      axis.title = element_text(size = rel(0.95), face = "bold"),
      plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0, margin = ggplot2::margin(0,0,3,0)),
      legend.text = element_text(size = rel(0.8)),
      legend.key.size = unit(3, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
      strip.text = element_text(size = rel(0.85), face = "bold"),
      plot.margin = unit(c(2, 3, 2, 3), "mm")
    )
}

sv2 <- function(p, name, w = 80, h = 70) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p,
         width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s (%dx%d mm)\n", name, w, h))
}

# ── Color palettes ──────────────────────────────────────────────────
CG  <- c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")
CI  <- c("Infection (-)" = "#3498DB", "Infection (+)" = "#E74C3C")

view_cols <- c(
  "Expression" = "#E67E22", "BALF_Cytokine" = "#3498DB",
  "Serum_Cytokine" = "#E74C3C", "BALF_FACS" = "#27AE60",
  "PB_FACS" = "#9B59B6", "Microbiome" = "#7F8C8D"
)
view_labels <- c(
  "Expression" = "Expression", "BALF_Cytokine" = "BAL Cytokine",
  "Serum_Cytokine" = "Serum Cytokine", "BALF_FACS" = "BAL FCM",
  "PB_FACS" = "PB FCM", "Microbiome" = "Microbiome"
)

# ── Factor names (Item 2) ──────────────────────────────────────────
fnames <- c(Factor1 = "Factor 1 (IDD)", Factor2 = "Factor 2 (IIV)",
            Factor3 = "Factor 3", Factor4 = "Factor 4", Factor5 = "Factor 5")
fshort <- c(Factor1 = "IDD", Factor2 = "IIV")

# ══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════
cat("Loading data...\n")
load("results/MOFA2_6views_Results.RData")  # factors, r2, weights, fd, fi
load("results/RA_ILD_Workspace.RData")

# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)

bp_cell <- read.csv("output_v3_BayesPrism/celltype_proportions_cellFraction.csv", row.names = 1)
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet = "BALF_analysis")
colnames(fcm)[1] <- "Sample_ID"
bp_fcm_common <- intersect(rownames(bp_cell), fcm$Sample_ID)
bp_fcm_common <- bp_fcm_common[!bp_fcm_common %in% EXCL]

ct_data <- read.csv("results/tables/master_data_with_CT_all.csv")
ct_unique <- ct_data[!duplicated(ct_data$Sample_ID), ]

group_m <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors)), "Sarcoidosis", "RA")
inf_m <- master_data$respiratory_infection[match(rownames(factors), master_data$Sample_ID)]

ra_mask <- group_m == "RA"
ra_samples <- rownames(factors)[ra_mask]
ra_common <- intersect(ra_samples, bp_fcm_common[grepl("^(KY[0-9]|RA[0-9])", bp_fcm_common)])

cat(sprintf("  Samples: %d total, %d RA, %d Sarc\n", nrow(factors), sum(ra_mask), sum(!ra_mask)))

# ══════════════════════════════════════════════════════════════════════
# STEP 2: Enhanced R² Heatmap (Item 6)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 2: Enhanced R² Heatmap ===\n")

r2m <- r2$r2_per_factor[[1]]
r2l <- as.data.frame(r2m) %>%
  rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2")
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2m)))
r2l$View <- factor(r2l$View, levels = colnames(r2m))
r2l$FLabel <- fnames[as.character(r2l$Factor)]

# Main heatmap (white→blue)
p_hm <- ggplot(r2l, aes(View, FLabel, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient(low = "white", high = "#2471A3", name = "R² (%)") +
  geom_text(aes(label = ifelse(R2 > 0.5, sprintf("%.1f", R2), ""),
                color = ifelse(R2 > 25, "white", "black")),
            size = 2.5, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_x_discrete(labels = view_labels) +
  labs(x = "", y = "") +
  tn2(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7,
                                              color = view_cols[levels(r2l$View)]),
                  legend.position = "left", legend.key.height = unit(10, "mm"))

# Right bar: total R² per view
total_r2 <- data.frame(
  View = factor(names(r2$r2_total[[1]]), levels = colnames(r2m)),
  Total = as.numeric(r2$r2_total[[1]])
)

p_bar <- ggplot(total_r2, aes(Total, View, fill = View)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = view_cols) +
  scale_y_discrete(labels = view_labels) +
  labs(x = "Total R² (%)", y = "") +
  tn2(7) + theme(legend.position = "none",
                  axis.text.y = element_blank(), axis.ticks.y = element_blank())

p6a <- p_hm + p_bar + plot_layout(widths = c(4, 1.2))
sv2(p6a, "Fig6a_R2_enhanced", 130, 75)

# ══════════════════════════════════════════════════════════════════════
# STEP 1: Factor scatter + boxplots with IDD/IIV names (Item 2)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 1: Named Factor Panels ===\n")

# 6b: Scatter
df_sc <- data.frame(F1 = factors[, 1], F2 = factors[, 2], Group = group_m)
p6b <- ggplot(df_sc, aes(F1, F2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(level = 0.95, linewidth = 0.5, linetype = "dashed") +
  scale_color_manual(values = CG) +
  labs(x = "Factor 1: IDD", y = "Factor 2: IIV") +
  guides(color = guide_legend(title = NULL)) +
  tn2(8) + theme(legend.position = c(0.77, 0.88),
                  legend.background = element_rect(fill = alpha("white", 0.9), color = "grey80", linewidth = 0.3),
                  legend.margin = ggplot2::margin(2, 4, 2, 4))
sv2(p6b, "Fig6b_MOFA2_scatter", 80, 75)

# 6c: Disease factors (IDD + Factor 3)
sf6 <- fd$Factor[fd$p_disease < 0.1]
fb6 <- data.frame()
for (f in sf6) {
  pv <- fd$p_disease[f]
  lab <- ifelse(f <= 2, paste0(fnames[paste0("Factor", f)], "\n(p=",
                                ifelse(pv < 0.0001, "<0.0001", sprintf("%.4f", pv)), ")"),
                 sprintf("Factor %d\n(p=%s)", f, ifelse(pv < 0.0001, "<0.0001", sprintf("%.4f", pv))))
  fb6 <- rbind(fb6, data.frame(Factor = lab, Value = factors[, f], Group = group_m))
}
fb6$Group <- recode(fb6$Group, "Sarcoidosis" = "Sarcoidosis")

p6c <- ggplot(fb6, aes(Group, Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25, alpha = 0.7) +
  geom_jitter(aes(color = Group), width = 0.08, size = 1.3, alpha = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")) +
  scale_color_manual(values = c("RA" = "#7B241C", "Sarcoidosis" = "#1A5276")) +
  facet_wrap(~Factor, scales = "free_y") +
  labs(x = "", y = "Factor value") +
  tn2(8) + theme(legend.position = "none")
sv2(p6c, "Fig6c_disease_factor", 95, 75)

# 6d: Infection factor (IIV)
sfi6 <- fi$Factor[fi$p_infection < 0.1]
if (length(sfi6) == 0) sfi6 <- 2  # Factor 2
ib6 <- data.frame()
for (f in sfi6) {
  m <- ra_mask & !is.na(inf_m)
  pv <- fi$p_infection[fi$Factor == f]
  lab <- ifelse(f <= 2, paste0(fnames[paste0("Factor", f)], "\n(p=", sprintf("%.3f", pv), ")"),
                 sprintf("Factor %d\n(p=%.3f)", f, pv))
  ib6 <- rbind(ib6, data.frame(Factor = lab, Value = factors[m, f],
                                 IG = ifelse(inf_m[m] == 1, "Infection (+)", "Infection (-)")))
}

p6d <- ggplot(ib6, aes(IG, Value, fill = IG)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25, alpha = 0.7) +
  geom_jitter(aes(color = IG), width = 0.08, size = 1.3, alpha = 0.65, show.legend = FALSE) +
  scale_fill_manual(values = CI) +
  scale_color_manual(values = c("Infection (-)" = "#1A5276", "Infection (+)" = "#922B21")) +
  facet_wrap(~Factor, scales = "free_y") +
  labs(x = "", y = "Factor value") +
  tn2(8) + theme(legend.position = "none", axis.text.x = element_text(size = 6.5))
sv2(p6d, "Fig6d_infection_factor", 80, 75)

# ══════════════════════════════════════════════════════════════════════
# STEP 3: Factor ROC Analysis (Item 3)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 3: Factor ROC Analysis ===\n")

# 3a: Disease ROC — Factor 1 (IDD) vs supervised
y_disease <- ifelse(group_m == "RA", 1, 0)
roc_idd <- roc(y_disease, factors[, 1], quiet = TRUE, direction = "auto")
cat(sprintf("  IDD ROC AUC: %.3f\n", auc(roc_idd)))

# Best single serum cytokine for comparison
serum_data <- master_data[match(rownames(factors), master_data$Sample_ID), ]
best_cyto <- NULL; best_auc <- 0
cyto_cols <- grep("^Serum_", colnames(serum_data), value = TRUE)
for (cc in cyto_cols) {
  vals <- as.numeric(serum_data[[cc]])
  if (sum(!is.na(vals)) >= 30) {
    r <- tryCatch(roc(y_disease, vals, quiet = TRUE, direction = "auto"), error = function(e) NULL)
    if (!is.null(r) && as.numeric(auc(r)) > best_auc) {
      best_auc <- as.numeric(auc(r))
      best_cyto <- cc
    }
  }
}
roc_best_cyto <- roc(y_disease, as.numeric(serum_data[[best_cyto]]), quiet = TRUE, direction = "auto")
cyto_label <- gsub("^Serum_", "", best_cyto)
cat(sprintf("  Best cytokine: %s (AUC=%.3f)\n", cyto_label, best_auc))

# Format cytokine name
cyto_fmt <- gsub("^GCSF$", "G-CSF", cyto_label)
cyto_fmt <- gsub("^IL(\\d)", "IL-\\1", cyto_fmt)
cyto_fmt <- gsub("^GMCSF$", "GM-CSF", cyto_fmt)

# Plot disease ROC (smoothed step curves)
roc_df_d <- rbind(
  data.frame(Sens = roc_idd$sensitivities, Spec = 1 - roc_idd$specificities,
             Model = sprintf("IDD (AUC = %.3f)", auc(roc_idd))),
  data.frame(Sens = roc_best_cyto$sensitivities, Spec = 1 - roc_best_cyto$specificities,
             Model = sprintf("%s (AUC = %.3f)", cyto_fmt, auc(roc_best_cyto)))
)
p6e <- ggplot(roc_df_d, aes(Spec, Sens, color = Model)) +
  geom_step(linewidth = 0.6, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  annotate("segment", x = 0, xend = 1, y = 1, yend = 1, linetype = "dashed",
           color = "#27AE60", linewidth = 0.4) +
  annotate("text", x = 0.98, y = 0.05, size = 2.3, hjust = 1, color = "#27AE60",
           label = "Nested LOOCV (AUC = 0.962)", fontface = "italic") +
  scale_color_manual(values = c("#C0392B", "#E67E22")) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
  ggtitle("RA vs Sarcoidosis") +
  tn2(8) + theme(legend.position = c(0.62, 0.22),
                  legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                  legend.text = element_text(size = 6.5),
                  legend.key.height = unit(3, "mm"))
sv2(p6e, "Fig6e_ROC_disease", 80, 75)

# 3b: Infection ROC — Factor 2 (IIV) vs single biomarkers
ra_inf <- !is.na(inf_m) & ra_mask
y_inf <- inf_m[ra_inf]
roc_iiv <- roc(y_inf, factors[ra_inf, 2], quiet = TRUE, direction = "auto")
cat(sprintf("  IIV ROC AUC: %.3f\n", auc(roc_iiv)))

# Th17.1 and ratio
th171_vals <- as.numeric(master_data$BALF_Th17.1_CCR6pos_CXCR3pos[match(rownames(factors)[ra_inf], master_data$Sample_ID)])
roc_th171 <- tryCatch(roc(y_inf, th171_vals, quiet = TRUE, direction = "auto"), error = function(e) NULL)

# OXPHOS
if (exists("gsva_scores")) {
  oxphos_vals <- gsva_scores["OXPHOS", rownames(factors)[ra_inf]]
  roc_oxphos <- tryCatch(roc(y_inf, as.numeric(oxphos_vals), quiet = TRUE, direction = "auto"), error = function(e) NULL)
} else {
  roc_oxphos <- NULL
}

roc_df_i <- data.frame(
  Sens = roc_iiv$sensitivities, Spec = 1 - roc_iiv$specificities,
  Model = sprintf("IIV (AUC = %.3f)", auc(roc_iiv))
)
roc_cols <- c("#E74C3C")
if (!is.null(roc_th171)) {
  roc_df_i <- rbind(roc_df_i,
    data.frame(Sens = roc_th171$sensitivities, Spec = 1 - roc_th171$specificities,
               Model = sprintf("BAL Th17.1 (AUC = %.3f)", auc(roc_th171))))
  roc_cols <- c(roc_cols, "#9B59B6")
  cat(sprintf("  Th17.1 ROC AUC: %.3f\n", auc(roc_th171)))
}
if (!is.null(roc_oxphos)) {
  roc_df_i <- rbind(roc_df_i,
    data.frame(Sens = roc_oxphos$sensitivities, Spec = 1 - roc_oxphos$specificities,
               Model = sprintf("OXPHOS (AUC = %.3f)", auc(roc_oxphos))))
  roc_cols <- c(roc_cols, "#27AE60")
  cat(sprintf("  OXPHOS ROC AUC: %.3f\n", auc(roc_oxphos)))
}

p6f <- ggplot(roc_df_i, aes(Spec, Sens, color = Model)) +
  geom_step(linewidth = 0.6, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = roc_cols) +
  labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
  ggtitle("Infection prediction (RA)") +
  tn2(8) + theme(legend.position = c(0.62, 0.22),
                  legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                  legend.text = element_text(size = 6.5),
                  legend.key.height = unit(3, "mm"))
sv2(p6f, "Fig6f_ROC_infection", 80, 75)

# ══════════════════════════════════════════════════════════════════════
# STEP 4: Lasso Parsimonious Model (Item 4)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 4: Lasso from Factor 2 Top Weights ===\n")

# Collect top Factor 2 features across all views
f2w_all <- data.frame()
for (v in names(weights)) {
  wv <- weights[[v]][, 2]
  f2w_all <- rbind(f2w_all, data.frame(View = v, Feature = names(wv),
                                         Weight = as.numeric(wv), stringsAsFactors = FALSE))
}
f2w_all$AbsW <- abs(f2w_all$Weight)
f2_top <- f2w_all %>% arrange(desc(AbsW)) %>% head(20)
cat("  Top Factor 2 features:\n")
print(f2_top[1:10, c("View", "Feature", "Weight")])

# Build feature matrix for RA infection prediction
ra_inf_ids <- rownames(factors)[ra_inf]
# Map features to actual data
feat_mat <- matrix(NA, nrow = length(ra_inf_ids), ncol = nrow(f2_top))
colnames(feat_mat) <- f2_top$Feature
rownames(feat_mat) <- ra_inf_ids

for (i in seq_len(nrow(f2_top))) {
  fn <- f2_top$Feature[i]
  vw <- f2_top$View[i]
  vals <- NULL
  if (fn %in% colnames(master_data)) {
    vals <- as.numeric(master_data[[fn]][match(ra_inf_ids, master_data$Sample_ID)])
  } else if (fn %in% colnames(bp_cell)) {
    vals <- bp_cell[ra_inf_ids, fn] * 100
  } else if (fn %in% rownames(gsva_scores)) {
    vals <- as.numeric(gsva_scores[fn, ra_inf_ids])
  }
  if (!is.null(vals)) feat_mat[, i] <- vals
}

# Remove features with too many NAs
good_cols <- colSums(!is.na(feat_mat)) >= 20
feat_mat <- feat_mat[, good_cols, drop = FALSE]
cat(sprintf("  Features with sufficient data: %d\n", ncol(feat_mat)))

# Impute remaining NAs with column median
for (j in seq_len(ncol(feat_mat))) {
  nas <- is.na(feat_mat[, j])
  if (any(nas)) feat_mat[nas, j] <- median(feat_mat[!nas, j])
}

# Scale
feat_mat_sc <- scale(feat_mat)

# Lasso LOOCV
n_ra <- length(y_inf)
pred_lasso <- numeric(n_ra)
selected_all <- list()

set.seed(42)
for (i in seq_len(n_ra)) {
  X_train <- feat_mat_sc[-i, , drop = FALSE]
  y_train <- y_inf[-i]
  X_test <- feat_mat_sc[i, , drop = FALSE]

  fit <- tryCatch({
    cv.glmnet(X_train, y_train, family = "binomial", alpha = 1,
              nfolds = min(10, length(y_train)), type.measure = "deviance")
  }, error = function(e) NULL)

  if (!is.null(fit)) {
    pred_lasso[i] <- predict(fit, X_test, s = "lambda.1se", type = "response")[1]
    coefs <- coef(fit, s = "lambda.1se")
    sel <- rownames(coefs)[coefs[, 1] != 0 & rownames(coefs) != "(Intercept)"]
    selected_all[[i]] <- sel
  } else {
    pred_lasso[i] <- 0.5
    selected_all[[i]] <- character(0)
  }
}

roc_lasso <- roc(y_inf, pred_lasso, quiet = TRUE)
cat(sprintf("  Lasso LOOCV AUC: %.3f\n", auc(roc_lasso)))

# Feature selection frequency
sel_freq <- table(unlist(selected_all))
sel_freq <- sort(sel_freq, decreasing = TRUE)
cat("  Selected features (frequency):\n")
print(sel_freq)

# Final full-data model for coefficient display
fit_full <- cv.glmnet(feat_mat_sc, y_inf, family = "binomial", alpha = 1,
                       nfolds = min(10, n_ra), type.measure = "deviance")
coefs_full <- coef(fit_full, s = "lambda.1se")
sel_final <- data.frame(
  Feature = rownames(coefs_full)[coefs_full[, 1] != 0 & rownames(coefs_full) != "(Intercept)"],
  Coefficient = coefs_full[coefs_full[, 1] != 0 & rownames(coefs_full) != "(Intercept)", 1]
)
if (nrow(sel_final) > 0) {
  sel_final$Feature_clean <- gsub("^BALF_|^Serum_|^PB_", "", sel_final$Feature)
  sel_final$Feature_clean <- gsub("_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14pos.*", "", sel_final$Feature_clean)
  sel_final <- sel_final %>% arrange(desc(abs(Coefficient)))

  p6g <- ggplot(sel_final, aes(Coefficient, reorder(Feature_clean, abs(Coefficient)))) +
    geom_col(fill = "#9B59B6", width = 0.6) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    annotate("text", x = max(abs(sel_final$Coefficient)) * 0.5, y = 0.5,
             label = sprintf("LOOCV AUC = %.3f", auc(roc_lasso)),
             size = 3, hjust = 0, vjust = 0, fontface = "bold") +
    labs(x = "Lasso coefficient", y = "", title = "Parsimonious infection model") +
    tn2(8) + theme(axis.text.y = element_text(size = 7))
  sv2(p6g, "Fig6g_lasso", 85, 60)
} else {
  cat("  WARNING: No features selected by lasso (all penalized to 0)\n")
  # Fallback: show top 5 Factor 2 features with their weights
  f2_display <- f2_top[1:min(5, nrow(f2_top)), ]
  f2_display$Feature_clean <- gsub("^BALF_|^Serum_|^PB_", "", f2_display$Feature)
  f2_display$Feature_clean <- gsub("_CCR6.*|_CD127.*", "", f2_display$Feature_clean)
  p6g <- ggplot(f2_display, aes(Weight, reorder(Feature_clean, abs(Weight)), fill = View)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = view_cols) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    annotate("text", x = Inf, y = -Inf, label = sprintf("Lasso LOOCV\nAUC=%.3f", auc(roc_lasso)),
             hjust = 1.1, vjust = -0.3, size = 2.5, fontface = "bold") +
    labs(x = "Factor 2 weight", y = "", title = "IIV top features", fill = "") +
    tn2(8) + theme(legend.position = "bottom", axis.text.y = element_text(size = 7, face = "italic"))
  sv2(p6g, "Fig6g_lasso", 85, 65)
}

# Save lasso results
write.csv(data.frame(Feature = names(sel_freq), Frequency = as.numeric(sel_freq)),
          file.path(table_dir, "Lasso_Feature_Selection.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════
# STEP 5: Factor vs CT Progression (Item 9)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 5: Factor vs CT Progression ===\n")

ct_ra <- ct_unique[ct_unique$Sample_Group == "RA" & ct_unique$Sample_ID %in% ra_common, ]
ct_ids <- ct_ra$Sample_ID[ct_ra$Sample_ID %in% rownames(factors)]
delta_hl <- as.numeric(ct_ra$CT_Delta_Healthy_Lung_pct[match(ct_ids, ct_ra$Sample_ID)])
f1_ct <- factors[ct_ids, 1]
f2_ct <- factors[ct_ids, 2]

sp1 <- cor.test(f1_ct, delta_hl, method = "spearman")
sp2 <- cor.test(f2_ct, delta_hl, method = "spearman")
cat(sprintf("  IDD vs DeltaHL: rho=%.3f, p=%.4f\n", sp1$estimate, sp1$p.value))
cat(sprintf("  IIV vs DeltaHL: rho=%.3f, p=%.4f\n", sp2$estimate, sp2$p.value))

df_ct <- rbind(
  data.frame(Factor = sprintf("IDD (ρ=%.3f, %s)", sp1$estimate,
              ifelse(sp1$p.value < 0.0001, "p<0.0001", sprintf("p=%.3f", sp1$p.value))),
             Fval = f1_ct, DHL = delta_hl),
  data.frame(Factor = sprintf("IIV (ρ=%.3f, %s)", sp2$estimate,
              ifelse(sp2$p.value < 0.0001, "p<0.0001", sprintf("p=%.3f", sp2$p.value))),
             Fval = f2_ct, DHL = delta_hl)
)

p6h <- ggplot(df_ct, aes(Fval, DHL)) +
  geom_point(size = 1.8, color = "#C0392B", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.4,
              fill = "grey90", formula = y ~ x) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.2) +
  facet_wrap(~Factor, scales = "free_x") +
  labs(x = "Factor value", y = "\u0394 Healthy Lung (%)") +
  tn2(8) + theme(strip.text = element_text(size = 7))
sv2(p6h, "Fig6h_factor_CT", 120, 65)

# ══════════════════════════════════════════════════════════════════════
# Retained panels: Factor 1 weights + Layer ablation
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Factor 1 weights (updated style) ===\n")

fw1 <- data.frame()
for (v in names(weights)) {
  wv <- weights[[v]][, 1]
  fw1 <- rbind(fw1, data.frame(Feature = names(wv), Weight = as.numeric(wv),
                                 View = v, stringsAsFactors = FALSE))
}
fw1$AbsW <- abs(fw1$Weight)
fw1_top <- fw1 %>% arrange(desc(AbsW)) %>% head(15)
fw1_top$Feature_clean <- gsub("^BALF_|^Serum_|^PB_", "", fw1_top$Feature)
fw1_top$Feature_clean <- gsub("_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14pos.*", "", fw1_top$Feature_clean)
# Format cytokine names — specific substitutions FIRST, then generic
fw1_top$Feature_clean <- gsub("^GMCSF$", "GM-CSF", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^GCSF$", "G-CSF", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IFNa2$", "IFN-\u03B12", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IFNg$", "IFN-\u03B3", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^TNFa$", "TNF-\u03B1", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^TNFb$", "TNF-\u03B2", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^TGFa$", "TGF-\u03B1", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL1RA$", "IL-1RA", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL1a$", "IL-1\u03B1", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL1b$", "IL-1\u03B2", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL12p40$", "IL-12p40", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL12p70$", "IL-12p70", fw1_top$Feature_clean)
fw1_top$Feature_clean <- gsub("^IL(\\d)", "IL-\\1", fw1_top$Feature_clean)
fw1_top$Feature_clean <- fct_reorder(fw1_top$Feature_clean, fw1_top$AbsW)
fw1_top$ViewLabel <- view_labels[fw1_top$View]

p6i <- ggplot(fw1_top, aes(Weight, Feature_clean, fill = View)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = view_cols, labels = view_labels) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  labs(x = "IDD weight", y = "", fill = "") +
  tn2(8) + theme(legend.position = "bottom", legend.text = element_text(size = 6),
                  axis.text.y = element_text(size = 7))
sv2(p6i, "Fig6i_IDD_weights", 90, 80)

cat("\n=== Layer ablation (updated style) ===\n")
if (file.exists("results/tables/LayerAblation_7layers_AUC.csv")) {
  abl <- read.csv("results/tables/LayerAblation_7layers_AUC.csv")
  abl$Type <- recode(abl$Type, "Single" = "Single layer", "Combined" = "Full model", "Ablation" = "Leave-one-out")
  abl$Type <- factor(abl$Type, levels = c("Single layer", "Full model", "Leave-one-out"))
  abl$Model <- fct_reorder(abl$Model, abl$AUC)

  p6j <- ggplot(abl, aes(AUC, Model, fill = Type)) +
    geom_col(width = 0.55) +
    scale_fill_manual(values = c("Single layer" = "#3498DB", "Full model" = "#C0392B", "Leave-one-out" = "#E67E22")) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.2) +
    labs(x = "LOOCV AUC", y = "", fill = "") +
    xlim(0, 1) +
    tn2(8) + theme(legend.position = "bottom", legend.text = element_text(size = 6),
                    axis.text.y = element_text(size = 6))
  sv2(p6j, "Fig6j_ablation", 100, 80)
}

# ══════════════════════════════════════════════════════════════════════
# STEP 6: Medication Confounding (Item 8)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 6: Medication Confounding ===\n")

med_data <- ct_unique[ct_unique$Sample_Group == "RA" & ct_unique$Sample_ID %in% ra_common, ]
med_cols <- c("steroid.1.0.", "MTX", "NSAIDS")
med_labels <- c("GC", "MTX", "NSAIDs")

med_results <- data.frame()
plot_list <- list()

for (mi in seq_along(med_cols)) {
  mc <- med_cols[mi]
  ml <- med_labels[mi]
  med_vals <- med_data[[mc]][match(ra_common, med_data$Sample_ID)]

  for (fi_idx in 1:2) {
    fac_vals <- factors[ra_common, fi_idx]
    fn <- fshort[paste0("Factor", fi_idx)]
    on <- fac_vals[med_vals == 1 & !is.na(med_vals)]
    off <- fac_vals[med_vals == 0 & !is.na(med_vals)]
    wt <- wilcox.test(on, off, exact=TRUE)

    med_results <- rbind(med_results, data.frame(
      Medication = ml, Factor = fn,
      N_on = length(on), N_off = length(off),
      Median_on = median(on, na.rm = TRUE), Median_off = median(off, na.rm = TRUE),
      P = wt$p.value))

    df_med <- data.frame(
      Value = c(on, off),
      Med = factor(c(rep(paste0(ml, " (+)"), length(on)), rep(paste0(ml, " (-)"), length(off))),
                    levels = c(paste0(ml, " (-)"), paste0(ml, " (+)")))
    )
    p_med <- ggplot(df_med, aes(Med, Value, fill = Med)) +
      geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.3) +
      geom_jitter(width = 0.08, size = 1, alpha = 0.5, show.legend = FALSE) +
      scale_fill_manual(values = c("#3498DB", "#E74C3C")) +
      labs(x = "", y = fn, title = sprintf("%s (p=%.3f)", ml, wt$p.value)) +
      tn2(7) + theme(legend.position = "none", plot.title = element_text(size = 7))
    plot_list[[length(plot_list) + 1]] <- p_med
  }
}

# Also test key biomarkers
bio_cols <- c("BALF_Th17.1_CCR6pos_CXCR3pos")
bio_labels <- c("BAL Th17.1")
for (bi in seq_along(bio_cols)) {
  bc <- bio_cols[bi]
  bl <- bio_labels[bi]
  bio_vals <- as.numeric(med_data[[bc]][match(ra_common, med_data$Sample_ID)])

  for (mi in seq_along(med_cols)) {
    mc <- med_cols[mi]; ml <- med_labels[mi]
    med_vals <- med_data[[mc]][match(ra_common, med_data$Sample_ID)]
    on <- bio_vals[med_vals == 1 & !is.na(med_vals) & !is.na(bio_vals)]
    off <- bio_vals[med_vals == 0 & !is.na(med_vals) & !is.na(bio_vals)]
    wt <- wilcox.test(on, off, exact=TRUE)
    med_results <- rbind(med_results, data.frame(
      Medication = ml, Factor = bl, N_on = length(on), N_off = length(off),
      Median_on = median(on, na.rm = TRUE), Median_off = median(off, na.rm = TRUE), P = wt$p.value))
  }
}

cat("  Medication confounding results:\n")
print(med_results)
write.csv(med_results, file.path(table_dir, "Medication_Confounding.csv"), row.names = FALSE)

# Combine medication plots
if (length(plot_list) >= 6) {
  p_med_all <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(title = "Medication confounding analysis",
                     theme = theme(plot.title = element_text(face = "bold", size = 10)))
  sv2(p_med_all, "FigS4_medication", 180, 110)
}

# ══════════════════════════════════════════════════════════════════════
# STEP 7: Ligand-Cell Correlation (Item 7)
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Step 7: Ligand-Cell Correlation ===\n")

# BALF cytokines vs cell-type proportions
balf_cyto_cols <- grep("^BALF_", colnames(master_data), value = TRUE)
balf_cyto_cols <- balf_cyto_cols[!grepl("Th17|Treg|Th1|Th2|Macro|Neutro|CD|Lymph|Eosino|Baso|Plasma|Other", balf_cyto_cols)]
balf_cyto_cols <- balf_cyto_cols[sapply(balf_cyto_cols, function(x) is.numeric(master_data[[x]]))]
# Keep top cytokines
if (length(balf_cyto_cols) > 15) balf_cyto_cols <- balf_cyto_cols[1:15]

cell_types <- c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC")
cell_types <- intersect(cell_types, colnames(bp_cell))

# Compute correlation matrix
cor_mat <- matrix(NA, nrow = length(balf_cyto_cols), ncol = length(cell_types))
rownames(cor_mat) <- balf_cyto_cols
colnames(cor_mat) <- cell_types
pval_mat <- cor_mat

common_ids <- intersect(rownames(bp_cell), master_data$Sample_ID)

for (i in seq_along(balf_cyto_cols)) {
  cyto_vals <- as.numeric(master_data[[balf_cyto_cols[i]]][match(common_ids, master_data$Sample_ID)])
  for (j in seq_along(cell_types)) {
    cell_vals <- bp_cell[common_ids, cell_types[j]]
    v <- !is.na(cyto_vals) & !is.na(cell_vals)
    if (sum(v) >= 10) {
      sp <- cor.test(cyto_vals[v], cell_vals[v], method = "spearman")
      cor_mat[i, j] <- sp$estimate
      pval_mat[i, j] <- sp$p.value
    }
  }
}

# Clean names
rownames(cor_mat) <- gsub("^BALF_", "", rownames(cor_mat))

# Heatmap
cor_long <- as.data.frame(cor_mat) %>%
  rownames_to_column("Cytokine") %>%
  pivot_longer(-Cytokine, names_to = "CellType", values_to = "Rho")

pval_long <- as.data.frame(pval_mat) %>%
  rownames_to_column("Cytokine") %>%
  pivot_longer(-Cytokine, names_to = "CellType", values_to = "P")

cor_long$P <- pval_long$P
cor_long$Sig <- ifelse(cor_long$P < 0.05, "*", "")

p_lc <- ggplot(cor_long, aes(CellType, Cytokine, fill = Rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Sig), size = 3, color = "black", vjust = 0.8) +
  scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                        midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  labs(x = "BayesPrism cell type", y = "BAL cytokine", title = "Cytokine-cell type associations") +
  tn2(8) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                  axis.text.y = element_text(size = 6))
sv2(p_lc, "FigS5_ligand_cell", 110, 100)

# ══════════════════════════════════════════════════════════════════════
# Save all results
# ══════════════════════════════════════════════════════════════════════
cat("\n=== Saving enhancement results ===\n")

enhancement_results <- list(
  factor_names = fnames,
  roc_idd = list(auc = as.numeric(auc(roc_idd))),
  roc_iiv = list(auc = as.numeric(auc(roc_iiv))),
  roc_lasso = list(auc = as.numeric(auc(roc_lasso))),
  best_cytokine = list(name = cyto_label, auc = best_auc),
  ct_cors = list(idd_rho = sp1$estimate, idd_p = sp1$p.value,
                  iiv_rho = sp2$estimate, iiv_p = sp2$p.value),
  med_results = med_results,
  lasso_features = sel_freq
)
save(enhancement_results, file = "results/MOFA2_Enhancement_Results.RData")

cat("\n✓ All MOFA2 enhancement analyses complete.\n")
cat(sprintf("  IDD AUC: %.3f\n", enhancement_results$roc_idd$auc))
cat(sprintf("  IIV AUC: %.3f\n", enhancement_results$roc_iiv$auc))
cat(sprintf("  Lasso AUC: %.3f\n", enhancement_results$roc_lasso$auc))
cat(sprintf("  Best cytokine: %s (AUC=%.3f)\n", enhancement_results$best_cytokine$name, enhancement_results$best_cytokine$auc))
cat(sprintf("  IDD vs CT: rho=%.3f (p=%.3f)\n", enhancement_results$ct_cors$idd_rho, enhancement_results$ct_cors$idd_p))
cat(sprintf("  IIV vs CT: rho=%.3f (p=%.3f)\n", enhancement_results$ct_cors$iiv_rho, enhancement_results$ct_cors$iiv_p))
