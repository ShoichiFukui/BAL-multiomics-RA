# panels_infection_markers.R -- Figure 3c-f (BAL Th17.1 and Treg proportions, BAL/PB Th17.1 ratio and ROC curves for subsequent
#   respiratory infection, RA n = 24) and Supplementary Fig. 3b (BAL microbiome genera, infection + vs -).
# Produces: <ART_FIGS>/n31/panels/Fig3c-f.png/.pdf, SFig3b.png/.pdf and infection_markers_values.csv (printed P and AUC).
# Inputs: planA_data_RA24.csv, results/tables/master_data_with_CT_all.csv, microbiome/SD4_microbiome_genus_exact.csv (18).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(pROC) })
TB <- file.path(ART_TABLES, "n31"); panel_dir <- file.path(ART_FIGS, "n31", "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          legend.text = element_text(size = rel(0.8)), legend.key.size = unit(3, "mm"),
          legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33), axis.line = element_blank(),
          plot.margin = unit(c(2, 3, 2, 3), "mm"))
}
sv <- function(p, name, w = 80, h = 70) {   # jitter drawn from set.seed(42) in both files
  set.seed(42); ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  set.seed(42); ggsave(file.path(panel_dir, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))
NEG <- "Infection (−)"; POS <- "Infection (+)"
CI <- setNames(c("#3498DB", "#E74C3C"), c(NEG, POS)); CD <- setNames(c("#1A5276", "#922B21"), c(NEG, POS))

ra <- read.csv(file.path(TB, "planA/planA_data_RA24.csv"), check.names = FALSE)
ms <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE)); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
ig <- factor(ifelse(ra$respiratory_infection == 1, POS, NEG), levels = c(NEG, POS)); stopifnot(sum(ig == POS) == 6)
th171 <- as.numeric(ms$BALF_Th17.1_CCR6pos_CXCR3pos); treg <- as.numeric(ms$BALF_Treg_CD127neg_CD25pos)
ratio <- th171 / (as.numeric(ms$PB_Th17.1_CCR6pos_CXCR3pos) + 0.01)   # PB Th17.1 is 0% in two patients
vals <- list()
wt <- function(v) wilcox.test(v[ig == POS], v[ig == NEG], exact = TRUE)$p.value
bxi <- function(v, yl, nm) {
  df <- data.frame(V = v, IG = ig)[!is.na(v), ]; p <- wt(v); vals[[nm]] <<- c(P = p)
  ymx <- max(df$V); ysp <- (ymx - min(df$V)) * 0.12
  ggplot(df, aes(IG, V, fill = IG)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25, alpha = 0.7) +
    geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65, show.legend = FALSE) +
    scale_fill_manual(values = CI) + scale_color_manual(values = CD) + labs(x = "", y = yl) +
    annotate("text", x = 1.5, y = ymx + ysp * 1.2, label = pf(p), size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    tn(8) + theme(legend.position = "none", axis.text.x = element_text(size = 9.5))
}
sv(bxi(th171, "BAL Th17.1 (%)", "Fig3c"), "Fig3c")
sv(bxi(treg, "BAL Treg (%)", "Fig3d"), "Fig3d")

# Figure 3e -- BAL/PB Th17.1 ratio (log10 axis)
dfr <- data.frame(V = ratio, IG = ig); pr <- wt(ratio); vals$Fig3e <- c(P = pr)
sv(ggplot(dfr, aes(IG, V, fill = IG)) +
     geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25, alpha = 0.7) +
     geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65, show.legend = FALSE) +
     scale_fill_manual(values = CI) + scale_color_manual(values = CD) +
     scale_y_log10(labels = function(b) parse(text = paste0("10^", round(log10(b))))) + labs(x = "", y = "Th17.1 BAL/PB (log10)") +
     annotate("text", x = 1.5, y = max(dfr$V) * 2, label = pf(pr), size = 2.5) +
     tn(8) + theme(legend.position = "none", axis.text.x = element_text(size = 9.5)), "Fig3e")

# Figure 3f -- ROC curves for subsequent infection (AUC in the legend)
y <- ra$respiratory_infection; roc_df <- data.frame()
for (bm in list(list(v = th171, label = "BAL Th17.1"), list(v = treg, label = "BAL Treg"), list(v = ratio, label = "Th17.1 BAL/PB"))) {
  ok <- !is.na(bm$v); ro <- roc(y[ok], bm$v[ok], quiet = TRUE, direction = "auto"); a <- as.numeric(auc(ro)); vals[[paste("Fig3f", bm$label)]] <- c(AUC = a)
  x <- 1 - ro$specificities; s <- ro$sensitivities; o <- order(x, s); x <- x[o]; s <- s[o]
  if (x[1] != 0 || s[1] != 0) { x <- c(0, x); s <- c(0, s) }
  if (tail(x, 1) != 1 || tail(s, 1) != 1) { x <- c(x, 1); s <- c(s, 1) }
  roc_df <- rbind(roc_df, data.frame(Spec = x, Sens = s, Model = sprintf("%s (%.3f)", bm$label, a)))
}
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))
sv(ggplot(roc_df, aes(Spec, Sens, color = Model)) + geom_step(linewidth = 0.6, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#E74C3C", "#3498DB", "#27AE60")) + labs(x = "1 – Specificity", y = "Sensitivity", color = "") + coord_equal() +
     tn(8) + theme(legend.position = c(0.62, 0.22), legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm")), "Fig3f")
stopifnot(abs(vals$Fig3c[["P"]] - 0.0084) < 5e-4, abs(vals[["Fig3f BAL Th17.1"]][["AUC"]] - 0.870) < 5e-4)   # values of the article

# Supplementary Fig. 3b -- the eight genera with the lowest P (exact Wilcoxon, Supplementary Data 4)
mi <- read.csv(file.path(TB, "microbiome/SD4_microbiome_genus_exact.csv"))
tmi <- head(mi[order(mi$P_value), ], 8) %>% arrange(Log2FC); tmi$Genus <- fct_inorder(tmi$Genus)
tmi$sg <- ifelse(tmi$P_value < 0.05, "p < 0.05", "p >= 0.05")
sv(ggplot(tmi, aes(Log2FC, Genus, fill = sg)) + geom_col(width = 0.5) +
     scale_fill_manual(values = setNames(c("#C0392B", "#95A5A6"), c("p < 0.05", "p >= 0.05"))) + geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "log2FC(Infection+ vs −)", y = "", fill = "") +
     tn(7) + theme(axis.text.y = element_text(size = 6, face = "italic"), legend.position = "bottom", legend.text = element_text(size = 6)), "SFig3b", 90, 70)
write.csv(data.frame(Item = rep(names(vals), lengths(vals)), Stat = unlist(lapply(vals, names)), Value = unlist(vals)), file.path(panel_dir, "infection_markers_values.csv"), row.names = FALSE)
