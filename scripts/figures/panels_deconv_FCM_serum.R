# panels_deconv_FCM_serum.R -- Figure 1f (BayesPrism cell-fraction stacked bar), Figure 1g-i (deconvolution vs BAL differential counts),
#   Figure 2f (BayesPrism plasma cell fraction), Supplementary Fig. 1 (Bland-Altman) and Supplementary Fig. 2 (top-10 serum
#   cytokines vs the healthy 95% reference interval), n = 31 (RA 24, sarcoidosis 7).
# Produces: <ART_FIGS>/n31/panels/Fig1f, Fig1g-i, Fig2f, SFig1, SFig2 .png/.pdf and a value-check table.
# Inputs: bayesprism_n31/celltype_proportions_cellFraction_n31.csv (01b_BayesPrism_n31.R), Fig1ghi_SFig1_deconv_vs_differential.csv and
#   Fig2f_bp_cell_plasma_Bcell.csv (11_deconvolution_FCM_validation_n31.R), ST1_serum_vs_healthy.csv (= Supplementary
#   Table 1 values; 08e), FCM_integrated_data_transformed.xlsx and 41plex_healthy_control.xlsx (controlled access, under <ART_DIR>),
#   scripts/figures/panel_map.json (pixel sizes of the panels; shipped with the code).
# Every number printed on a panel is recomputed and checked against the tables (the script stops on a mismatch).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(patchwork); library(DescTools)
  library(jsonlite); library(magick)
})
cfg <- set_cohort("n31"); BPD <- file.path(ART_TABLES, "n31", "bayesprism_n31")

out_dir <- file.path(cfg$figs, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pmap <- fromJSON("scripts/figures/panel_map.json")   # pixel sizes of the figure panels (shipped with the code)

# ---- value log: every number printed on a panel ------------------------------
vlog <- list()
logv <- function(panel, item, label, table_value, recomputed) {
  ok <- isTRUE(all.equal(as.numeric(table_value), as.numeric(recomputed), tolerance = 1e-9))
  vlog[[length(vlog) + 1]] <<- data.frame(Panel = panel, Item = item, Label = label,
    Table_value = table_value, Recomputed = recomputed, Match = ok)
  if (!ok) stop(sprintf("%s %s: table %.10g vs recomputed %.10g", panel, item, table_value, recomputed))
}

# ---- save at the exact pixel size / dpi of the article panel -----------------
sv_px <- function(p, name, dpi) {
  px <- pmap[[name]]$px
  w <- px[1] / dpi; h <- px[2] / dpi
  png_f <- file.path(out_dir, paste0(name, ".png"))
  ggsave(png_f, p, width = w, height = h, units = "in", dpi = dpi, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", device = cairo_pdf)
  im <- image_read(png_f); inf <- image_info(im)
  if (inf$width != px[1] || inf$height != px[2]) {
    # the PNG must have the pixel size of the article panel (checked, not corrected)
    stop(sprintf("%s: %dx%d, expected %dx%d", name, inf$width, inf$height, px[1], px[2]))
  }
  
}

# ---- data --------------------------------------------------------------------
bp_cell <- read.csv(file.path(BPD, "celltype_proportions_cellFraction_n31.csv"), row.names = 1)
bp_cell <- bp_cell[!rownames(bp_cell) %in% cfg$EXCL, ]
check_n(rownames(bp_cell), cfg, "bp_cell")
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet = "BALF_analysis")
colnames(fcm)[1] <- "Sample_ID"

t_val  <- read.csv(file.path(BPD, "Fig1ghi_SFig1_deconv_vs_differential.csv"))
t_2g   <- read.csv(file.path(BPD, "Fig2f_bp_cell_plasma_Bcell.csv"))
t_st2  <- read.csv(file.path(cfg$tables, "ST1_serum_vs_healthy.csv"), check.names = FALSE)

CG <- c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")
CC <- c(Macrophage = "#E74C3C", T_cell = "#3498DB", NK = "#2ECC71",
        B_cell = "#9B59B6", Plasma = "#E67E22", Neutrophil = "#F1C40F",
        DC = "#8B4513", Epithelial = "#E91E63", Mast = "#95A5A6",
        Other = "#BDC3C7")

# Theme as run for the manuscript (Fig1g-i, Fig1f)
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
      legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
      axis.line = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3),
      strip.text = element_text(size = rel(0.85), face = "bold"),
      plot.margin = unit(c(2, 3, 2, 3), "mm"))
}

# Fig1f — stacked bar (as run for the manuscript)
tn_1f <- function(bs = 8) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"),
        axis.title = element_text(size = rel(0.95), face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33), axis.line = element_blank(),
        plot.margin = unit(c(1, 2, 1, 2), "mm"))
# Calibration (measured on the Fig1f.png of the article): that image has 1/2 mm
# plot margins and 2.1 mm legend keys. Re-rendered with ggplot2 4.0.2 these settings
# reproduce its frame position (x 166-1629 px) and legend rows within 5 px.
bp_l <- bp_cell %>% rownames_to_column("S") %>%
  mutate(G = sample_group(S), S_d = public_label(S)) %>%   # de-identified labels (RA#, Sarcoidosis#) from cohort_definition.csv
  pivot_longer(-c(S, G, S_d), names_to = "CT", values_to = "Fr") %>%
  mutate(S_d = factor(S_d, levels = unique(S_d[order(G == "Sarcoidosis", as.integer(gsub("\\D", "", S_d)))])))  # numeric order within RA, then sarcoidosis (RA5 before RA10)

p1f <- ggplot(bp_l %>% mutate(CT = factor(CT, levels = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"))), aes(S_d, Fr, fill = CT)) + geom_col(width = 0.85) + scale_fill_manual(values = CC, breaks = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"), labels = function(x) sub("_", " ", x)) + guides(fill = guide_legend(byrow = TRUE)) +
  labs(x = "", y = "Proportion", fill = "") +
  tn_1f(7) + scale_x_discrete(breaks = { lv <- levels(bp_l$S_d); nr <- sum(grepl("^RA", lv)); ns <- length(lv) - nr; lv[c(ceiling(nr / 2), nr + ceiling(ns / 2))] }, labels = c("RA", "Sarcoidosis")) +
  theme(axis.text.x = element_text(size = 9.9, color = "black", margin = margin(t = 2)), axis.ticks.x = element_blank(),
                   legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = element_text(size = 9),
                   legend.key.size = unit(3.15, "mm"), legend.key.spacing.y = unit(1.6, "mm"), legend.margin = margin(-2, 0, 0, 0))
sv_px(p1f, "Fig1f", 300)

# Fig1g-i — deconvolution vs the BAL differential counts (columns BALF_*_Percent of the FCM workbook; file names keep 'FCM')
cf <- intersect(rownames(bp_cell), fcm$Sample_ID)
fi_idx <- match(cf, fcm$Sample_ID)
lc <- intersect(c("T_cell", "NK", "B_cell", "Plasma"), colnames(bp_cell))

fmt_rho <- function(r) sprintf("%.3f", r)
fmt_p4  <- function(p) ifelse(p < 0.0001, "p<0.0001", sprintf("p=%.4f", p))
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))  # figure format = manuscript text
rf <- function(r) sprintf("ρ = %.2f", r)

mkf <- function(bv, fv, lb, ct, panel) {
  v <- !is.na(fv) & !is.na(bv)
  sp <- suppressWarnings(cor.test(bv[v], fv[v], method = "spearman", exact = TRUE))
  tr <- t_val[t_val$CellType == ct, ]
  logv(panel, "N", as.character(sum(v)), tr$N, sum(v))
  logv(panel, "Spearman rho", paste0("ρ=", fmt_rho(sp$estimate)), tr$Spearman_rho, unname(sp$estimate))
  logv(panel, "Spearman p", fmt_p4(sp$p.value), tr$Spearman_p, sp$p.value)
  df <- data.frame(F = fv[v], B = bv[v], G = sample_group(cf[v]))
  ggplot(df, aes(F, B, color = G)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.3 / 0.902,
                linetype = "solid", formula = y ~ x, inherit.aes = FALSE, aes(x = F, y = B)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    scale_color_manual(values = CG) +
    annotate("text", x = if (panel == "Fig1h") Inf else -Inf, y = if (panel == "Fig1h") -Inf else Inf, hjust = if (panel == "Fig1h") 1.05 else -0.05, vjust = if (panel == "Fig1h") -0.3 else 1.3, size = 3,
             label = paste0(rf(tr$Spearman_rho), "\n", pf(tr$Spearman_p))) +
    labs(x = paste("Differential", lb, "(%)"), y = paste("Deconv.", lb, "(%)"), color = "") +
    tn(7) + theme(legend.position = "bottom", legend.direction = "horizontal",
                  legend.text = element_text(size = 9.5),
                  # calibration: default 1.2-line keys reproduce the frame of the article panel
                  # (y 25-571 px) exactly under ggplot2 4.0.2
                  legend.key.size = unit(1.2, "lines"))
}
sets <- list(
  Fig1g = list(B = bp_cell[cf, "Macrophage"] * 100, F = as.numeric(fcm$BALF_Macrophage_Percent[fi_idx]), lb = "Mac", ct = "Macrophage"),
  Fig1h = list(B = rowSums(bp_cell[cf, lc]) * 100, F = as.numeric(fcm$BALF_Lymphocyte_Percent[fi_idx]), lb = "Lym", ct = "Lymphocyte"),
  Fig1i = list(B = bp_cell[cf, "Neutrophil"] * 100, F = as.numeric(fcm$BALF_Neutrophil_Percent[fi_idx]), lb = "Neu", ct = "Neutrophil"))
for (nm in names(sets)) {
  s <- sets[[nm]]
  sv_px(mkf(s$B, s$F, s$lb, s$ct, nm), nm, 300)
}

# Fig2f — plasma cells (as run for the manuscript)
tn_2g <- function(bs = 7.5) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(.9), color = "black"),
          axis.title = element_text(size = rel(.95), face = "bold"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = .3),
          axis.line = element_blank(),
          plot.margin = unit(c(1, 2, 1, 2), "mm"))
}
ga <- sample_group(rownames(bp_cell))
dp <- data.frame(P = bp_cell[, "Plasma"] * 100, G = ga)
wtp <- wilcox.test(dp$P[dp$G == "RA"], dp$P[dp$G == "Sarcoidosis"], exact = TRUE)
tp <- t_2g[t_2g$CellType == "Plasma", ]
p_label <- pf(tp$p)
logv("Fig2f", "Wilcoxon p (plasma, bp_cell)", p_label, tp$p, wtp$p.value)
logv("Fig2f", "n RA", as.character(sum(ga == "RA")), tp$n_RA, sum(ga == "RA"))
logv("Fig2f", "n Sarc", as.character(sum(ga == "Sarcoidosis")), tp$n_Sarc, sum(ga == "Sarcoidosis"))
set.seed(42)  # geom_jitter
p2g <- ggplot(dp, aes(G, P + 1e-6, fill = G)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.2) +
  geom_jitter(width = 0.07, size = 1.2, alpha = 0.5, show.legend = FALSE) +
  scale_y_log10(labels = function(b) parse(text = paste0("10^", round(log10(b))))) + scale_fill_manual(values = CG) +
  labs(x = "", y = "Plasma cell fraction (%)") +
  annotate("text", x = 1.5, y = max(dp$P) * 3, label = p_label, size = 2.5) +
  tn_2g(8) + theme(legend.position = "none")
sv_px(p2g, "Fig2f", 600)

# SFig1 — Bland-Altman a-c (as run for the manuscript)
tn_s1 <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.1), face = "bold", hjust = 0,
                                    margin = ggplot2::margin(0, 0, 3, 0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(),
          plot.margin = unit(c(2, 3, 2, 3), "mm"))
}
ba_panel <- function(D, F, lb, ct, tag) {
  v <- !is.na(D) & !is.na(F); D <- D[v]; F <- F[v]
  m <- (D + F) / 2; d <- D - F
  bias <- mean(d); s <- sd(d); lo <- bias - 1.96 * s; hi <- bias + 1.96 * s
  tr <- t_val[t_val$CellType == ct, ]
  ccc <- DescTools::CCC(D, F)$rho.c$est
  logv("SFig1", paste(tag, ct, "CCC"), sprintf("CCC = %.3f", tr$CCC), tr$CCC, ccc)
  # not printed, but drawn as lines: check against the table too
  logv("SFig1", paste(tag, ct, "BA bias (line)"), "(dashed line)", tr$BA_bias, bias)
  logv("SFig1", paste(tag, ct, "LoA lower (line)"), "(dotted line)", tr$BA_LoA_lower, lo)
  logv("SFig1", paste(tag, ct, "LoA upper (line)"), "(dotted line)", tr$BA_LoA_upper, hi)
  df <- data.frame(m = m, d = d, G = sample_group(cf[v]))
  ggplot(df, aes(m, d, color = G)) +
    geom_hline(yintercept = bias, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_hline(yintercept = c(lo, hi), linetype = "dotted", color = "#C0392B", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.8) +
    scale_color_manual(values = c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")) +
    # a = bottom-left just above the lower LoA line; c = just below the upper LoA line
    annotate("text", x = -Inf, y = switch(tag, a = lo, c = hi, Inf), hjust = -0.05, vjust = switch(tag, a = -0.6, c = 1.5, 1.3), size = 2.8,
             label = sprintf("CCC = %.3f", tr$CCC)) +
    labs(title = paste0(toupper(tag), "   ", lb), x = "Mean of methods (%)",
         y = "Deconv. − differential (%)", color = "") +
    tn_s1(7) + theme(legend.position = "none")
}
pa <- ba_panel(sets$Fig1g$B, sets$Fig1g$F, "Macrophage", "Macrophage", "a")
pb <- ba_panel(sets$Fig1h$B, sets$Fig1h$F, "Lymphocyte", "Lymphocyte", "b")
pc <- ba_panel(sets$Fig1i$B, sets$Fig1i$F, "Neutrophil", "Neutrophil", "c")
sv_px(pa + pb + pc + plot_layout(nrow = 1), "SFig1", 300)

# SFig2 — serum cytokines vs healthy 95% CI (as run for the manuscript)
tn_s2 <- function(bs = 7) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.0), face = "bold", hjust = 0.5,
                                    margin = ggplot2::margin(0, 0, 2, 0)),
          legend.text = element_text(size = rel(0.8)),
          legend.key.size = unit(3, "mm"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(),
          plot.margin = unit(c(2, 2, 2, 2), "mm"))
}
fmt_p_s2 <- function(p) if (p < 0.0001) "p < 0.0001" else sprintf("p = %.4f", p)
cyto_lab <- c(GCSF = "G-CSF", IL6 = "IL-6", IL13 = "IL-13", IL17A = "IL-17A", IFNa2 = "IFN-α2",
              IL22 = "IL-22", IL4 = "IL-4", IL2 = "IL-2", GROa = "GRO-α", FGF2 = "FGF-2",
              MCP3 = "MCP-3", TGFa = "TGF-α")
t10 <- t_st2[t_st2$in_top10_this_cohort, ]
t10 <- t10[order(-abs(t10$cliff_delta), t10$RA_vs_Sarc_p), ]   # same rule as the screen script
sel <- t10$Cytokine

stopifnot(length(sel) == 10)

e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
md <- e$master_data; rm(e); invisible(gc())
md <- md[!md$Sample_ID %in% cfg$EXCL, ]
check_n(md$Sample_ID, cfg, "master_data")
col_ra <- "#C0392B"; col_sarc <- "#2980B9"; col_healthy <- "#27AE60"

plots_s2 <- list()
for (cy in sel) {
  vals <- as.numeric(md[[paste0("Serum_", cy)]])
  grp <- ifelse(md$Sample_Group == "RA", "RA", "Sarcoidosis")
  ok <- !is.na(vals)
  df <- data.frame(Value = vals[ok], Group = factor(grp[ok], levels = c("RA", "Sarcoidosis")))
  wt <- suppressWarnings(wilcox.test(Value ~ Group, data = df, exact = TRUE))
  tr <- t_st2[t_st2$Cytokine == cy, ]
  p_lab <- pf(tr$RA_vs_Sarc_p)
  logv("SFig2", paste(cy, "RA vs Sarc p"), p_lab, tr$RA_vs_Sarc_p, wt$p.value)
  logv("SFig2", paste(cy, "n RA"), "(n)", tr$RA_n, sum(df$Group == "RA"))
  logv("SFig2", paste(cy, "n Sarc"), "(n)", tr$Sarc_n, sum(df$Group == "Sarcoidosis"))
  h_med <- tr$Healthy_median; h_lo <- tr[["Healthy_2.5"]]; h_hi <- tr[["Healthy_97.5"]]
  ymax <- max(df$Value, h_hi, na.rm = TRUE) * 1.25
  set.seed(42)  # geom_jitter
  p <- ggplot(df, aes(Group, Value, fill = Group)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = h_lo, ymax = h_hi, fill = col_healthy, alpha = 0.15) +
    geom_hline(yintercept = h_med, linetype = "dashed", color = col_healthy, linewidth = 0.4) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = c("RA" = col_ra, "Sarcoidosis" = col_sarc)) +
    labs(x = "", y = "pg/mL", title = cyto_lab[[cy]]) +
    annotate("text", x = 1.5, y = ymax, label = p_lab, size = 2.2) +
    coord_cartesian(ylim = c(0, ymax * 1.05), clip = "off") +
    tn_s2(7) + theme(legend.position = "none", axis.text.x = element_text(size = 6, angle = 0))
  plots_s2[[cy]] <- p
}
p_s2 <- wrap_plots(plots_s2, ncol = 5) +
  plot_annotation(theme = theme(plot.margin = ggplot2::margin(2, 2, 2, 2)))
sv_px(p_s2, "SFig2", 300)

# Healthy band check: ST2 percentiles vs the raw healthy xlsx (independent recomputation)
hc <- read_excel(file.path(ART_DIR, "41plex_healthy_control.xlsx"), sheet = "FIXED")
hmap <- c(GCSF = "G-CSF", IL6 = "IL-6", IL13 = "IL-13", IL17A = "IL-17A", IFNa2 = "IFNa2", IL22 = "IL-22",
          IL4 = "IL-4", IL2 = "IL-2", GROa = "GROa", FGF2 = "FGF-2", MCP3 = "MCP-3", TGFa = "TGFa")
for (cy in sel) {
  hv <- as.numeric(hc[[hmap[[cy]]]]); hv <- hv[!is.na(hv)]
  q <- unname(quantile(hv, c(0.025, 0.975))); tr <- t_st2[t_st2$Cytokine == cy, ]
  logv("SFig2", paste(cy, "healthy 2.5th pct (band)"), "(band)", tr[["Healthy_2.5"]], q[1])
  logv("SFig2", paste(cy, "healthy 97.5th pct (band)"), "(band)", tr[["Healthy_97.5"]], q[2])
  logv("SFig2", paste(cy, "healthy median (line, raw)"), "(line)", tr$Healthy_median, median(hv))
}

# ---- value report ------------------------------------------------------------
vl <- do.call(rbind, vlog)
write.csv(vl, file.path(out_dir, "panels_deconv_FCM_serum_values.csv"), row.names = FALSE)

print(vl[vl$Label != "(line)" & vl$Label != "(band)" & vl$Label != "(n)" &
           !grepl("^\\(", vl$Label), c("Panel", "Item", "Label", "Match")], row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "panels_deconv_FCM_serum_sessionInfo.txt"))
