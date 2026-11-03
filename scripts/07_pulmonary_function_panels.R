# 07_pulmonary_function_panels.R -- Figure 4 and Figure 5g-i panels.
# Produces: panels/Fig4*, Fig5g-i .png and .pdf, panel_data_pft.csv, panel_values_pft.csv.
# Inputs: planA_data_RA24.csv, master_data_with_CT_all.csv, 06, 05 and 09 outputs.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(pROC) }); set.seed(42)
R0 <- file.path(ART_TABLES, "n31"); PD <- file.path(ART_FIGS, "n31", "panels"); dir.create(PD, recursive = TRUE, showWarnings = FALSE)
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); num <- function(x) suppressWarnings(as.numeric(x))
ms <- read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
ms <- check_english_columns(ms)
lym <- num(ms[["BAL_lymphocytes_per_mL"]]); stopifnot(all(abs(lym - num(ms[["BAL_total_cells_per_mL"]]) * num(ms$BALF_Lymphocyte_Percent) / 100) < 1))
cd4c <- lym * num(ms$CD4) / 100
d <- data.frame(id = ra$Sample_ID, dFVC = ra$d_pFVC, dDLCO = ra$d_pDLCO,
  t171c = cd4c * num(ms$`BALF_Th17.1_CCR6pos_CXCR3pos`) / 100, t171p = num(ms$`BALF_Th17.1_CCR6pos_CXCR3pos`),
  th1c = cd4c * num(ms$BALF_Th1_CCR6neg_CXCR3pos) / 100, pbth17 = num(ms$PB_Th17_CCR6pos_CXCR3neg), pbth171 = num(ms$`PB_Th17.1_CCR6pos_CXCR3pos`))
f1 <- read.csv(file.path(R0, "mofa_RA/SourceData_Fig5_n31.csv")); d$F1 <- f1$RA_Factor1[match(d$id, f1$Sample_ID)]   # RA Factor 1 of 05_RA_only_MOFA2.R
d$decline <- ifelse(is.finite(d$dFVC), as.integer(d$dFVC <= -5), NA)
d$relFVC <- 100 * (ra$PFT2_pFVC - ra$PFT1_pFVC) / ra$PFT1_pFVC; stopifnot(all(abs((ra$PFT2_pFVC - ra$PFT1_pFVC) - d$dFVC) < 1e-6, na.rm = TRUE))
d$inbuild <- ifelse(is.finite(d$relFVC), as.integer(d$relFVC <= -10), NA)
op <- read.csv(file.path(R0, "opacity/opacity_per_patient.csv")); d$dOpa <- round(op$d_opa[match(d$id, op$Sample_ID)], 2)  # recorded to 2 decimals: keep real ties
d$T2opa <- round(100 - ra$T2_aer, 2); d$pFVC2 <- ra$PFT2_pFVC; d$T1opa <- round(100 - ra$T1_aer, 2); d$pFVC1 <- ra$PFT1_pFVC   # baseline (Figure 4g)
 d$dDense <- round(ra$CT_Delta_Dense_Tissue_pct, 2); stopifnot(sum(is.finite(d$dOpa)) == 24)
write.csv(d, file.path(PD, "panel_data_pft.csv"), row.names = FALSE)
SZ <- list(F4a1 = c(520, 456), F4a2 = c(520, 456), F4a3 = c(520, 456), F4b = c(944, 826), F4c = c(944, 826), F4d = c(582, 560),
           F4e = c(1165, 561),
           F4g = c(771, 677), F4h = c(771, 677), F4i = c(771, 677), F4j = c(771, 677),
           F5g = c(858, 801), F5i = c(944, 885), F5h = c(858, 801))
LWS <- list(F4a1 = 0.3, F4a2 = 0.3, F4a3 = 0.3, F4b = 0.3 / 0.520, F4c = 0.3 / 0.520, F4g = 0.3 / 0.767, F4h = 0.3 / 0.767, F4i = 0.3 / 0.767, F4j = 0.3 / 0.767, F5g = 0.3 / 0.925, F5h = 0.3 / 0.925)   # regression line width
th <- function(bs = 13, bw = 0.4) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = bw), axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 2), "mm"))
sv <- function(p, nm) { s <- SZ[[nm]]; f <- sub("F", "Fig", nm)
  ggsave(file.path(PD, paste0(f, ".png")), p, width = s[1], height = s[2], units = "px", dpi = 300, bg = "white")
  ggsave(file.path(PD, paste0(f, ".pdf")), p, width = s[1] / 300 * 25.4, height = s[2] / 300 * 25.4, units = "mm", device = cairo_pdf) }
sp <- function(x, y) { ok <- is.finite(x) & is.finite(y); r <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = TRUE)); c(rho = unname(r$estimate), p = r$p.value, n = sum(ok)) }
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))  # figure format = manuscript text
rf <- function(r) sub("-", "\u2212", sprintf("ρ = %.2f", ifelse(abs(r) < 0.005, 0, r)), fixed = TRUE)  # true minus sign; no "-0.00"
lab <- function(s) paste0(rf(s["rho"]), "\n", pf(s["p"]))  # manuscript format
RED <- "#C0392B"; BLU <- "#2980B9"; ORA <- "#E67E22"; GRN <- "#27AE60"
chk <- list()
scp <- function(x, y, xl, yl, nm, bs = 13, logx = FALSE, right = FALSE, pcol = RED, br = FALSE, brk = NULL, xh = 0.5, ps = 2.4, k = 1, hz = TRUE, lw = LWS[[nm]]) { s <- sp(d[[x]], d[[y]]); chk[[nm]] <<- s
  dd <- d[is.finite(d[[x]]) & is.finite(d[[y]]), ]
  p <- ggplot(dd, aes(.data[[x]], .data[[y]])) + (if (hz) geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3 * k) else NULL) +  # hz = FALSE for level (not change) outcomes
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey30", linewidth = lw) +
    geom_point(color = pcol, size = ps, alpha = 0.75) +
    annotate("text", x = if (right || br) Inf else if (logx) min(dd[[x]]) else -Inf, y = if (br) -Inf else Inf, label = lab(s), hjust = if (right || br) 1.08 else if (logx) 0 else -0.08, vjust = if (br) -0.25 else 1.15, size = bs / 3.2, lineheight = 0.9) + labs(x = xl, y = yl) + th(bs, 0.4 * k)  # -Inf is NaN on a log axis
  if (logx) p <- p + (if (is.null(brk)) scale_x_log10(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)) else scale_x_log10(breaks = brk, labels = parse(text = paste0("10^", log10(brk)))))
  if (bs != 13) p <- p + theme(axis.text = element_text(size = 11 * bs / 13, color = "black"), axis.title = element_text(size = 11 * bs / 13, face = "bold"))
  if (xh != 0.5) p <- p + theme(axis.title.x = element_text(hjust = xh))  # keep a long x title inside a small panel
  sv(p, nm) }
scp("t171c", "dFVC", "BAL Th17.1 (cells/mL)", "Δ%FVC (points)", "F4a1", bs = 8.02, ps = 1.48, k = 0.617, logx = TRUE)
scp("t171p", "dFVC", "BAL Th17.1 (% of CD4+)", "Δ%FVC (points)", "F4a2", bs = 8.02, ps = 1.48, k = 0.617, br = TRUE)
scp("th1c", "dFVC", "BAL Th1 (cells/mL)", "Δ%FVC (points)", "F4a3", bs = 8.02, ps = 1.48, k = 0.617, logx = TRUE, brk = c(100, 1e3, 1e4, 1e5))  # label top-left (bottom-right touched a point)
scp("pbth171", "dFVC", "PB Th17.1 (% of CD4+)", "Δ%FVC (points)", "F4b", br = TRUE)
scp("t171c", "dDLCO", "BAL Th17.1 (cells/mL)", "Δ%DLCO (points)", "F4c", logx = TRUE)
rocs <- function(vars, cols, nm, title, y = "decline", bs = 12, lt = 7, ts = 11, lw = 0.9, lpos = c(0.98, 0.03), textleg = FALSE) { ok <- !is.na(d[[y]]); L <- list(); rl <- list()
  for (v in names(vars)) { r <- roc(d[[y]][ok], d[[vars[v]]][ok], levels = c(0, 1), direction = "auto", quiet = TRUE); set.seed(42)
    ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, quiet = TRUE); chk[[paste0(nm, "_", v)]] <<- c(AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3], dir = r$direction == "<")
    L[[v]] <- sprintf("%s (AUC = %.3f)", v, auc(r)); rl[[v]] <- data.frame(v = v, fpr = 1 - r$specificities, tpr = r$sensitivities) }
  rd <- do.call(rbind, rl); rd$v <- factor(rd$v, levels = names(vars))
  p <- ggplot(rd, aes(fpr, tpr, color = v)) + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") + geom_path(linewidth = lw) +
    scale_color_manual(values = cols, labels = unlist(L[names(vars)]), name = NULL) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 – Specificity", y = "Sensitivity", title = title) + th(bs) +
    theme(legend.position = lpos, legend.justification = c(1, 0), legend.text = element_text(size = lt), legend.key.size = unit(if (lt < 6) 2 else 3, "mm"), legend.key.spacing.y = unit(0, "mm"),
          legend.background = element_rect(fill = alpha("white", 0.9), color = NA), plot.title = element_text(face = "bold", size = ts, hjust = 0.5), plot.margin = unit(if (bs < 9) c(0.5, 1, 0.5, 0.5) else c(2, 3, 2, 2), "mm"))
  if (bs < 9) p <- p + theme(axis.text = element_text(size = bs, color = "black"), axis.title = element_text(size = bs + 0.5, face = "bold"))  # explicit: th replaces the text element
  if (textleg) {
    k <- names(vars); lab2 <- sub(" \\(AUC = ", "\nAUC = ", sub("\\)$", "", unlist(L[k])))
    p <- p + guides(color = "none") + scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) + scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      annotate("text", x = 0.99, y = 0.02 + (rev(seq_along(k)) - 1) * 0.15, label = lab2, color = cols[seq_along(k)], hjust = 1, vjust = 0, size = lt / .pt, lineheight = 0.85, fontface = "bold")
  }
  sv(p, nm) }
rocs(c(`BAL Th17.1 count` = "t171c", `BAL Th17.1 %` = "t171p", `PB Th17` = "pbth17"), c(ORA, RED, BLU), "F4d", "Relative %FVC decline ≥ 10%", y = "inbuild", bs = 6, lt = 4.5, ts = 6, lw = 0.55, textleg = TRUE)  # INBUILD threshold; 5 points in the text only
S <- read.csv(file.path(R0, "common_panel/dFVC_panel.csv")); S <- S[is.finite(S$P), ]   # 113 of the 114 panel markers tested
S$nl <- -log10(S$P); S$cat <- S$Category; S$cat[grepl("^BAL FCM|BAL differential|BAL lymphocyte|BAL CD4", S$cat)] <- "BAL cells"; S$cat[S$cat %in% c("Microbiome diversity")] <- "Microbiome"
nice <- c(BAL_Th17.1_count = "BAL Th17.1 count", BAL_Th1_count = "BAL Th1 count", Serum_IL10 = "Serum IL-10", Serum_TNFb = "Serum TNF-β", BALF_IFNa2 = "BAL IFN-α2", Serum_IL22 = "Serum IL-22")
S <- S[order(S$P), ]; S$lab <- ifelse(S$Marker %in% c("BAL_Th17.1_count", "BAL_Th1_count", "Serum_IL10", "Serum_TNFb", "BALF_IFNa2"), nice[S$Marker], NA); chk$F4e <- c(n = nrow(S), minBH = min(S$BH), nominal = sum(S$P < 0.05))
p <- ggplot(S, aes(rho, nl)) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) + geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = cat), size = 1.4, alpha = 0.8) + geom_text_repel(data = S[S$rho < 0, ], aes(label = lab), size = 2.1, na.rm = TRUE, seed = 42, min.segment.length = 0, box.padding = 0.4, direction = "y", nudge_x = 0.22, xlim = c(-0.45, -0.02), max.overlaps = Inf) +
  geom_text_repel(data = S[S$rho > 0, ], aes(label = lab), size = 2.1, na.rm = TRUE, seed = 42, min.segment.length = 0, box.padding = 0.4, direction = "y", nudge_x = -0.25, xlim = c(0.02, 0.55), max.overlaps = Inf) +
  scale_color_manual(values = c(`BAL cells` = "#1B9E77", `BAL cytokine` = "#D95F02", Microbiome = "#7570B3", `PB FCM` = "#66A61E", `Serum cytokine` = "#E6AB02"), name = NULL) + scale_x_continuous(limits = c(-0.65, 0.7), breaks = c(-0.5, 0, 0.5)) + annotate("text", x = 0.69, y = 0.33, hjust = 1, vjust = 0.5, size = 2.1, label = { stopifnot(nrow(S) == 113); sprintf("%d markers tested\nlowest BH-\nadjusted p = %.2f", nrow(S), min(S$BH)) }, lineheight = 0.85) +
  labs(x = "Spearman ρ with Δ%FVC", y = expression(bold(-log[10]~p))) + th(11) + theme(axis.text = element_text(size = 6.8, color = "black"), axis.title = element_text(size = 7, face = "bold"), legend.position = "right", legend.text = element_text(size = 6), legend.key.size = unit(2.5, "mm"), legend.key.spacing.y = unit(0.3, "mm"), plot.margin = unit(c(1, 1, 1, 1), "mm")) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.8)))
sv(p, "F4e")
scp("T1opa", "pFVC1", "Total opacity at baseline (% of lung)", "%FVC at baseline", "F4g", bs = 10.46, ps = 1.48, k = 0.617, pcol = BLU, right = TRUE, hz = FALSE)
scp("T2opa", "pFVC2", "Total opacity at follow-up (% of lung)", "%FVC at follow-up", "F4h", bs = 10.46, ps = 1.48, k = 0.617, pcol = BLU, right = TRUE, hz = FALSE)
scp("dDense", "dFVC", "ΔDense Tissue (% of lung)", "Δ%FVC (points)", "F4i", bs = 10.46, ps = 1.48, k = 0.617, pcol = BLU, right = TRUE)
scp("t171c", "dOpa", "BAL Th17.1 (cells/mL)", "ΔTotal opacity (% of lung)", "F4j", bs = 10.46, ps = 1.48, k = 0.617, logx = TRUE, br = TRUE)
scp("F1", "dFVC", "RA Factor 1", "Δ%FVC (points)", "F5g", right = TRUE)
scp("F1", "dOpa", "RA Factor 1", "ΔTotal opacity (% of lung)", "F5h", right = TRUE)
rocs(c(`BAL Th17.1 count` = "t171c", `RA Factor 1` = "F1"), c(ORA, BLU), "F5i", "Relative %FVC decline ≥ 10% (RA)", y = "inbuild")  # INBUILD criterion as Figure 4d
# Results (no panel): BAL Th17.1 count for an absolute %FVC decline of 5 points or more (progressive pulmonary fibrosis criterion)
ok <- !is.na(d$decline) & is.finite(d$t171c); r5 <- roc(d$decline[ok], d$t171c[ok], levels = c(0, 1), direction = "auto", quiet = TRUE); set.seed(42)
ci5 <- ci.auc(r5, method = "bootstrap", boot.n = 2000, quiet = TRUE); chk[["Decline5_BAL Th17.1 count"]] <- c(AUC = as.numeric(auc(r5)), lo = ci5[1], hi = ci5[3], n_decline = sum(d$decline[ok] == 1), n_other = sum(d$decline[ok] == 0))
stopifnot(abs(as.numeric(auc(r5)) - 0.83) < 5e-3, sum(d$decline[ok] == 1) == 10)   # value of the article
out <- do.call(rbind, lapply(names(chk), function(k) data.frame(panel = k, stat = names(chk[[k]]), value = as.numeric(chk[[k]]))))
write.csv(out, file.path(PD, "panel_values_pft.csv"), row.names = FALSE); print(out, row.names = FALSE)
# Figure 4a is one image: the three sub-panels side by side (520 px each)
a <- lapply(1:3, function(i) png::readPNG(file.path(PD, sprintf("Fig4a%d.png", i))))
strip <- array(0, c(dim(a[[1]])[1], sum(sapply(a, function(z) dim(z)[2])), dim(a[[1]])[3])); x0 <- 0
for (z in a) { strip[, x0 + seq_len(dim(z)[2]), ] <- z; x0 <- x0 + dim(z)[2] }
png::writePNG(strip, file.path(PD, "Fig4a.png"))
