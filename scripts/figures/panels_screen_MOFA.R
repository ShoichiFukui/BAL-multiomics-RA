# panels_screen_MOFA.R -- Figure 2e, 2h-j and 5a-c panels.
# Produces: panels/Fig2e, Fig2h-j, Fig5a-c .png and .pdf (panel letters of the article).
# Inputs: 08c and 04 outputs, results/RA_ILD_Workspace.RData, Fig2h_BAL_Th_subsets.csv, Fig2ij_BAL_vs_PB.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
suppressPackageStartupMessages({
  library(tidyverse); library(pROC); library(patchwork); library(jsonlite)
})
source("scripts/00_cohort_config.R")
cfg <- set_cohort()
TB <- cfg$tables
panel_dir <- file.path(cfg$figs, "panels"); dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
pmap <- fromJSON("scripts/figures/panel_map.json")   # pixel sizes of the figure panels (shipped with the code)

rd <- function(f) { alt <- file.path(TB, "supp_tables", f); read.csv(if (file.exists(alt)) alt else file.path(TB, f), check.names = FALSE) }

# ---- value check log --------------------------------------------------------
chk <- data.frame()
log_val <- function(panel, item, printed, table_value, source) {
  ok <- identical(as.character(printed), as.character(table_value))
  chk <<- rbind(chk, data.frame(panel = panel, item = item, printed = printed,
                                table_value_formatted = table_value, source = source,
                                match = ok))
  if (!ok) stop(sprintf("%s %s: printed '%s' != table '%s'", panel, item, printed, table_value))
}
same <- function(a, b, what, tol = 1e-9) {
  if (abs(a - b) > tol) stop(sprintf("%s: recomputed %.10g != table %.10g", what, a, b))
}

# ---- styling (verbatim from the source scripts) ------------------------------
tn_A <- function(bs = 8) {
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
tn_I <- function(bs = 7.5) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(.9), color = "black"),
          axis.title = element_text(size = rel(.95), face = "bold"),
          plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0),
          legend.text = element_text(size = rel(.78)),
          legend.key.size = unit(2.5, "mm"),
          legend.background = element_rect(fill = alpha("white", .95), color = NA),
          panel.border = element_rect(fill = NA, color = "black", linewidth = .3),
          axis.line = element_blank(),
          strip.background = element_rect(fill = "grey95", color = "black", linewidth = .2),
          strip.text = element_text(size = rel(.85), face = "bold"),
          plot.margin = unit(c(1, 2, 1, 2), "mm"))
}
tn_L <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
          axis.line = element_blank(),
          plot.margin = ggplot2::margin(2, 3, 2, 3, "mm"))
}
tn_5 <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0,
                                    margin = ggplot2::margin(0, 0, 3, 0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 3), "mm"))
}
sv <- function(p, name, w, h) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p, width = w, height = h,
         units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p, width = w, height = h,
         units = "mm", dpi = 300, bg = "white")
  px <- dim(png::readPNG(file.path(panel_dir, paste0(name, ".png"))))[2:1]
  tgt <- as.integer(pmap[[name]]$px)
  
  if (!all(px == tgt)) stop(name, ": pixel size differs from the base panel")
}
CG <- c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")
CUD_ORANGE <- "#FF9900"; CUD_BLUE <- "#0041FF"
add_ep <- function(spec, sens) {
  ord <- order(spec, sens); spec <- spec[ord]; sens <- sens[ord]
  if (spec[1] != 0 || sens[1] != 0) { spec <- c(0, spec); sens <- c(0, sens) }
  if (tail(spec, 1) != 1 || tail(sens, 1) != 1) { spec <- c(spec, 1); sens <- c(sens, 1) }
  list(spec = spec, sens = sens)
}
fmt4 <- function(p) ifelse(p < 0.0001, "p<0.0001", sprintf("p=%.4f", p))
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))  # figure format = manuscript text
rf <- function(r) sprintf("ρ = %.2f", r)

# ---- raw data (workspace master_data) ----
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
md <- e$master_data; rm(e); invisible(gc())
md <- md[!md$Sample_ID %in% cfg$EXCL, ]
check_n(md$Sample_ID, cfg, "master_data")

nice2f <- function(v) {
  pre <- ifelse(grepl("^BALF_", v), "BAL ", ifelse(grepl("^PB_", v), "PB ", ifelse(grepl("^Serum_", v), "Serum ", "")))
  x <- gsub("^BALF_|^PB_|^Serum_", "", v)
  cyt <- c(IL6 = "IL-6", GCSF = "G-CSF", MCP3 = "MCP-3", GROa = "GRO-α", IL13 = "IL-13", IL22 = "IL-22", IL4 = "IL-4", IFNa2 = "IFN-α2", IL17A = "IL-17A",
           TGFa = "TGF-α", MIP1a = "MIP-1α", GMCSF = "GM-CSF", VEGFA = "VEGF-A", IL2 = "IL-2", IL1a = "IL-1α", IL1b = "IL-1β", IL8 = "IL-8", IL27 = "IL-27",
           IL10 = "IL-10", IL17F = "IL-17F", IL5 = "IL-5", IL7 = "IL-7", FGF2 = "FGF-2", IL12p70 = "IL-12p70", TNFb = "TNF-β", FLT3L = "FLT-3L", IL1RA = "IL-1RA",
           IFNg = "IFN-γ", TNFa = "TNF-α", MIP1b = "MIP-1β", IL12p40 = "IL-12p40", IL15 = "IL-15", IL18 = "IL-18", IP10 = "IP-10", MCP1 = "MCP-1", PDGFAA = "PDGF-AA")
  fcm <- c(`Th17.1_CCR6pos_CXCR3pos` = "Th17.1", Th2_CCR6neg_CXCR3neg = "Th2-enriched", Th1_CCR6neg_CXCR3pos = "Th1", Th17_CCR6pos_CXCR3neg = "Th17",
           Treg_CD127neg_CD25pos = "Treg", Neutrophil_Percent = "neutrophils (%)", Neutrophil_CD66bpos_CD206neg = "CD66b+CD206− neutrophils",
           Macrophage_CD14pos_CD16neg = "CD14+CD16− monocytes", Activated_Macrophage_CD86pos_CD14pos = "CD86+CD14+ monocytes")
  out <- ifelse(x %in% names(cyt), cyt[x], ifelse(x %in% names(fcm), fcm[x], x))
  out[v == "BALF_Macrophage_CD14pos_CD16neg"] <- "CD14+CD16− cells"; out[v == "BALF_Activated_Macrophage_CD86pos_CD14pos"] <- "CD86+CD14+ cells"
  paste0(pre, out) }
es <- rd("SD2_effect_sizes_RA_vs_Sarc.csv")
src_lab <- function(cat) case_when(
  cat == "Serum_Cytokine" ~ "Serum", cat == "BALF_Cytokine" ~ "BAL Cytokine",
  cat == "BALF_FCM" ~ "BAL Flow Cytometry", cat == "PB_FCM" ~ "PB Flow Cytometry", cat == "BAL_differential_count" ~ "BAL Differential Count",
  grepl("Deconv", cat) ~ "Deconv", TRUE ~ "Other")
ses_all <- es %>% filter(!is.na(p_adjusted), p_adjusted < 0.05) %>%
  mutate(L = nice2f(Variable), Src = src_lab(Category)) %>%
  arrange(cliff_delta)
stopifnot(nrow(ses_all) == 32)
ses_serum <- ses_all %>% filter(Src == "Serum") %>% tail(10)
srm <- sort(ses_all$cliff_delta[ses_all$Src == "Serum"], decreasing = TRUE)
stopifnot(srm[10] > srm[11])
ses <- rbind(ses_all %>% filter(Src != "Serum"), ses_serum) %>% arrange(cliff_delta)
stopifnot(!any(duplicated(ses$Variable)))

write.csv(ses[, c("Variable", "Category", "Src", "cliff_delta", "p_value", "p_adjusted")],
          file.path(panel_dir, "panels_screen_MOFA_Fig2e_variables.csv"),
          row.names = FALSE)
ses$Variable <- factor(ses$Variable, levels = ses$Variable)
labs_y <- setNames(ses$L, ses$Variable)
pal2f <- c(Serum = "#E74C3C", "BAL Cytokine" = "#3498DB",
           "BAL Flow Cytometry" = "#9B59B6", "PB Flow Cytometry" = "#27AE60", "BAL Differential Count" = "#7F7F7F",
           Deconv = "#E67E22", Other = "#95A5A6")
pal2f <- pal2f[names(pal2f) %in% ses$Src]
p2f <- ggplot(ses, aes(cliff_delta, Variable, fill = Src)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = pal2f) +
  scale_y_discrete(labels = labs_y) +
  { if (Sys.getenv("F2F_BREAKS", "101") == "101") scale_x_continuous(breaks = c(-1, 0, 1)) else NULL } +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(x = "Cliff's delta (RA vs Sarcoidosis)", y = "", fill = "") +
  tn_I(7) + guides(fill = guide_legend(ncol = 2)) + theme(legend.position = "bottom", legend.location = "plot", legend.key.size = unit(3, "mm"), legend.margin = margin(-2, 0, 0, 0),
                  axis.text.y = element_text(size = 5.5),
                  legend.text = element_text(size = 9))
log_val("Fig2e", "number of BH<0.05 variables (pool)", nrow(ses_all),
        sum(es$p_adjusted < 0.05, na.rm = TRUE), "SD2_effect_sizes_RA_vs_Sarc.csv")
log_val("Fig2e", "variables shown (non-serum + top10 serum)", nrow(ses),
        sum(es$p_adjusted < 0.05 & es$Category != "Serum_Cytokine", na.rm = TRUE) + 10,
        "SD2_effect_sizes_RA_vs_Sarc.csv")
sv(p2f, "Fig2e", 100, 90)

f2i <- rd("Fig2h_BAL_Th_subsets.csv")
tc <- c("BALF_Th1_CCR6neg_CXCR3pos", "BALF_Th2_CCR6neg_CXCR3neg",
        "BALF_Th17_CCR6pos_CXCR3neg", "BALF_Th17.1_CCR6pos_CXCR3pos",
        "BALF_Treg_CD127neg_CD25pos")
tn2_labels <- c("Th1", "Th2-enriched", "Th17", "Th17.1", "Treg")   # axis labels of the article
stopifnot(identical(f2i$Variable, tc))
td <- data.frame()
for (i in seq_along(tc)) {
  v <- as.numeric(md[[tc[i]]]); g <- md$Sample_Group; ok <- !is.na(v) & !is.na(g)
  td <- rbind(td, data.frame(Sub = tn2_labels[i], V = v[ok],
                             G = ifelse(g[ok] == "Control", "Sarcoidosis", "RA")))
}
td$Sub <- factor(td$Sub, levels = tn2_labels)
pv_td <- sapply(tn2_labels, function(s) {
  d <- td[td$Sub == s, ]
  suppressWarnings(wilcox.test(d$V[d$G == "RA"], d$V[d$G == "Sarcoidosis"], exact = TRUE)$p.value)
})
for (i in seq_along(tc)) {
  same(pv_td[[i]], f2i$p[i], paste("Fig2h p", tn2_labels[i]))
  same(median(td$V[td$Sub == tn2_labels[i] & td$G == "RA"]), f2i$median_RA[i], "Fig2h median RA")
  same(median(td$V[td$Sub == tn2_labels[i] & td$G == "Sarcoidosis"]), f2i$median_Sarc[i], "Fig2h median Sarc")
  stopifnot(sum(td$Sub == tn2_labels[i] & td$G == "RA") == f2i$n_RA[i],
            sum(td$Sub == tn2_labels[i] & td$G == "Sarcoidosis") == f2i$n_Sarc[i])
}
lab_of <- function(p) ifelse(p < 0.05, paste0("*", pf(p)), "NS")
pv_lab <- data.frame(
  Sub = factor(tn2_labels, levels = tn2_labels),
  Y = sapply(tn2_labels, function(s) max(td$V[td$Sub == s], na.rm = TRUE) * 1.12),
  Lab = lab_of(pv_td))
for (i in seq_along(tc))
  log_val("Fig2h", paste(gsub("\n", "", tn2_labels[i]), "p label"), pv_lab$Lab[i], lab_of(f2i$p[i]),
          sprintf("Fig2h_BAL_Th_subsets.csv p = %.3g", f2i$p[i]))
p2i <- ggplot(td, aes(Sub, V, fill = G)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.2) +
  geom_jitter(aes(group = G),
              position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.5, seed = 42),
              size = 0.6, alpha = 0.4, show.legend = FALSE) +
  scale_fill_manual(values = CG) +
  geom_text(data = pv_lab, aes(Sub, Y, label = Lab), inherit.aes = FALSE,
            size = 3.75, color = "black") +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = "", y = "% of CD4+", fill = "") +
  tn_A(8) + theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 7.5),
                  legend.key.size = unit(3, "mm"), legend.margin = margin(-2, 0, 7.2, 0))
set.seed(42)
sv(p2i, "Fig2h", 100, 75)

jk <- rd("Fig2ij_BAL_vs_PB.csv")
bvp <- function(bc, pc, nm) {
  ba <- as.numeric(md[[bc]]); pb <- as.numeric(md[[pc]])
  gg <- ifelse(md$Sample_Group == "Control", "Sarcoidosis", "RA")
  v <- !is.na(ba) & !is.na(pb)
  sp <- suppressWarnings(cor.test(ba[v], pb[v], method = "spearman", exact = TRUE))
  r <- jk[jk$Subset == nm, ]
  same(unname(sp$estimate), r$Spearman_rho, paste("rho", nm)); same(sp$p.value, r$p, paste("p", nm))
  stopifnot(sum(v) == r$N)
  pan <- if (nm == "Th17.1") "Fig2i" else "Fig2j"
  log_val(pan, "rho", sprintf("%.3f", sp$estimate), sprintf("%.3f", r$Spearman_rho), "Fig2ij_BAL_vs_PB.csv")
  log_val(pan, "p", fmt4(sp$p.value), fmt4(r$p), sprintf("Fig2ij_BAL_vs_PB.csv p = %.3g", r$p))
  df <- data.frame(X = pb[v], Y = ba[v], G = gg[v])
  ggplot(df, aes(X, Y, color = G)) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    scale_color_manual(values = CG) +
    annotate("text", x = if (nm == "Treg") Inf else -Inf, y = Inf, hjust = if (nm == "Treg") 1.05 else -0.05, vjust = 1.3, size = 3,
             label = paste0(rf(sp$estimate), "\n", pf(sp$p.value))) +
    labs(x = paste("PB", nm, "(%)"), y = paste("BAL", nm, "(%)"), color = "") +
    tn_I(8) + theme(legend.position = "bottom", legend.direction = "horizontal",
                    legend.text = element_text(size = 6))
}
sv(bvp("BALF_Th17.1_CCR6pos_CXCR3pos", "PB_Th17.1_CCR6pos_CXCR3pos", "Th17.1"), "Fig2i", 80, 70)
sv(bvp("BALF_Treg_CD127neg_CD25pos", "PB_Treg_CD127neg_CD25pos", "Treg"), "Fig2j", 80, 70)

view_labels <- c(Expression = "Expression", BALF_Cytokine = "BAL Cytokine",
                 Serum_Cytokine = "Serum Cytokine", BALF_FACS = "BAL FCM",
                 PB_FACS = "PB FCM", Microbiome = "Microbiome")
view_cols <- c(Expression = "#E67E22", BALF_Cytokine = "#3498DB", Serum_Cytokine = "#E74C3C",
               BALF_FACS = "#27AE60", PB_FACS = "#9B59B6", Microbiome = "#95A5A6")
r2t <- rd("diseaseMOFA2_R2_per_factor_view.csv")
met <- rd("diseaseMOFA2_metrics.csv"); mv <- setNames(met$Value, met$Metric)
stopifnot(mv[["Disease factor (retrained index)"]] == "1")   # disease factor is Factor 1
em <- new.env(); suppressWarnings(load(file.path(TB, "diseaseMOFA2_model.RData"), envir = em))
r2m <- em$r2$r2_per_factor[[1]]; r2tot <- em$r2$r2_total[[1]]; rm(em)
r2c <- as.matrix(r2t[r2t$Factor != "Total", -1]); rownames(r2c) <- r2t$Factor[r2t$Factor != "Total"]
stopifnot(identical(dimnames(r2m), dimnames(r2c)), max(abs(round(r2m, 2) - r2c)) < 1e-9,
          max(abs(round(r2tot, 2) - unlist(r2t[r2t$Factor == "Total", -1]))) < 1e-9)
stopifnot(sprintf("%.1f", r2m["Factor1", "Serum_Cytokine"]) == mv[["Disease factor R2 Serum_Cytokine (%)"]],
          sprintf("%.1f", r2m["Factor1", "BALF_Cytokine"]) == mv[["Disease factor R2 BALF_Cytokine (%)"]])
fnames <- setNames(paste("Factor", 1:5), paste0("Factor", 1:5))
r2l <- as.data.frame(r2m) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2")
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2m)))
r2l$View <- factor(r2l$View, levels = colnames(r2m))
r2l$FLabel <- factor(fnames[as.character(r2l$Factor)], levels = rev(fnames[rownames(r2m)]))
for (k in which(r2l$R2 > 0.5))
  log_val("Fig5a", paste(r2l$Factor[k], r2l$View[k]), sprintf("%.1f", r2l$R2[k]),
          sprintf("%.1f", r2m[as.character(r2l$Factor[k]), as.character(r2l$View[k])]),
          sprintf("diseaseMOFA2_model.RData r2 = %.4f (CSV 2dp = %.2f)",
                  r2m[as.character(r2l$Factor[k]), as.character(r2l$View[k])],
                  r2c[as.character(r2l$Factor[k]), as.character(r2l$View[k])]))
p_hm <- ggplot(r2l, aes(View, FLabel, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient(low = "white", high = "#2471A3", name = "R² (%)") +
  geom_text(aes(label = ifelse(R2 > 0.5, sprintf("%.1f", R2), ""),
                color = ifelse(R2 > 25, "white", "black")), size = 2.5, fontface = "bold",
            show.legend = FALSE) +
  scale_color_identity() + scale_x_discrete(labels = view_labels) +
  labs(x = "", y = "") +
  tn_5(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7,
                                             color = view_cols[levels(r2l$View)]),
                  legend.position = "left", legend.key.height = unit(10, "mm"))
total_r2 <- data.frame(View = factor(names(r2tot), levels = colnames(r2m)), Total = as.numeric(r2tot))
p_bar <- ggplot(total_r2, aes(Total, View, fill = View)) +
  geom_col(width = 0.65) + scale_fill_manual(values = view_cols) +
  scale_y_discrete(labels = view_labels) + labs(x = "Total R² (%)", y = "") +
  scale_x_continuous(breaks = c(0, 40, 80)) +
  tn_5(7) + theme(legend.position = "none", axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
sv(p_hm + p_bar + plot_layout(widths = c(4, 1.2)), "Fig5a", 130, 75)
# the 'Total R² (%)' title of the bar chart sits 60 px (1.5 character heights) higher in the article image: the title pixels are moved up
f5a <- file.path(panel_dir, "Fig5a.png"); im <- png::readPNG(f5a)
reg <- im[701:885, 1101:1535, 1:3, drop = FALSE]; ink <- which(apply(reg < 1, c(1, 2), any), arr.ind = TRUE)
ys <- 700 + range(ink[, 1]); xs <- 1100 + range(ink[, 2]); lab <- im[ys[1]:ys[2], xs[1]:xs[2], , drop = FALSE]
im[ys[1]:ys[2], xs[1]:xs[2], 1:3] <- 1; im[(ys[1]:ys[2]) - 60, xs[1]:xs[2], ] <- lab; png::writePNG(im, f5a)

fs <- rd("diseaseMOFA2_factor_scores.csv")
check_n(fs$Sample_ID, cfg, "factor scores")
pfac <- rd("diseaseMOFA2_factor_disease_p.csv")
pval <- wilcox.test(fs$Disease_factor_aligned[fs$Group == "RA"],
                    fs$Disease_factor_aligned[fs$Group == "Sarcoidosis"], exact = TRUE)$p.value
same(pval, pfac$p_disease_exact[pfac$Factor == 1], "Fig5b p")
plab <- pf(pval)
log_val("Fig5b", "p", plab,
        pf(pfac$p_disease_exact[pfac$Factor == 1]),
        sprintf("diseaseMOFA2_factor_disease_p.csv p = %.6f", pfac$p_disease_exact[1]))
df5b <- data.frame(Group = factor(fs$Group, levels = c("RA", "Sarcoidosis")),
                   Value = fs$Disease_factor_aligned)
set.seed(42)
p5b <- ggplot(df5b, aes(Group, Value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.6, color = "grey25") +
  annotate("text", x = 1.5, y = max(df5b$Value), label = plab, size = 2.8,
           vjust = -0.2) +
  scale_fill_manual(values = c(RA = "#C0392B", Sarcoidosis = "#2980B9")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = "", y = "Factor 1 value") + ggtitle("Factor 1") +
  tn_5(8) + theme(legend.position = "none")
sv(p5b, "Fig5b", 62, 75)

y <- ifelse(fs$Group == "RA", 1, 0)
roc_f1 <- roc(y, fs$Disease_factor_aligned, quiet = TRUE, direction = "auto")
same(as.numeric(auc(roc_f1)), as.numeric(mv[["Disease factor AUC"]]), "Fig5c F1 AUC", 5e-4)
rcm <- rd("diseaseMOFA2_ROC_coords.csv")
stopifnot(all.equal(sort(rcm$Sens[rcm$Model == "DiseaseFactor"]), sort(roc_f1$sensitivities)))
ssa <- rd("diseaseMOFA2_single_serum_AUC.csv")
serum <- md[match(fs$Sample_ID, md$Sample_ID), ]
best <- ssa$Serum[1]
curves <- best
nice <- function(cc) { cf <- gsub("^Serum_", "", cc); cf <- gsub("^GCSF$", "G-CSF", cf)
  gsub("^IL(\\d)", "IL-\\1", cf) }
ep <- add_ep(1 - roc_f1$specificities, roc_f1$sensitivities)
f1lab <- sprintf("Factor 1 (AUC = %.3f)", auc(roc_f1))
log_val("Fig5c", "Factor 1 AUC", sprintf("%.3f", auc(roc_f1)), mv[["Disease factor AUC"]],
        "diseaseMOFA2_metrics.csv")
roc_df <- data.frame(Spec = ep$spec, Sens = ep$sens, Model = f1lab)
for (cc in curves) {
  rr <- roc(y, as.numeric(serum[[cc]]), quiet = TRUE, direction = "auto")
  same(as.numeric(auc(rr)), ssa$AUC[ssa$Serum == cc], paste("Fig5c AUC", cc))
  log_val("Fig5c", paste(nice(cc), "AUC"), sprintf("%.3f", auc(rr)),
          sprintf("%.3f", ssa$AUC[ssa$Serum == cc]), "diseaseMOFA2_single_serum_AUC.csv")
  ep <- add_ep(1 - rr$specificities, rr$sensitivities)
  roc_df <- rbind(roc_df, data.frame(Spec = ep$spec, Sens = ep$sens,
                                     Model = sprintf("%s (AUC = %.3f)", paste("Serum", nice(cc)), auc(rr))))  # compartment named
}
log_val("Fig5c", "best single serum cytokine", nice(best), "IL-6",
               "diseaseMOFA2_metrics.csv (Serum_IL6)")
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))
cols5c <- c(CUD_BLUE, CUD_ORANGE, "#35A16B")[seq_len(nlevels(roc_df$Model))]
lty5c <- c("solid", "solid", "22")[seq_len(nlevels(roc_df$Model))]
p5c <- ggplot(roc_df, aes(Spec, Sens, color = Model, linetype = Model)) +
  geom_step(linewidth = 0.6, direction = "vh") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = cols5c) + scale_linetype_manual(values = lty5c) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "on") +
  labs(x = "1 – Specificity", y = "Sensitivity", color = "", linetype = "") +
  ggtitle("RA vs Sarcoidosis") +
  tn_5(8) + theme(legend.position = c(0.62, if (nlevels(roc_df$Model) > 2) 0.21 else 0.18),
                  legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
                  legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm"))
sv(p5c, "Fig5c", 80, 75)

# ---- value log ----------------------------------------------------------------
write.csv(chk, file.path(panel_dir, "panels_screen_MOFA_values.csv"), row.names = FALSE)

print(chk[, c("panel", "item", "printed", "table_value_formatted", "match")], row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(panel_dir, "panels_screen_MOFA_sessionInfo.txt"))
