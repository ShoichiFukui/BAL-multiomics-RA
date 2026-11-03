# panels_fig3a_3b_3h.R -- Figure 3a, 3b and 3h panels.
# Produces: panels/Fig3a, Fig3b, Fig3h .png and .pdf.
# Inputs: 06, 17 (GSVA GO OXPHOS score; Figure 3h) and 02 outputs, planA_data_RA24.csv (infection outcome).
# Figure 3a PERMANOVA (R2 = 0.091, p = 0.027; 31-sample VST, RA columns) is read from PERMANOVA_infection_RA24.csv of
#   16_PERMANOVA_infection_n31.R.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(forcats) })
R0 <- file.path(ART_TABLES, "n31"); PD <- file.path(ART_FIGS, "n31", "panels"); dir.create(PD, recursive = TRUE, showWarnings = FALSE)
sv <- function(p, name, w, h, px) { ggsave(file.path(PD, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(PD, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  d <- dim(png::readPNG(file.path(PD, paste0(name, ".png")))); stopifnot(d[2] == px[1], d[1] == px[2]) }
tn <- function(bs = 8) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"), axis.title = element_text(size = rel(0.95), face = "bold"),
        plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 3, 0)), legend.text = element_text(size = rel(0.8)),
        legend.key.size = unit(3, "mm"), legend.background = element_rect(fill = alpha("white", 0.95), color = NA),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33), axis.line = element_blank(),
        strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.3), strip.text = element_text(size = rel(0.85), face = "bold"), plot.margin = unit(c(2, 3, 2, 3), "mm"))
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))
nice <- function(v) {
  pre <- ifelse(grepl("^BALF_", v), "BAL ", ifelse(grepl("^PB_", v), "PB ", ifelse(grepl("^Serum_", v), "Serum ", "")))
  x <- gsub("^BALF_|^PB_|^Serum_|^GSVA_|^Pathway_|^BAL_", "", v)
  cyt <- c(IL6 = "IL-6", GCSF = "G-CSF", MCP3 = "MCP-3", GROa = "GRO-α", IL13 = "IL-13", IL22 = "IL-22", IL4 = "IL-4", IFNa2 = "IFN-α2", IL17A = "IL-17A", TGFa = "TGF-α", MIP1a = "MIP-1α", GMCSF = "GM-CSF", VEGFA = "VEGF-A", IL2 = "IL-2", IL10 = "IL-10")
  fcm <- c(`Th17.1_CCR6pos_CXCR3pos` = "Th17.1", Th2_CCR6neg_CXCR3neg = "Th2-enriched", Th1_CCR6neg_CXCR3pos = "Th1", Th17_CCR6pos_CXCR3neg = "Th17", Treg_CD127neg_CD25pos = "Treg",
           CD3 = "CD3+ T cells", CD22 = "CD22+ B cells", CD56 = "CD56+ NK cells", Th2_count = "Th2-enriched count", Th17.1_count = "Th17.1 count", Th1_count = "Th1 count", Th17_count = "Th17 count", Treg_count = "Treg count")
  out <- ifelse(x %in% names(cyt), cyt[x], ifelse(x %in% names(fcm), fcm[x], x))
  pre[grepl("^BAL_", v)] <- "BAL "; paste0(pre, out) }
# ---- Figure 3b ----
ip <- read.csv(file.path(R0, "common_panel/infection_panel.csv"))
ti <- ip %>% filter(P < 0.05) %>% arrange(Cliff) %>% mutate(L = nice(Variable), Category = recode(Category, "BAL cytokine" = "Cytokine", "BAL FCM" = "FCM", "BAL lymphocyte subsets (clinical)" = "FCM", "BAL CD4+ subset count" = "FCM"))
stopifnot(nrow(ti) == 12, all(ti$Category %in% c("Cytokine", "FCM")));   # 12 of the 114 panel markers at nominal p < 0.05
 ti$L <- fct_inorder(ti$L); print(ti[, c("Variable", "L", "Category", "Cliff", "P")], row.names = FALSE)
write.csv(ti[, c("Variable", "L", "Category", "Median_pos", "Median_neg", "N_pos", "N_neg", "AUC", "Cliff", "P", "BH")], file.path(PD, "Fig3b_12markers.csv"), row.names = FALSE)
sv(ggplot(ti, aes(Cliff, L, fill = Category)) + geom_col(width = 0.5) + scale_fill_brewer(palette = "Set2") + geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "Cliff's delta (Infection + vs −)", y = "", fill = "") + tn(7) +
     theme(legend.position = "bottom", legend.direction = "horizontal", legend.location = "plot", legend.margin = margin(-2, 0, 0, 0), legend.key.size = unit(3.5, "mm"), axis.text.y = element_text(size = 8.5), legend.text = element_text(size = 9),
           legend.background = element_rect(fill = alpha("white", 0.9), color = NA)), "Fig3b", 95, 95, c(1122, 1122))
# ---- Figure 3h ----
gv <- read.csv(file.path(R0, "gsva_GO/gsva_GO_scores_n31.csv")); ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))   # 17_GSVA_GO_scores.R
ox <- subset(gv, Group == "RA"); stopifnot(nrow(ox) == 24)   # GSVA score of GO BP oxidative phosphorylation, RA rows in file order
dfo <- data.frame(V = ox$OXPHOS, IG = factor(ifelse(ra$respiratory_infection[match(ox$Sample_ID, ra$Sample_ID)] == 1, "Infection (+)", "Infection (−)"), levels = c("Infection (−)", "Infection (+)")))
wto <- wilcox.test(dfo$V[dfo$IG == "Infection (+)"], dfo$V[dfo$IG == "Infection (−)"], exact = TRUE); stopifnot(abs(wto$p.value - 0.00942) < 5e-5)
CI <- c("Infection (−)" = "#3498DB", "Infection (+)" = "#E74C3C"); set.seed(42)
sv(ggplot(dfo, aes(IG, V, fill = IG)) + geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25, alpha = 0.7) +
     geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65, show.legend = FALSE) + scale_fill_manual(values = CI) +
     scale_color_manual(values = c("Infection (−)" = "#1A5276", "Infection (+)" = "#922B21")) + labs(x = "", y = "OXPHOS GSVA score") +
     annotate("text", x = 1.5, y = max(dfo$V) * 1.08, label = pf(wto$p.value), size = 2.5) + tn(8) + theme(legend.position = "none", axis.text.x = element_text(size = 9.5)), "Fig3h", 80, 70, c(944, 826))

pca <- read.csv(file.path(R0, "transcriptome_PCA_scores.csv"))   # 02_transcriptome_RA_vs_sarcoidosis.R (all 31 samples, as Figure 2A)
pa3 <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); inf <- data.frame(Sample_ID = pa3$Sample_ID, Infection = pa3$respiratory_infection)
pd <- merge(pca, inf, by = "Sample_ID"); stopifnot(nrow(pd) == 24); pd <- pd[order(public_label(pd$Sample_ID)), ]   # drawing order = de-identified labels
pd$InfSt <- factor(ifelse(pd$Infection == 1, "Infection (+)", "Infection (−)"), levels = c("Infection (−)", "Infection (+)"))
v <- read.csv(file.path(R0, "transcriptome_PCA_variance.csv")); print(head(v, 2))
pc <- as.numeric(unlist(v[1:2, ncol(v)])); if (max(pc) <= 1) pc <- pc * 100
cen <- aggregate(cbind(PC1m = PC1, PC2m = PC2) ~ InfSt, pd, mean)
write.csv(pd[, c("Sample_ID", "Group", "PC1", "PC2", "Infection")], file.path(PD, "Fig3a_PCA_RA24_n31coords.csv"), row.names = FALSE)
pm <- read.csv(file.path(R0, "PERMANOVA_infection_RA24.csv")); pm <- pm[pm$Matrix == "VST_31", ]   # 16_PERMANOVA_infection_n31.R
perm_lab <- sprintf("PERMANOVA\nR² = %.3f, p = %.3f", pm$R2, pm$P); stopifnot(perm_lab == "PERMANOVA\nR² = 0.091, p = 0.027")
sv(ggplot(pd, aes(PC1, PC2, color = InfSt)) + stat_ellipse(level = 0.95, linewidth = 0.5) + geom_point(size = 2, alpha = 0.8) +
     geom_point(data = cen, aes(PC1m, PC2m), shape = 4, size = 2.4, stroke = 1.0, show.legend = FALSE) +
     scale_color_manual(values = c("Infection (−)" = "#3498DB", "Infection (+)" = "#E74C3C")) +
     annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.5, label = perm_lab) +
     labs(x = sprintf("PC1 (%.1f%%)", pc[1]), y = sprintf("PC2 (%.1f%%)", pc[2]), color = "") + tn(8) + theme(legend.position = "bottom"), "Fig3a", 80, 70, c(944, 826))
