#!/usr/bin/env Rscript
# =================================================================
# Comprehensive Figure Generation Script
# Hosogaya BAL — Nature Communications manuscript
# =================================================================
# Generates all figure panels from pre-computed analysis results.
# Output panel names match PPT v21 assembly.
#
# Prerequisites (run in order):
#   01_BayesPrism_Deconvolution.R -> output_v3_BayesPrism/
#   02_PostDeconvolution.R        -> results/RA_ILD_Workspace.RData
#   analysis_modules/Enhanced_Analysis.R -> results/tables/*.csv
#   analysis_modules/CT_Multiomics.R     -> results/CT_Multiomics_Results.RData
#   analysis_modules/Integration.R       -> results/Integration_Results.RData
#   03_Analysis.R (builds 6-view MOFA2)  -> results/MOFA2_6views_Results.RData
#   05_MOFA2_RA_only.R                   -> results/MOFA2_RA_only_Results.RData
#
# Manual panels not generated here:
#   Fig1a (study design), Fig4a (CT GMM images x3)
# =================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(patchwork)
  library(readxl); library(pROC); library(effsize); library(randomForest)
  library(grid); library(vegan); library(Seurat); library(ggrepel)
  library(glmnet)
})

BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
set.seed(42)

panel_dir <- "results/panels"
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

EXCL <- character(0)

# ─── Load data ────────────────────────────────────────────────────
cat("Loading data...\n")
load("results/RA_ILD_Workspace.RData")
load("results/CT_Multiomics_Results.RData")
load("results/MOFA2_6views_Results.RData")
load("results/MOFA2_RA_only_Results.RData")

master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if (exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)

bp_cell <- read.csv("output_v3_BayesPrism/celltype_proportions_cellFraction.csv", row.names = 1)
bp_cell <- bp_cell[!rownames(bp_cell) %in% EXCL, ]
bp_theta <- read.csv("output_v3_BayesPrism/celltype_proportions_BayesPrism.csv", row.names = 1)
bp_theta <- bp_theta[!rownames(bp_theta) %in% EXCL, ]
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet = "BALF_analysis")
colnames(fcm)[1] <- "Sample_ID"

deg <- read.csv("results/tables/DEG_ShrunkLFC.csv")
ggo <- read.csv("results/tables/GSEA_GO_BP.csv")
gkk <- read.csv("results/tables/GSEA_KEGG.csv")
es <- read.csv("results/tables/Effect_sizes_RA_vs_Control.csv")
inf_comp <- read.csv("results/tables/Infection_AllOmics_Comparison.csv")
bio_roc <- read.csv("results/tables/Infection_SingleBiomarker_ROC.csv")
deg_inf <- read.csv("results/tables/DEG_Infection_pos_vs_neg.csv")
ct_data <- read.csv("results/tables/master_data_with_CT_all.csv")
ct_unique <- ct_data[!duplicated(ct_data$Sample_ID), ]
s3_table <- read.csv("results/tables/S3_Healthy_Serum_Comparison.csv")

ra_md <- master_data[master_data$Sample_Group == "RA" &
                       !is.na(master_data$respiratory_infection), ]
ra_md$IG <- ifelse(ra_md$respiratory_infection == 1,
                   "Infection (+)", "Infection (-)")
ra_ext2 <- master_ext[master_ext$Sample_Group == "RA" &
                        master_ext$Sample_ID %in% common_samples, ]

samples <- common_samples
grp_all <- recode(as.character(master_ext$Subgroup), "Control" = "Sarcoidosis",
                  "RA_nonILD" = "RA-nonILD", "RA_ILD" = "RA-ILD")

# ─── Theme and helpers ────────────────────────────────────────────
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
      strip.background = element_rect(fill = "grey95", color = "black",
                                      linewidth = 0.3),
      strip.text = element_text(size = rel(0.85), face = "bold"),
      plot.margin = unit(c(2, 3, 2, 3), "mm")
    )
}

sv <- function(p, name, w = 80, h = 70) {
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p,
         width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name, ".png")), p,
         width = w, height = h, units = "mm", dpi = 300, bg = "white")
  cat(sprintf("  Saved: %s (%dx%d mm)\n", name, w, h))
}

fmt_p <- function(p) {
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.2f", p)
}

CG <- c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")
CS <- c("Sarcoidosis" = "#2980B9", "RA-nonILD" = "#27AE60", "RA-ILD" = "#C0392B")
CI <- c("Infection (-)" = "#3498DB", "Infection (+)" = "#E74C3C")
CP <- c("Stable" = "#3498DB", "Progressor" = "#E74C3C")
CC <- c(Macrophage = "#E74C3C", T_cell = "#3498DB", NK = "#2ECC71",
        B_cell = "#9B59B6", Plasma = "#E67E22", Neutrophil = "#F1C40F",
        DC = "#8B4513", Epithelial = "#E91E63", Mast = "#95A5A6",
        Other = "#BDC3C7")
col_ds <- c("GSE145926" = "#E74C3C", "GSE193782" = "#3498DB",
            "GSE184735" = "#27AE60")
view_cols <- c("Expression" = "#E67E22", "BALF_Cytokine" = "#3498DB",
               "Serum_Cytokine" = "#E74C3C", "BALF_FACS" = "#27AE60",
               "PB_FACS" = "#9B59B6", "Microbiome" = "#7F8C8D")
view_labels <- c("Expression" = "Expression", "BALF_Cytokine" = "BAL Cytokine",
                 "Serum_Cytokine" = "Serum Cytokine", "BALF_FACS" = "BAL FCM",
                 "PB_FACS" = "PB FCM", "Microbiome" = "Microbiome")

sct <- function(x, y, xl, yl, col = "#C0392B", bs = 8) {
  v <- !is.na(x) & !is.na(y)
  sp <- cor.test(x[v], y[v], method = "spearman", exact = TRUE)
  df <- data.frame(X = x[v], Y = y[v])
  ggplot(df, aes(X, Y)) +
    geom_point(size = 1.8, color = col, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.4,
                fill = "grey90", formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60",
               linewidth = 0.2) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3.5,
             label = paste0("\u03C1=", sprintf("%.3f", sp$estimate), "\n",
                            ifelse(sp$p.value < 0.0001, "p<0.0001",
                                   sprintf("p=%.4f", sp$p.value)))) +
    labs(x = xl, y = yl) + tn(bs) + theme(legend.position = "none")
}

bxi <- function(col, yl, rd, bs = 8) {
  v <- as.numeric(rd[[col]]); ig <- rd$IG
  ok <- !is.na(v) & !is.na(ig)
  df <- data.frame(V = v[ok], IG = ig[ok])
  wt <- wilcox.test(df$V[df$IG == "Infection (+)"],
                    df$V[df$IG == "Infection (-)"], exact=TRUE)
  ymx <- max(df$V, na.rm = TRUE); ymn <- min(df$V, na.rm = TRUE)
  ysp <- (ymx - ymn) * 0.12
  ggplot(df, aes(IG, V, fill = IG)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25, alpha = 0.7) +
    geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65,
                show.legend = FALSE) +
    scale_fill_manual(values = CI) +
    scale_color_manual(values = c("Infection (-)" = "#1A5276",
                                  "Infection (+)" = "#922B21")) +
    labs(x = "", y = yl) +
    annotate("text", x = 1.5, y = ymx + ysp * 1.2,
             label = ifelse(wt$p.value < 0.0001, "p<0.0001",
                            sprintf("p=%.4f", wt$p.value)), size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    tn(bs) + theme(legend.position = "none",
                   axis.text.x = element_text(size = 9.5))
}

bvp <- function(bc, pc, nm) {
  ba <- as.numeric(master_data[[bc]]); pb <- as.numeric(master_data[[pc]])
  gg <- ifelse(master_data$Sample_Group == "Control", "Sarcoidosis", "RA")
  v <- !is.na(ba) & !is.na(pb)
  sp <- cor.test(ba[v], pb[v], method = "spearman", exact = TRUE)
  df <- data.frame(X = pb[v], Y = ba[v], G = gg[v])
  ggplot(df, aes(X, Y, color = G)) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
    scale_color_manual(values = CG) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3,
             label = paste0("\u03C1=", sprintf("%.3f", sp$estimate), "\n",
                            ifelse(sp$p.value < 0.0001, "p<0.0001",
                                   sprintf("p=%.4f", sp$p.value)))) +
    labs(x = paste("PB", nm, "(%)"), y = paste("BAL", nm, "(%)"), color = "") +
    tn(8) + theme(legend.position = "bottom", legend.direction = "horizontal",
                  legend.text = element_text(size = 6))
}


# =================================================================
#  FIGURE 1: BAL fluid cell-type deconvolution
# =================================================================
cat("\n========== FIGURE 1 ==========\n")

# PCA setup
pca <- prcomp(t(expr_matrix), scale. = TRUE)
ve <- summary(pca)$importance[2, 1:5] * 100
pd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 SID = rownames(pca$x)) %>%
  mutate(Group = ifelse(grepl("^(KYC|Sarcoidosis)", SID), "Sarcoidosis", "RA"),
         Sub = recode(as.character(meta$Subgroup[match(SID, meta$Sample_ID)]),
                      "Control" = "Sarcoidosis", "RA_nonILD" = "RA-nonILD",
                      "RA_ILD" = "RA-ILD"),
         InfSt = recode(as.character(meta$Infection_Group[match(SID, meta$Sample_ID)]),
                        "Infection_Negative" = "Infection (-)",
                        "Infection_Positive" = "Infection (+)"))
xl <- sprintf("PC1 (%.1f%%)", ve[1])
yl <- sprintf("PC2 (%.1f%%)", ve[2])

# Fig1b — UMAP by dataset (paper Fig 1b)
ref <- readRDS("output_v3_BayesPrism/BAL_reference_author_annotated.rds")
uc <- Embeddings(ref, "umap"); ms <- ref@meta.data
ds_orig <- data.frame(U1 = uc[, 1], U2 = uc[, 2], DS = ms$dataset,
                      CT = ms$celltype)

# Dotplot data (before shuffling)
mks <- intersect(c("CD68","MARCO","FABP4","FCN1","CD3D","CD3E","NKG7","GNLY",
                    "CD79A","MS4A1","JCHAIN","MZB1","FCGR3B","FCER1A","EPCAM",
                    "KRT18"), rownames(ref))
ed <- GetAssayData(ref, layer = "data"); dd <- data.frame()
for (ct in unique(ds_orig$CT)) {
  cl <- which(ds_orig$CT == ct)
  if (length(cl) < 10) next
  for (mk in mks) {
    vv <- ed[mk, cl]
    dd <- rbind(dd, data.frame(CT = ct, M = mk,
                               Pct = sum(vv > 0) / length(vv) * 100,
                               Avg = mean(vv)))
  }
}

set.seed(42); ds <- ds_orig[sample(nrow(ds_orig)), ]

sv(ggplot(ds, aes(U1, U2, color = DS)) +
     geom_point(size = 0.06, alpha = 0.2) +
     scale_color_manual(values = col_ds,
                        labels = c("GSE145926" = "COVID-19",
                                   "GSE193782" = "Healthy",
                                   "GSE184735" = "Sarcoidosis")) +
     guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
     labs(x = "UMAP 1", y = "UMAP 2", color = "") +
     tn(8) + theme(legend.position = "bottom"),
   "Fig1c_UMAP_dataset", 85, 75)

# Fig1c — UMAP by cell type (paper Fig 1c)
sv(ggplot(ds, aes(U1, U2, color = CT)) +
     geom_point(size = 0.06, alpha = 0.2) +
     scale_color_manual(values = CC) +
     guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                 ncol = 5)) +
     labs(x = "UMAP 1", y = "UMAP 2", color = "") +
     tn(8) + theme(legend.position = "bottom",
                   legend.text = element_text(size = 6),
                   legend.key.size = unit(3, "mm")),
   "Fig1d_UMAP_celltype", 90, 85)

# Fig1d — Dotplot (paper Fig 1d)
dd$M <- factor(dd$M, levels = rev(mks))
sv(ggplot(dd, aes(CT, M, size = Pct, color = Avg)) +
     geom_point() +
     scale_size_continuous(range = c(0.5, 3.5), name = "%Exp") +
     scale_color_gradient(low = "grey90", high = "#C0392B", name = "Avg") +
     labs(x = "", y = "") +
     tn(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                   axis.text.y = element_text(size = 6, face = "italic")),
   "Fig1e_dotplot", 100, 80)

# Fig1e — Reference composition (paper Fig 1e)
cm <- ds %>% group_by(DS, CT) %>% summarise(N = n(), .groups = "drop") %>%
  group_by(DS) %>% mutate(P = N / sum(N))
cm$DS <- recode(cm$DS, "GSE145926" = "COVID-19", "GSE193782" = "Healthy",
                "GSE184735" = "Sarcoidosis")
sv(ggplot(cm, aes(DS, P, fill = CT)) +
     geom_col(width = 0.6) + scale_fill_manual(values = CC) +
     labs(x = "", y = "Proportion", fill = "") +
     tn(8) + theme(legend.position = "right",
                   legend.text = element_text(size = 7),
                   legend.key.size = unit(3.5, "mm")),
   "Fig1f_composition", 110, 75)

rm(ref, ed); gc()

# Fig1f — BayesPrism stacked bar (paper Fig 1f)
bp_l <- bp_cell %>% rownames_to_column("S") %>%
  mutate(G = ifelse(grepl("^(KYC|Sarcoidosis)", S), "Sarcoidosis", "RA"),
         S_d = ifelse(grepl("^(KYC|Sarcoidosis)", S),
                      paste0("Sarcoidosis", gsub("^(KYC0*|Sarcoidosis)", "", S)),
                      paste0("RA", gsub("^(KY0*|RA)", "", S)))) %>%
  pivot_longer(-c(S, G, S_d), names_to = "CT", values_to = "Fr") %>%
  mutate(S_d = fct_reorder(S_d, as.numeric(factor(G))))
sv(ggplot(bp_l, aes(S_d, Fr, fill = CT)) +
     geom_col(width = 0.85) + scale_fill_manual(values = CC) +
     labs(x = "", y = "Proportion", fill = "") +
     tn(7) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4),
                   legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = element_text(size = 6)),
   "Fig3a_stacked_bar", 140, 70)

# Fig1g-i — FCM validation (paper Fig 1g-i)
cf <- intersect(rownames(bp_cell), fcm$Sample_ID)
fi_idx <- match(cf, fcm$Sample_ID)
lc <- intersect(c("T_cell", "NK", "B_cell", "Plasma"), colnames(bp_cell))

mkf <- function(bv, fv, lb) {
  v <- !is.na(fv) & !is.na(bv)
  sp <- cor.test(bv[v], fv[v], method = "spearman", exact = TRUE)
  df <- data.frame(F = fv[v], B = bv[v],
                   G = ifelse(grepl("^(KYC|Sarcoidosis)", cf[v]), "Sarcoidosis", "RA"))
  ggplot(df, aes(F, B, color = G)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.3,
                linetype = "dashed", formula = y ~ x, inherit.aes = FALSE,
                aes(x = F, y = B)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                color = "grey60") +
    scale_color_manual(values = CG) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3,
             label = paste0("\u03C1=", sprintf("%.3f", sp$estimate), "\n",
                            ifelse(sp$p.value < 0.0001, "p<0.0001",
                                   sprintf("p=%.4f", sp$p.value)))) +
    labs(x = paste("FCM", lb, "(%)"), y = paste("Deconv.", lb, "(%)"),
         color = "") +
    tn(7) + theme(legend.position = "bottom", legend.direction = "horizontal",
                  legend.text = element_text(size = 6))
}

sv(mkf(bp_cell[cf, "Macrophage"] * 100,
       as.numeric(fcm$BALF_Macrophage_Percent[fi_idx]), "Mac"),
   "Fig3b_FCM_Mac")
sv(mkf(rowSums(bp_cell[cf, lc]) * 100,
       as.numeric(fcm$BALF_Lymphocyte_Percent[fi_idx]), "Lym"),
   "Fig3c_FCM_Lym")
sv(mkf(bp_cell[cf, "Neutrophil"] * 100,
       as.numeric(fcm$BALF_Neutrophil_Percent[fi_idx]), "Neu"),
   "Fig3d_FCM_Neu")


# =================================================================
#  FIGURE 2: RA vs Sarcoidosis profiling + classification
# =================================================================
cat("\n========== FIGURE 2 ==========\n")

# Fig2a — PCA by subgroup (paper Fig 2a)
sv(ggplot(pd %>% filter(!is.na(Sub)), aes(PC1, PC2, color = Sub)) +
     geom_point(size = 2, alpha = 0.8) +
     scale_color_manual(values = CS) +
     labs(x = xl, y = yl, color = "") +
     tn(8) + theme(legend.position = "bottom"),
   "Fig1b_PCA_subgroup")

# Fig2b — Volcano (paper Fig 2b)
deg$sig <- case_when(
  deg$padj < 0.05 & deg$log2FoldChange > 0.585 ~
    "Up (padj<0.05, |LFC|>0.585)",
  deg$padj < 0.05 & deg$log2FoldChange < -0.585 ~
    "Down (padj<0.05, |LFC|>0.585)",
  deg$padj < 0.05 ~ "Sig (padj<0.05)",
  TRUE ~ "NS (padj>=0.05)")
deg$lab <- ifelse(!is.na(deg$padj) & deg$padj < 0.05, deg$Gene_Symbol, NA)
sv(ggplot(deg, aes(log2FoldChange, -log10(padj), color = sig)) +
     geom_point(data = deg %>% filter(grepl("^NS", sig)),
                size = 1.0, alpha = 0.05) +
     geom_point(data = deg %>% filter(!grepl("^NS", sig)),
                size = 1.8, alpha = 0.9) +
     geom_text_repel(data = deg %>% filter(!is.na(lab)), aes(label = lab),
                     size = 2.2, fontface = "italic", color = "black",
                     show.legend = FALSE, max.overlaps = 20,
                     segment.size = 0.3, segment.color = "grey40",
                     min.segment.length = 0, box.padding = 0.5,
                     point.padding = 0.3, force = 2, force_pull = 0.5) +
     scale_color_manual(values = c(
       "Up (padj<0.05, |LFC|>0.585)" = "#C0392B",
       "Down (padj<0.05, |LFC|>0.585)" = "#2980B9",
       "Sig (padj<0.05)" = "#E67E22",
       "NS (padj>=0.05)" = "grey70"),
       guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                color = "grey50", linewidth = 0.2) +
     geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed",
                color = "grey50", linewidth = 0.2) +
     xlim(-2, 2) + labs(x = "Shrunken log2FC", y = "-log10(padj)", color = "") +
     tn(8) + theme(legend.position = "right",
                   legend.text = element_text(size = 7.5),
                   legend.key.size = unit(3.5, "mm")),
   "Fig2a_volcano", 115, 90)

# Fig2c — GSEA GO (paper Fig 2c)
tgo <- head(ggo[order(ggo$p.adjust), ], 10) %>%
  mutate(D = str_trunc(Description, 44))
tgo <- tgo[order(tgo$NES), ]; tgo$D <- fct_inorder(tgo$D)
sv(ggplot(tgo, aes(NES, D, fill = NES > 0)) +
     geom_col(width = 0.6) +
     scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9"),
                       labels = c("TRUE" = "Up in RA", "FALSE" = "Down in RA")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "NES", y = "", fill = "") +
     tn(7) + theme(axis.text.y = element_text(size = 8),
                   legend.position = "bottom",
                   legend.text = element_text(size = 7),
                   plot.margin = unit(c(1, 2, 4, 2), "mm")),
   "Fig2b_GSEA_GO", 110, 85)

# Fig2d — GSEA KEGG (paper Fig 2d)
tkk <- head(gkk[order(gkk$p.adjust), ], 10) %>%
  mutate(D = str_trunc(Description, 44))
tkk <- tkk[order(tkk$NES), ]; tkk$D <- fct_inorder(tkk$D)
sv(ggplot(tkk, aes(NES, D, fill = NES > 0)) +
     geom_col(width = 0.6, show.legend = FALSE) +
     scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "NES", y = "") + tn(7) + theme(axis.text.y = element_text(size = 8)),
   "Fig2c_GSEA_KEGG", 110, 75)

# Fig2e — GSVA heatmap (paper Fig 2e)
gm <- gsva_scores[, common_samples]
so <- common_samples[order(meta$Sample_Group[match(common_samples, meta$Sample_ID)])]
glg <- as.data.frame(gm) %>% rownames_to_column("PW") %>%
  pivot_longer(-PW, names_to = "S", values_to = "Sc")
glg$S_d <- ifelse(grepl("^(KYC|Sarcoidosis)", glg$S),
                  paste0("Sarcoidosis", gsub("^(KYC0*|Sarcoidosis)", "", glg$S)),
                  paste0("RA", gsub("^(KY0*|RA)", "", glg$S)))
so_d <- ifelse(grepl("^(KYC|Sarcoidosis)", so),
               paste0("Sarcoidosis", gsub("^(KYC0*|Sarcoidosis)", "", so)),
               paste0("RA", gsub("^(KY0*|RA)", "", so)))
glg$S_d <- factor(glg$S_d, levels = so_d)
sv(ggplot(glg, aes(S_d, PW, fill = Sc)) +
     geom_tile() +
     scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                          midpoint = 0, name = "Score") +
     annotate("rect", xmin = which(grepl("^Sarcoidosis", so_d))[1] - 0.5,
              xmax = tail(which(grepl("^Sarcoidosis", so_d)), 1) + 0.5,
              ymin = nrow(gm) + 0.55, ymax = nrow(gm) + 1.20,
              fill = "#2980B9", alpha = 0.8) +
     annotate("rect", xmin = which(grepl("^RA", so_d))[1] - 0.5,
              xmax = tail(which(grepl("^RA", so_d)), 1) + 0.5,
              ymin = nrow(gm) + 0.55, ymax = nrow(gm) + 1.20,
              fill = "#C0392B", alpha = 0.8) +
     annotate("text", x = median(which(grepl("^Sarcoidosis", so_d))),
              y = nrow(gm) + 0.8, label = "Sarcoidosis", size = 2.5,
              color = "white", fontface = "bold") +
     annotate("text", x = median(which(grepl("^RA", so_d))),
              y = nrow(gm) + 0.8, label = "RA", size = 2.5,
              color = "white", fontface = "bold") +
     coord_cartesian(clip = "off") +
     labs(x = "", y = "") +
     tn(6) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 3),
                   axis.text.y = element_text(size = 6),
                   plot.margin = unit(c(10, 5, 1, 4), "mm")),
   "Fig2d_GSVA_heatmap", 120, 85)

# Fig2f — Effect sizes (paper Fig 2f)
ses_all <- es %>% filter(!is.na(p_adjusted), p_adjusted < 0.05) %>%
  mutate(L = gsub("^BALF_|^Serum_|^PB_", "", Variable),
         Src = case_when(
           grepl("Serum", Category) ~ "Serum",
           grepl("BALF_Cytokine", Category) ~ "BAL Cytokine",
           grepl("FCM", Category) ~ "Flow Cytometry",
           grepl("Deconv", Category) ~ "Deconv",
           grepl("GSVA", Category) ~ "GSVA", TRUE ~ "Other")) %>%
  arrange(cliff_delta)
ses_serum <- ses_all %>% filter(Src == "Serum") %>% tail(10)
ses_other <- ses_all %>% filter(Src != "Serum")
ses <- rbind(ses_other, ses_serum) %>% arrange(cliff_delta)
ses$L <- fct_inorder(ses$L)
sv(ggplot(ses, aes(cliff_delta, L, fill = Src)) +
     geom_col(width = 0.5) +
     scale_fill_manual(values = c(Serum = "#E74C3C", "BAL Cytokine" = "#3498DB",
                                  "Flow Cytometry" = "#9B59B6",
                                  Deconv = "#E67E22", GSVA = "#27AE60",
                                  Other = "#95A5A6")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "Cliff's delta (RA vs Sarcoidosis)", y = "", fill = "") +
     tn(7) + theme(legend.position = "right",
                   axis.text.y = element_text(size = 5.5),
                   legend.text = element_text(size = 6)),
   "Fig2f_effect_sizes", 100, 90)

# Fig2g — Plasma cells (paper Fig 2g) — uses bp_cell (cell fraction)
ga_pl <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(bp_cell)), "Sarcoidosis", "RA")
dp <- data.frame(P = bp_cell[, "Plasma"] * 100, G = ga_pl)
wtp <- wilcox.test(dp$P[dp$G == "RA"], dp$P[dp$G == "Sarcoidosis"], exact=TRUE)
sv(ggplot(dp, aes(G, P + 1e-6, fill = G)) +
     geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.2) +
     geom_jitter(width = 0.07, size = 1.2, alpha = 0.5, show.legend = FALSE) +
     scale_y_log10() + scale_fill_manual(values = CG) +
     labs(x = "", y = "Plasma cell fraction (%, log10)") +
     annotate("text", x = 1.5, y = max(dp$P) * 3,
              label = ifelse(wtp$p.value < 0.0001, "p<0.0001",
                             sprintf("p=%.4f", wtp$p.value)), size = 2.5) +
     tn(8) + theme(legend.position = "none"),
   "Fig2g_plasma")

# Fig2h — Ig genes (paper Fig 2h)
ig <- deg[grepl("^IGHG|^IGHM$|^IGKC$|^IGHA[12]$|^JCHAIN$", deg$Gene_Symbol) &
            !is.na(deg$padj), ] %>%
  arrange(desc(abs(log2FoldChange))) %>% head(7)
ig$Gene_Symbol <- fct_reorder(ig$Gene_Symbol, ig$log2FoldChange)
ig$sg <- ifelse(ig$padj < 0.05, "padj < 0.05", "padj >= 0.05")
sv(ggplot(ig, aes(log2FoldChange, Gene_Symbol, fill = sg)) +
     geom_col(width = 0.5) +
     scale_fill_manual(values = c("padj < 0.05" = "#C0392B",
                                  "padj >= 0.05" = "#95A5A6")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "log2FC (RA vs Sarcoidosis)", y = "", fill = "") +
     tn(8) + theme(legend.position = "bottom",
                   legend.text = element_text(size = 6),
                   axis.text.y = element_text(face = "italic")),
   "Fig2h_Ig_genes")

# Fig2i — T cell subsets (paper Fig 2i)
tc <- c("BALF_Th17_CCR6pos_CXCR3neg", "BALF_Th17.1_CCR6pos_CXCR3pos",
        "BALF_Th1_CCR6neg_CXCR3pos", "BALF_Th2_CCR6neg_CXCR3neg",
        "BALF_Treg_CD127neg_CD25pos")
tn2_labels <- c("Th1", "Th2", "Th17", "Th17.1", "Treg")
td <- data.frame()
for (i in seq_along(tc)) {
  if (tc[i] %in% colnames(master_data)) {
    v <- as.numeric(master_data[[tc[i]]])
    g <- master_data$Sample_Group
    ok <- !is.na(v) & !is.na(g)
    td <- rbind(td, data.frame(
      Sub = tn2_labels[i], V = v[ok],
      G = ifelse(g[ok] == "Control", "Sarcoidosis", "RA")))
  }
}
td$Sub <- factor(td$Sub, levels = tn2_labels)
pv_td <- sapply(tn2_labels, function(s) {
  d <- td[td$Sub == s, ]
  tryCatch(wilcox.test(d$V[d$G == "RA"], d$V[d$G == "Sarcoidosis"], exact=TRUE)$p.value,
           error = function(e) NA)
})
pv_lab <- data.frame(
  Sub = factor(tn2_labels, levels = tn2_labels),
  Y = sapply(tn2_labels, function(s) max(td$V[td$Sub == s], na.rm = TRUE) * 1.12),
  Lab = ifelse(pv_td < 0.05,
               ifelse(pv_td < 0.0001, "*p<0.0001",
                      sprintf("*p=%.4f", pv_td)), "NS"))
sv(ggplot(td, aes(Sub, V, fill = G)) +
     geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.2) +
     geom_jitter(aes(group = G),
                 position = position_jitterdodge(jitter.width = 0.08,
                                                 dodge.width = 0.5),
                 size = 0.6, alpha = 0.4, show.legend = FALSE) +
     scale_fill_manual(values = CG) +
     geom_text(data = pv_lab, aes(Sub, Y, label = Lab), inherit.aes = FALSE,
               size = 2, color = "grey30") +
     labs(x = "", y = "% of CD4+", fill = "") +
     tn(8) + theme(legend.position = "right"),
   "Fig3e_Th_subsets", 100, 75)

# Fig2j-k — BAL vs PB (paper Fig 2j-k)
sv(bvp("BALF_Th17.1_CCR6pos_CXCR3pos", "PB_Th17.1_CCR6pos_CXCR3pos", "Th17.1"),
   "Fig3f_BAL_PB_Th171")
sv(bvp("BALF_Treg_CD127neg_CD25pos", "PB_Treg_CD127neg_CD25pos", "Treg"),
   "Fig3g_BAL_PB_Treg")

# Fig2l — Nested LOOCV ROC (paper Fig 2l)
lo_s <- intersect(intersect(rownames(expr_mat), rownames(cyto_mat)),
                  rownames(facs_mat))
lo_d <- data.frame(expr_mat[lo_s, 1:min(50, ncol(expr_mat))],
                   cyto_mat[lo_s, ], facs_mat[lo_s, ],
                   Group = meta$Sample_Group[match(lo_s, meta$Sample_ID)])
lo_d <- lo_d[!is.na(lo_d$Group), ]; lo_d$Group <- factor(lo_d$Group)
for (cc in setdiff(colnames(lo_d), "Group"))
  lo_d[[cc]][is.na(lo_d[[cc]]) | is.infinite(lo_d[[cc]])] <-
    median(lo_d[[cc]], na.rm = TRUE)
fvv <- apply(lo_d[, setdiff(colnames(lo_d), "Group"), drop = FALSE], 2,
             var, na.rm = TRUE)
lo_d <- lo_d[, c(names(fvv[!is.na(fvv) & fvv > 0]), "Group")]
nl <- nrow(lo_d)
set.seed(42)
np <- matrix(NA, nl, 2, dimnames = list(NULL, levels(lo_d$Group)))
for (i in 1:nl) {
  tr <- lo_d[-i, ]
  ft <- setdiff(colnames(tr), "Group")
  pv <- apply(tr[, ft], 2, function(x)
    tryCatch(wilcox.test(x ~ tr$Group, exact=TRUE)$p.value, error = function(e) 1))
  tp <- names(sort(pv))[1:min(50, length(pv))]
  rf <- randomForest(x = tr[, tp], y = tr$Group, ntree = 500)
  np[i, ] <- stats::predict(rf, lo_d[i, tp], type = "prob")
}
roc1 <- roc(lo_d$Group, np[, "RA"], quiet = TRUE)
ci1 <- ci.auc(roc1, method = "bootstrap", boot.n = 2000, quiet = TRUE)
rd1 <- data.frame(S = roc1$sensitivities, Sp = 1 - roc1$specificities)
rd1 <- rd1[order(rd1$Sp, rd1$S), ]
sv(ggplot(rd1, aes(Sp, S)) +
     geom_step(color = "#C0392B", linewidth = 0.8, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                 color = "grey60") +
     annotate("text", x = 0.55, y = 0.1, size = 2.8, hjust = 0,
              label = sprintf("AUC=%.3f\n95%%CI:%.3f-%.3f",
                              as.numeric(auc(roc1)), ci1[1], ci1[3])) +
     labs(x = "1-Specificity", y = "Sensitivity") +
     coord_equal() + tn(8),
   "Fig2e_ROC_LOOCV")


# =================================================================
#  FIGURE 3: Infection prediction
# =================================================================
cat("\n========== FIGURE 3 ==========\n")

# Fig3a — PCA infection (paper Fig 3a)
pd_inf <- pd %>% filter(!is.na(InfSt))
cen <- pd_inf %>% group_by(InfSt) %>%
  summarise(PC1m = mean(PC1), PC2m = mean(PC2), .groups = "drop")
expr_inf_m <- t(expr_matrix[, pd_inf$SID])
dist_inf <- dist(expr_inf_m)
set.seed(42)
suppressMessages(perm_inf <- adonis2(dist_inf ~ pd_inf$InfSt, permutations = 999))
sv(ggplot(pd_inf, aes(PC1, PC2, color = InfSt)) +
     stat_ellipse(level = 0.95, linewidth = 0.5) +
     geom_point(size = 2, alpha = 0.8) +
     geom_point(data = cen, aes(PC1m, PC2m), shape = 4, size = 2.4,
                stroke = 1.0, show.legend = FALSE) +
     scale_color_manual(values = c("Infection (-)" = "#3498DB",
                                   "Infection (+)" = "#E74C3C")) +
     annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.5,
              label = sprintf("PERMANOVA\nR2=%.3f, p=%.4f",
                              perm_inf$R2[1], perm_inf$`Pr(>F)`[1])) +
     labs(x = xl, y = yl, color = "") +
     tn(8) + theme(legend.position = "bottom"),
   "Fig1g_PCA_infection")

# Fig3b — Cliff's delta infection (paper Fig 3b)
ti <- inf_comp %>% filter(P_value < 0.05) %>%
  mutate(L = gsub("BALF_|PB_|GSVA_|Deconv_|_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14.*",
                  "", Variable)) %>%
  arrange(Cliff_delta)
ti$L <- fct_inorder(ti$L)
ti$Category <- recode(ti$Category, "BALF_Cytokine" = "Cytokine",
                      "GSVA_Pathway" = "GSVA")
sv(ggplot(ti, aes(Cliff_delta, L, fill = Category)) +
     geom_col(width = 0.5) + scale_fill_brewer(palette = "Set2") +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "Cliff's delta (Infection + vs \u2212)", y = "", fill = "") +
     tn(7) + theme(legend.position = c(0.82, 0.15),
                   legend.key.size = unit(3, "mm"),
                   axis.text.y = element_text(size = 6),
                   legend.text = element_text(size = 7),
                   legend.background = element_rect(fill = alpha("white", 0.9),
                                                    color = NA)),
   "Fig4b_Cliff_delta", 95, 95)

# Fig3c-e — Infection boxplots (paper Fig 3c-e)
sv(bxi("BALF_Th17.1_CCR6pos_CXCR3pos", "BALF Th17.1 (%)", ra_md),
   "Fig4d_Th171_infection")
sv(bxi("BALF_Treg_CD127neg_CD25pos", "BALF Treg (%)", ra_md),
   "Fig4e_Treg_infection")

ra_md$Ratio_T171 <- as.numeric(ra_md[["BALF_Th17.1_CCR6pos_CXCR3pos"]]) /
  (as.numeric(ra_md[["PB_Th17.1_CCR6pos_CXCR3pos"]]) + 0.01)
dfr <- data.frame(V = ra_md$Ratio_T171, IG = ra_md$IG)
dfr <- dfr[!is.na(dfr$V) & !is.na(dfr$IG), ]
wtr2 <- wilcox.test(dfr$V[dfr$IG == "Infection (+)"],
                    dfr$V[dfr$IG == "Infection (-)"], exact=TRUE)
sv(ggplot(dfr, aes(IG, V, fill = IG)) +
     geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25,
                  alpha = 0.7) +
     geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65,
                 show.legend = FALSE) +
     scale_fill_manual(values = CI) +
     scale_color_manual(values = c("Infection (-)" = "#1A5276",
                                   "Infection (+)" = "#922B21")) +
     scale_y_log10() + labs(x = "", y = "Th17.1 BAL/PB (log10)") +
     annotate("text", x = 1.5, y = max(dfr$V, na.rm = TRUE) * 2,
              label = ifelse(wtr2$p.value < 0.0001, "p<0.0001",
                             sprintf("p=%.4f", wtr2$p.value)), size = 2.5) +
     tn(8) + theme(legend.position = "none",
                   axis.text.x = element_text(size = 9.5)),
   "Fig4f_ratio_infection")

# Fig3f — FCM infection ROC (paper Fig 3f)
y_inf <- ra_md$respiratory_infection
roc_biomarkers <- list(
  list(col = "BALF_Th17.1_CCR6pos_CXCR3pos", label = "BALF Th17.1",
       color = "#E74C3C"),
  list(col = "BALF_Treg_CD127neg_CD25pos", label = "BALF Treg",
       color = "#3498DB"),
  list(col = "ratio", label = "Th17.1 BAL/PB", color = "#27AE60")
)
roc_df <- data.frame()
for (bm in roc_biomarkers) {
  if (bm$col == "ratio") {
    vals <- as.numeric(ra_md[["BALF_Th17.1_CCR6pos_CXCR3pos"]]) /
      (as.numeric(ra_md[["PB_Th17.1_CCR6pos_CXCR3pos"]]) + 0.01)
  } else {
    vals <- as.numeric(ra_md[[bm$col]])
  }
  ok <- !is.na(vals) & !is.na(y_inf)
  ro <- roc(y_inf[ok], vals[ok], quiet = TRUE, direction = "auto")
  x <- 1 - ro$specificities; y <- ro$sensitivities
  ord <- order(x, y); x <- x[ord]; y <- y[ord]
  if (x[1] != 0 || y[1] != 0) { x <- c(0, x); y <- c(0, y) }
  if (tail(x, 1) != 1 || tail(y, 1) != 1) { x <- c(x, 1); y <- c(y, 1) }
  roc_df <- rbind(roc_df, data.frame(
    Spec = x, Sens = y,
    Model = sprintf("%s (%.3f)", bm$label, as.numeric(auc(ro)))))
  cat(sprintf("  %s: AUC = %.3f\n", bm$label, as.numeric(auc(ro))))
}
roc_df$Model <- factor(roc_df$Model, levels = unique(roc_df$Model))
sv(ggplot(roc_df, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.6, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                 color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#E74C3C", "#3498DB", "#27AE60")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     coord_equal() +
     tn(8) + theme(legend.position = c(0.62, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig3f_infection_FCM_ROC", 80, 70)

# Fig3g — Surfactant genes (paper Fig 3g)
sf <- c("SFTPB", "ETV5", "LPCAT1", "PGC", "SFTA3", "CLDN18",
        "SFTPD", "SFTPA1", "SFTPA2")
sfd <- deg_inf[deg_inf$Gene_Symbol %in% sf, ] %>% arrange(log2FoldChange)
sfd$Gene_Symbol <- fct_inorder(sfd$Gene_Symbol)
sfd$sg <- ifelse(!is.na(sfd$padj) & sfd$padj < 0.05, "padj<0.05",
                 ifelse(!is.na(sfd$padj) & sfd$padj < 0.1, "padj<0.1", "NS"))
sv(ggplot(sfd, aes(log2FoldChange, Gene_Symbol, fill = sg)) +
     geom_col(width = 0.5) +
     scale_fill_manual(values = c("padj<0.05" = "#C0392B",
                                  "padj<0.1" = "#E67E22", "NS" = "#95A5A6")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     scale_x_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 0.5)) +
     labs(x = "log2FC (Infection + vs \u2212)", y = "", fill = "") +
     tn(8) + theme(legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = element_text(size = 6),
                   legend.key.size = unit(2.5, "mm"),
                   axis.text.y = element_text(face = "italic")),
   "Fig4c_surfactant", 90, 90)

# Fig3h — OXPHOS infection (paper Fig 3h)
ra_ids_inf <- intersect(ra_md$Sample_ID, colnames(gsva_scores))
dfo <- data.frame(V = gsva_scores["OXPHOS", ra_ids_inf],
                  IG = ra_md$IG[match(ra_ids_inf, ra_md$Sample_ID)])
dfo <- dfo[!is.na(dfo$V) & !is.na(dfo$IG), ]
wto <- wilcox.test(dfo$V[dfo$IG == "Infection (+)"],
                   dfo$V[dfo$IG == "Infection (-)"], exact=TRUE)
sv(ggplot(dfo, aes(IG, V, fill = IG)) +
     geom_boxplot(outlier.shape = NA, width = 0.4, linewidth = 0.25,
                  alpha = 0.7) +
     geom_jitter(aes(color = IG), width = 0.07, size = 1.3, alpha = 0.65,
                 show.legend = FALSE) +
     scale_fill_manual(values = CI) +
     scale_color_manual(values = c("Infection (-)" = "#1A5276",
                                   "Infection (+)" = "#922B21")) +
     labs(x = "", y = "OXPHOS score") +
     annotate("text", x = 1.5, y = max(dfo$V) * 1.08,
              label = ifelse(wto$p.value < 0.0001, "p<0.0001",
                             sprintf("p=%.4f", wto$p.value)), size = 2.5) +
     tn(8) + theme(legend.position = "none",
                   axis.text.x = element_text(size = 9.5)),
   "Fig4g_OXPHOS_infection")


# =================================================================
#  FIGURE 4: CT-quantified lung disease progression
# =================================================================
cat("\n========== FIGURE 4 ==========\n")
# Fig4a: manual CT images (not generated)

# Fig4b — ILD vs GGO (paper Fig 4b)
sv(sct(as.numeric(ra_ext2$ILD_Score),
       as.numeric(ra_ext2$CT_Delta_GGO_pct),
       "ILD score", "Delta GGO (%)"),
   "Fig5b_ILD_vs_GGO")

# Fig4c — PB Th17.1 vs DeltaHL (paper Fig 4c)
sv(sct(as.numeric(ra_ext2$`PB_Th17.1_CCR6pos_CXCR3pos`),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "PB Th17.1 (%)", "Delta Healthy Lung (%)"),
   "Fig5c_PBTh171_vs_DeltaHL")

# Fig4d — PB Th17 vs DeltaHL (paper Fig 4d)
sv(sct(as.numeric(ra_ext2$`PB_Th17_CCR6pos_CXCR3neg`),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "PB Th17 (%)", "Delta Healthy Lung (%)"),
   "Fig5d_PBTh17_vs_DeltaHL")

# Fig4e — Neutrophil vs DeltaHL (paper Fig 4e)
sv(sct(as.numeric(ra_ext2$BALF_Neutrophil_Percent),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "BALF Neutrophil (%)", "Delta Healthy Lung (%)"),
   "Fig5e_Neutrophil_vs_DeltaHL")

# Fig4f — CT progression omics (paper Fig 4f)
ct_ra <- ct_unique[ct_unique$Sample_Group == "RA" &
                     !ct_unique$Sample_ID %in% EXCL, ]
dhl <- as.numeric(ct_ra$CT_Delta_Healthy_Lung_pct)
prog_st <- ifelse(dhl < -10, "Progressor", "Stable")
ra_ct <- data.frame(
  Sample_ID = ct_ra$Sample_ID,
  Prog = factor(prog_st, levels = c("Stable", "Progressor")),
  PB_Th17 = as.numeric(ra_md$PB_Th17_CCR6pos_CXCR3neg[
    match(ct_ra$Sample_ID, ra_md$Sample_ID)]),
  Serum_IL27 = as.numeric(ra_md$Serum_IL27[
    match(ct_ra$Sample_ID, ra_md$Sample_ID)]),
  Serum_IL8 = as.numeric(ra_md$Serum_IL8[
    match(ct_ra$Sample_ID, ra_md$Sample_ID)]),
  Serum_MDC = as.numeric(ra_md$Serum_MDC[
    match(ct_ra$Sample_ID, ra_md$Sample_ID)])
)

bm_info <- list(
  list(col = "PB_Th17", label = "PB Th17 (%)"),
  list(col = "Serum_IL27", label = "Serum IL-27 (pg/mL)"),
  list(col = "Serum_IL8", label = "Serum IL-8 (pg/mL)"),
  list(col = "Serum_MDC", label = "Serum MDC (pg/mL)")
)
plots_4f <- list()
for (bi in seq_along(bm_info)) {
  info <- bm_info[[bi]]
  vals <- as.numeric(ra_ct[[info$col]])
  ok <- !is.na(vals) & !is.na(ra_ct$Prog)
  df <- data.frame(Value = vals[ok], Prog = ra_ct$Prog[ok])
  wt <- wilcox.test(df$Value[df$Prog == "Progressor"],
                    df$Value[df$Prog == "Stable"], exact=TRUE)
  ymx <- max(df$Value, na.rm = TRUE)
  ysp <- (ymx - min(df$Value, na.rm = TRUE)) * 0.12
  p <- ggplot(df, aes(Prog, Value, fill = Prog)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25,
                 alpha = 0.7) +
    geom_jitter(aes(color = Prog), width = 0.08, size = 1.3, alpha = 0.6,
                show.legend = FALSE) +
    scale_fill_manual(values = CP) +
    scale_color_manual(values = c("Stable" = "#1A5276",
                                  "Progressor" = "#922B21")) +
    labs(x = "", y = info$label) +
    annotate("text", x = 1.5, y = ymx + ysp * 1.5,
             label = sprintf("p=%.3f", wt$p.value), size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    tn(8) + theme(legend.position = "none",
                  axis.text.x = element_text(size = 8))
  plots_4f[[bi]] <- p
  cat(sprintf("  %s: p=%.3f\n", info$label, wt$p.value))
}
sv(wrap_plots(plots_4f, nrow = 1),
   "Fig4f_CT_progression_omics", 210, 65)

# Fig4g — Sensitivity heatmap (paper Fig 4g)
thresholds <- c(-5, -10, -15, -20)
bm_sens <- c("PB_Th17_CCR6pos_CXCR3neg", "PB_Th17.1_CCR6pos_CXCR3pos",
             "BALF_Neutrophil_Percent", "Serum_IL27", "Serum_IL8", "Serum_MDC")
bm_lab <- c("PB Th17", "PB Th17.1", "BAL Neutrophil", "Serum IL-27",
            "Serum IL-8", "Serum MDC")
sens_mat <- matrix(NA, nrow = length(bm_sens), ncol = length(thresholds))
rownames(sens_mat) <- bm_lab
colnames(sens_mat) <- paste0("<", thresholds, "%")
pval_mat_s <- sens_mat

ra_ids_ct <- ct_ra$Sample_ID[ct_ra$Sample_ID %in% ra_md$Sample_ID]
dhl_ct <- as.numeric(ct_ra$CT_Delta_Healthy_Lung_pct[
  match(ra_ids_ct, ct_ra$Sample_ID)])

for (bi in seq_along(bm_sens)) {
  vals <- as.numeric(ra_md[[bm_sens[bi]]][match(ra_ids_ct, ra_md$Sample_ID)])
  for (ti in seq_along(thresholds)) {
    prog_bi <- ifelse(dhl_ct < thresholds[ti], 1, 0)
    ok <- !is.na(vals) & !is.na(prog_bi)
    if (sum(ok) >= 10 && length(unique(prog_bi[ok])) == 2) {
      wt <- wilcox.test(vals[ok & prog_bi == 1], vals[ok & prog_bi == 0], exact=TRUE)
      pval_mat_s[bi, ti] <- wt$p.value
      sens_mat[bi, ti] <- -log10(wt$p.value)
    }
  }
}

sens_long <- as.data.frame(sens_mat) %>% rownames_to_column("Biomarker") %>%
  pivot_longer(-Biomarker, names_to = "Threshold", values_to = "NegLogP")
pval_long_s <- as.data.frame(pval_mat_s) %>% rownames_to_column("Biomarker") %>%
  pivot_longer(-Biomarker, names_to = "Threshold", values_to = "P")
sens_long$P <- pval_long_s$P
sens_long$Sig <- ifelse(!is.na(sens_long$P) & sens_long$P < 0.05, "*", "")
sens_long$Biomarker <- factor(sens_long$Biomarker, levels = rev(bm_lab))
sens_long$Threshold <- factor(sens_long$Threshold,
                              levels = paste0("<", thresholds, "%"))

sv(ggplot(sens_long, aes(Threshold, Biomarker, fill = NegLogP)) +
     geom_tile(color = "white", linewidth = 0.5) +
     geom_text(aes(label = Sig), size = 4, color = "black", vjust = 0.8) +
     scale_fill_gradient(low = "white", high = "#C0392B",
                         name = "-log10(p)", na.value = "grey90") +
     geom_vline(xintercept = 1.5, linetype = "dashed", color = "grey70",
                linewidth = 0.3) +
     labs(x = "\u0394 Healthy Lung threshold", y = "") +
     tn(8) + theme(axis.text.y = element_text(size = 7)),
   "Fig4g_sensitivity_heatmap", 100, 75)

# Fig4h — Progression ROC (paper Fig 4h)
ra_ids_mofa <- rownames(ra_factors)
dhl_mofa <- as.numeric(ct_unique$CT_Delta_Healthy_Lung_pct[
  match(ra_ids_mofa, ct_unique$Sample_ID)])
prog_mofa <- ifelse(dhl_mofa < -10, 1, 0)

pb_th17_vals <- as.numeric(master_data[["PB_Th17_CCR6pos_CXCR3neg"]][
  match(ra_ids_mofa, master_data$Sample_ID)])
ok_th17 <- !is.na(prog_mofa) & !is.na(pb_th17_vals)
roc_th17 <- roc(prog_mofa[ok_th17], pb_th17_vals[ok_th17],
                quiet = TRUE, direction = "auto")

il27_vals <- as.numeric(master_data[["Serum_IL27"]][
  match(ra_ids_mofa, master_data$Sample_ID)])
ok_il27 <- !is.na(prog_mofa) & !is.na(il27_vals)
roc_il27 <- roc(prog_mofa[ok_il27], il27_vals[ok_il27],
                quiet = TRUE, direction = "auto")

il8_vals <- as.numeric(master_data[["Serum_IL8"]][
  match(ra_ids_mofa, master_data$Sample_ID)])
ok_il8 <- !is.na(prog_mofa) & !is.na(il8_vals)
roc_il8 <- roc(prog_mofa[ok_il8], il8_vals[ok_il8],
               quiet = TRUE, direction = "auto")

add_ep <- function(spec, sens) {
  ord <- order(spec, sens); spec <- spec[ord]; sens <- sens[ord]
  if (spec[1] != 0 || sens[1] != 0) { spec <- c(0, spec); sens <- c(0, sens) }
  if (tail(spec, 1) != 1 || tail(sens, 1) != 1) {
    spec <- c(spec, 1); sens <- c(sens, 1) }
  list(spec = spec, sens = sens)
}

ep1 <- add_ep(1 - roc_th17$specificities, roc_th17$sensitivities)
ep2 <- add_ep(1 - roc_il27$specificities, roc_il27$sensitivities)
ep3 <- add_ep(1 - roc_il8$specificities, roc_il8$sensitivities)

roc_prog <- rbind(
  data.frame(Sens = ep1$sens, Spec = ep1$spec,
             Model = sprintf("PB Th17 (%.3f)", auc(roc_th17))),
  data.frame(Sens = ep2$sens, Spec = ep2$spec,
             Model = sprintf("IL-27 (%.3f)", auc(roc_il27))),
  data.frame(Sens = ep3$sens, Spec = ep3$spec,
             Model = sprintf("IL-8 (%.3f)", auc(roc_il8)))
)
roc_prog$Model <- factor(roc_prog$Model, levels = unique(roc_prog$Model))
cat(sprintf("  PB Th17 AUC: %.3f\n", auc(roc_th17)))
cat(sprintf("  IL-27 AUC: %.3f\n", auc(roc_il27)))
cat(sprintf("  IL-8 AUC: %.3f\n", auc(roc_il8)))

sv(ggplot(roc_prog, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.6, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                 color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#C0392B", "#E67E22", "#3498DB")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     coord_equal() +
     tn(8) + theme(legend.position = c(0.62, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig4h_ROC_progression", 80, 70)


# =================================================================
#  FIGURE 5: MOFA2 integration
# =================================================================
cat("\n========== FIGURE 5 ==========\n")

group_m <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors)), "Sarcoidosis", "RA")
inf_m <- master_data$respiratory_infection[match(rownames(factors),
                                                  master_data$Sample_ID)]

fnames <- c(Factor1 = "Factor 1 (IDD)", Factor2 = "Factor 2 (IIV)",
            Factor3 = "Factor 3", Factor4 = "Factor 4", Factor5 = "Factor 5")

# Fig5a — Enhanced R2 heatmap (paper Fig 5a)
r2m <- r2$r2_per_factor[[1]]
r2l <- as.data.frame(r2m) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2")
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2m)))
r2l$View <- factor(r2l$View, levels = colnames(r2m))
r2l$FLabel <- fnames[as.character(r2l$Factor)]

p_hm <- ggplot(r2l, aes(View, FLabel, fill = R2)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient(low = "white", high = "#2471A3", name = "R\u00B2 (%)") +
  geom_text(aes(label = ifelse(R2 > 0.5, sprintf("%.1f", R2), ""),
                color = ifelse(R2 > 25, "white", "black")),
            size = 2.5, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_x_discrete(labels = view_labels) +
  labs(x = "", y = "") +
  tn(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 7,
                                            color = view_cols[levels(r2l$View)]),
                legend.position = "left", legend.key.height = unit(10, "mm"))

total_r2 <- data.frame(
  View = factor(names(r2$r2_total[[1]]), levels = colnames(r2m)),
  Total = as.numeric(r2$r2_total[[1]]))
p_bar <- ggplot(total_r2, aes(Total, View, fill = View)) +
  geom_col(width = 0.65) + scale_fill_manual(values = view_cols) +
  scale_y_discrete(labels = view_labels) +
  labs(x = "Total R\u00B2 (%)", y = "") +
  tn(7) + theme(legend.position = "none", axis.text.y = element_blank(),
                axis.ticks.y = element_blank())
sv(p_hm + p_bar + plot_layout(widths = c(4, 1.2)),
   "Fig6a_R2_enhanced", 130, 75)

# Fig5b — Infection ROC (paper Fig 5b)
ra_mask <- group_m == "RA"
ra_inf <- !is.na(inf_m) & ra_mask
y_inf_m <- inf_m[ra_inf]
roc_iiv <- roc(y_inf_m, factors[ra_inf, 2], quiet = TRUE, direction = "auto")
th171_vals_m <- as.numeric(master_data$BALF_Th17.1_CCR6pos_CXCR3pos[
  match(rownames(factors)[ra_inf], master_data$Sample_ID)])
roc_th171_m <- roc(y_inf_m, th171_vals_m, quiet = TRUE, direction = "auto")

roc_df_i <- rbind(
  data.frame(Sens = roc_iiv$sensitivities, Spec = 1 - roc_iiv$specificities,
             Model = sprintf("IIV (AUC = %.3f)", auc(roc_iiv))),
  data.frame(Sens = roc_th171_m$sensitivities,
             Spec = 1 - roc_th171_m$specificities,
             Model = sprintf("BAL Th17.1 (AUC = %.3f)", auc(roc_th171_m)))
)
sv(ggplot(roc_df_i, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.6, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                 color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#E74C3C", "#9B59B6")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     ggtitle("Infection prediction (RA)") +
     tn(8) + theme(legend.position = c(0.62, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig6f_ROC_infection", 80, 75)

# Fig5c — Disease factors (paper Fig 5c)
sf6 <- fd$Factor[fd$p_disease < 0.1]
fb6 <- data.frame()
for (f in sf6) {
  pv <- fd$p_disease[f]
  lab <- ifelse(f <= 2,
                paste0(fnames[paste0("Factor", f)], "\n(p=",
                       ifelse(pv < 0.0001, "<0.0001",
                              sprintf("%.4f", pv)), ")"),
                sprintf("Factor %d\n(p=%s)", f,
                        ifelse(pv < 0.0001, "<0.0001",
                               sprintf("%.4f", pv))))
  fb6 <- rbind(fb6, data.frame(Factor = lab, Value = factors[, f],
                                Group = group_m))
}
sv(ggplot(fb6, aes(Group, Value, fill = Group)) +
     geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25,
                  alpha = 0.7) +
     geom_jitter(aes(color = Group), width = 0.08, size = 1.3, alpha = 0.65,
                 show.legend = FALSE) +
     scale_fill_manual(values = CG) +
     scale_color_manual(values = c("RA" = "#7B241C",
                                   "Sarcoidosis" = "#1A5276")) +
     facet_wrap(~Factor, scales = "free_y") +
     labs(x = "", y = "Factor value") +
     tn(8) + theme(legend.position = "none"),
   "Fig6c_disease_factor", 95, 75)

# Fig5d — Disease ROC (paper Fig 5d)
y_disease <- ifelse(group_m == "RA", 1, 0)
roc_idd <- roc(y_disease, factors[, 1], quiet = TRUE, direction = "auto")
serum_data <- master_data[match(rownames(factors), master_data$Sample_ID), ]
best_cyto <- NULL; best_auc <- 0
for (cc in grep("^Serum_", colnames(serum_data), value = TRUE)) {
  vals <- as.numeric(serum_data[[cc]])
  if (sum(!is.na(vals)) >= 30) {
    r <- tryCatch(roc(y_disease, vals, quiet = TRUE, direction = "auto"),
                  error = function(e) NULL)
    if (!is.null(r) && as.numeric(auc(r)) > best_auc) {
      best_auc <- as.numeric(auc(r)); best_cyto <- cc
    }
  }
}
roc_best <- roc(y_disease, as.numeric(serum_data[[best_cyto]]),
                quiet = TRUE, direction = "auto")
cyto_fmt <- gsub("^Serum_", "", best_cyto)
cyto_fmt <- gsub("^GCSF$", "G-CSF", cyto_fmt)
cyto_fmt <- gsub("^IL(\\d)", "IL-\\1", cyto_fmt)

roc_df_d <- rbind(
  data.frame(Sens = roc_idd$sensitivities, Spec = 1 - roc_idd$specificities,
             Model = sprintf("IDD (AUC = %.3f)", auc(roc_idd))),
  data.frame(Sens = roc_best$sensitivities, Spec = 1 - roc_best$specificities,
             Model = sprintf("%s (AUC = %.3f)", cyto_fmt, auc(roc_best)))
)
sv(ggplot(roc_df_d, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.6, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                 color = "grey60", linewidth = 0.3) +
     annotate("segment", x = 0, xend = 1, y = 1, yend = 1,
              linetype = "dashed", color = "#27AE60", linewidth = 0.4) +
     annotate("text", x = 0.98, y = 0.05, size = 2.3, hjust = 1,
              color = "#27AE60",
              label = "Nested LOOCV (AUC = 0.962)", fontface = "italic") +
     scale_color_manual(values = c("#C0392B", "#E67E22")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     ggtitle("RA vs Sarcoidosis") +
     tn(8) + theme(legend.position = c(0.62, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig6e_ROC_disease", 80, 75)

# Fig5e — RA-only MOFA2 R2 heatmap (paper Fig 5e)
r2m_ra <- ra_r2$r2_per_factor[[1]]
r2l_ra <- as.data.frame(r2m_ra) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2")
r2l_ra$View <- recode(r2l_ra$View,
                      "BALF_FACS" = "BAL FCM", "PB_FACS" = "PB FCM",
                      "BALF_Cytokine" = "BAL Cytokine",
                      "Serum_Cytokine" = "Serum Cytokine",
                      "Expression" = "Expression")
r2l_ra$Factor <- factor(r2l_ra$Factor, levels = rev(rownames(r2m_ra)))
r2l_ra$View <- factor(r2l_ra$View,
                      levels = c("Expression", "BAL Cytokine",
                                 "Serum Cytokine", "BAL FCM", "PB FCM"))
sv(ggplot(r2l_ra, aes(View, Factor, fill = R2)) +
     geom_tile(color = "white", linewidth = 0.5) +
     scale_fill_gradient(low = "white", high = "#2171B5",
                         name = "R\u00B2 (%)") +
     geom_text(aes(label = ifelse(R2 >= 1, sprintf("%.1f", R2), "")),
               size = 2.5, color = "white", fontface = "bold") +
     geom_text(aes(label = ifelse(R2 > 0 & R2 < 1, sprintf("%.1f", R2), "")),
               size = 2.5, color = "grey30") +
     labs(x = "", y = "") +
     tn(8) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                   legend.position = "right",
                   legend.key.height = unit(12, "mm"),
                   legend.key.width = unit(3, "mm"),
                   legend.title = element_text(size = 7),
                   legend.text = element_text(size = 6)),
   "Fig5e_RA_MOFA2_R2", 100, 65)

# Fig5f — RA Factor 1 by infection (paper Fig 5f)
ra_ids_m <- rownames(ra_factors)
inf_ra <- master_data$respiratory_infection[match(ra_ids_m,
                                                   master_data$Sample_ID)]
ok_inf <- !is.na(inf_ra)
df_5f <- data.frame(
  Value = ra_factors[ok_inf, 1],
  IG = factor(ifelse(inf_ra[ok_inf] == 1, "Infection (+)", "Infection (-)"),
              levels = c("Infection (-)", "Infection (+)"))
)
wt_5f <- wilcox.test(df_5f$Value[df_5f$IG == "Infection (+)"],
                     df_5f$Value[df_5f$IG == "Infection (-)"], exact=TRUE)
cat(sprintf("  RA Factor 1 infection: p=%.4f\n", wt_5f$p.value))
ymx <- max(df_5f$Value); ysp <- (ymx - min(df_5f$Value)) * 0.12
sv(ggplot(df_5f, aes(IG, Value, fill = IG)) +
     geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25,
                  alpha = 0.7) +
     geom_jitter(aes(color = IG), width = 0.08, size = 1.3, alpha = 0.65,
                 show.legend = FALSE) +
     scale_fill_manual(values = CI) +
     scale_color_manual(values = c("Infection (-)" = "#1A5276",
                                   "Infection (+)" = "#922B21")) +
     labs(x = "", y = "RA Factor 1") +
     annotate("text", x = 1.5, y = ymx + ysp * 1.5,
              label = sprintf("p = %.3f", wt_5f$p.value), size = 2.8) +
     scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
     tn(8) + theme(legend.position = "none"),
   "Fig5f_RA_factor_infection")

# Fig5g — RA Factor 1 by progression (paper Fig 5g)
dhl_ra <- as.numeric(ct_unique$CT_Delta_Healthy_Lung_pct[
  match(ra_ids_m, ct_unique$Sample_ID)])
prog_ra <- ifelse(dhl_ra < -10, "Progressor", "Stable")
ok_prog <- !is.na(prog_ra)
df_5g <- data.frame(
  Value = ra_factors[ok_prog, 1],
  Prog = factor(prog_ra[ok_prog], levels = c("Stable", "Progressor"))
)
wt_5g <- wilcox.test(df_5g$Value[df_5g$Prog == "Progressor"],
                     df_5g$Value[df_5g$Prog == "Stable"], exact=TRUE)
cat(sprintf("  RA Factor 1 progression: p=%.4f\n", wt_5g$p.value))
ymx <- max(df_5g$Value); ysp <- (ymx - min(df_5g$Value)) * 0.12
sv(ggplot(df_5g, aes(Prog, Value, fill = Prog)) +
     geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25,
                  alpha = 0.7) +
     geom_jitter(aes(color = Prog), width = 0.08, size = 1.3, alpha = 0.65,
                 show.legend = FALSE) +
     scale_fill_manual(values = CP) +
     scale_color_manual(values = c("Stable" = "#1A5276",
                                   "Progressor" = "#922B21")) +
     labs(x = "", y = "RA Factor 1") +
     annotate("text", x = 1.5, y = ymx + ysp * 1.5,
              label = sprintf("p = %.3f", wt_5g$p.value), size = 2.8) +
     scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
     tn(8) + theme(legend.position = "none"),
   "Fig5g_RA_factor_progression")

# Fig5h — RA Factor 1 vs DeltaHL (paper Fig 5h)
ok_ct <- !is.na(ra_factors[, 1]) & !is.na(dhl_ra)
sp_5h <- cor.test(ra_factors[ok_ct, 1], dhl_ra[ok_ct], method = "spearman", exact = TRUE)
cat(sprintf("  RA Factor 1 vs DeltaHL: rho=%.3f, p=%.4f\n",
            sp_5h$estimate, sp_5h$p.value))
df_5h <- data.frame(Fval = ra_factors[ok_ct, 1], DHL = dhl_ra[ok_ct])
sv(ggplot(df_5h, aes(Fval, DHL)) +
     geom_point(size = 1.8, color = "#C0392B", alpha = 0.7) +
     geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.4,
                 fill = "grey90", formula = y ~ x) +
     geom_hline(yintercept = 0, linetype = "dotted", color = "grey60",
                linewidth = 0.2) +
     annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
              size = 3.5,
              label = paste0("\u03C1=", sprintf("%.3f", sp_5h$estimate),
                             "\np=", sprintf("%.4f", sp_5h$p.value))) +
     labs(x = "RA Factor 1", y = "\u0394 Healthy Lung (%)") +
     tn(8) + theme(legend.position = "none"),
   "Fig5h_RA_factor_CT")

# Fig5i — RA Factor 1 weights (paper Fig 5i)
fw_ra <- data.frame()
for (v in names(ra_weights)) {
  wv <- ra_weights[[v]][, 1]
  fw_ra <- rbind(fw_ra, data.frame(Feature = names(wv),
                                    Weight = as.numeric(wv),
                                    View = v, stringsAsFactors = FALSE))
}
fw_ra$AbsW <- abs(fw_ra$Weight)
fw_ra_top <- fw_ra %>% arrange(desc(AbsW)) %>% head(15)
fw_ra_top$Feature_clean <- gsub("^BALF_|^Serum_|^PB_", "", fw_ra_top$Feature)
fw_ra_top$Feature_clean <- gsub("_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14pos.*",
                                "", fw_ra_top$Feature_clean)
fw_ra_top$Feature_clean <- gsub("^GCSF$", "G-CSF", fw_ra_top$Feature_clean)
fw_ra_top$Feature_clean <- gsub("^IL(\\d)", "IL-\\1", fw_ra_top$Feature_clean)
fw_ra_top$Feature_clean <- gsub("^IFNa2$", "IFN-\u03B12",
                                fw_ra_top$Feature_clean)
fw_ra_top$Feature_clean <- fct_reorder(fw_ra_top$Feature_clean,
                                        fw_ra_top$AbsW)
fw_ra_top$View <- recode(fw_ra_top$View,
                         "BALF_FACS" = "BAL FCM", "PB_FACS" = "PB FCM",
                         "BALF_Cytokine" = "BAL Cytokine",
                         "Serum_Cytokine" = "Serum Cytokine")
sv(ggplot(fw_ra_top, aes(Weight, Feature_clean, fill = View)) +
     geom_col(width = 0.6) +
     scale_fill_manual(values = c("Expression" = "#E67E22",
                                  "BAL Cytokine" = "#3498DB",
                                  "Serum Cytokine" = "#E74C3C",
                                  "BAL FCM" = "#27AE60",
                                  "PB FCM" = "#9B59B6")) +
     geom_vline(xintercept = 0, linewidth = 0.3) +
     labs(x = "RA Factor 1 weight", y = "", fill = "") +
     tn(8) + theme(legend.position = "bottom",
                   legend.text = element_text(size = 6),
                   axis.text.y = element_text(size = 7)),
   "Fig5i_RA_factor_weights", 90, 80)

# Fig5j — RA progression ROC (paper Fig 5j)
ok_prog_m <- !is.na(prog_mofa) & !is.na(pb_th17_vals)
roc_th17_5j <- roc(prog_mofa[ok_prog_m], pb_th17_vals[ok_prog_m],
                   quiet = TRUE, direction = "auto")
fv_best <- ra_factors[, 1]
ok_fb <- !is.na(prog_mofa) & !is.na(fv_best)
roc_fac_5j <- roc(prog_mofa[ok_fb], fv_best[ok_fb],
                  quiet = TRUE, direction = "auto")
cat(sprintf("  PB Th17 AUC: %.3f, RA Factor 1 AUC: %.3f\n",
            auc(roc_th17_5j), auc(roc_fac_5j)))

ep1 <- add_ep(1 - roc_th17_5j$specificities, roc_th17_5j$sensitivities)
ep2 <- add_ep(1 - roc_fac_5j$specificities, roc_fac_5j$sensitivities)
roc_5j <- rbind(
  data.frame(Sens = ep1$sens, Spec = ep1$spec,
             Model = sprintf("PB Th17 (AUC = %.3f)", auc(roc_th17_5j))),
  data.frame(Sens = ep2$sens, Spec = ep2$spec,
             Model = sprintf("RA Factor 1 (AUC = %.3f)", auc(roc_fac_5j)))
)
sv(ggplot(roc_5j, aes(Spec, Sens, color = Model)) +
     geom_step(linewidth = 0.7, direction = "vh") +
     geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                 color = "grey60", linewidth = 0.3) +
     scale_color_manual(values = c("#C0392B", "#9B59B6")) +
     labs(x = "1 \u2013 Specificity", y = "Sensitivity", color = "") +
     ggtitle("CT Progression (RA)") +
     tn(8) + theme(legend.position = c(0.65, 0.22),
                   legend.background = element_rect(
                     fill = alpha("white", 0.95), color = NA),
                   legend.text = element_text(size = 6.5),
                   legend.key.height = unit(3, "mm")),
   "Fig5j_RA_progression_ROC", 80, 75)


# =================================================================
#  SUPPLEMENTARY FIGURES
# =================================================================
cat("\n========== SUPPLEMENTARY FIGURES ==========\n")

# SFig1a-b — scRNA-seq UMAPs (generated in 03_Analysis.R)
# These require the Seurat reference object and are already saved in
# results/figures/. Copy to panels/ if needed.
if (file.exists("results/figures/FigS1c_UMAP_fine.pdf") &&
    !file.exists(file.path(panel_dir, "FigS1c_UMAP_fine.pdf"))) {
  file.copy("results/figures/FigS1c_UMAP_fine.pdf",
            file.path(panel_dir, "FigS1c_UMAP_fine.pdf"))
  file.copy("results/figures/FigS1c_UMAP_fine.png",
            file.path(panel_dir, "FigS1c_UMAP_fine.png"))
  cat("  Copied FigS1c_UMAP_fine\n")
}
if (file.exists("results/figures/FigS1d_UMAP_source.pdf") &&
    !file.exists(file.path(panel_dir, "FigS1d_UMAP_source.pdf"))) {
  file.copy("results/figures/FigS1d_UMAP_source.pdf",
            file.path(panel_dir, "FigS1d_UMAP_source.pdf"))
  file.copy("results/figures/FigS1d_UMAP_source.png",
            file.path(panel_dir, "FigS1d_UMAP_source.png"))
  cat("  Copied FigS1d_UMAP_source\n")
}

# SFig1c — Bland-Altman (already in panels/)
cat("  FigS1c_BlandAltman: pre-generated (01_BayesPrism)\n")

# SFig2a — GSEA infection (paper SFig 2a)
ginf <- tryCatch(read.csv("results/tables/GSEA_Infection_GO_BP.csv"),
                 error = function(e) NULL)
if (!is.null(ginf)) {
  tgi <- head(ginf[order(ginf$p.adjust), ], 8) %>%
    mutate(D = str_trunc(Description, 42))
  tgi <- tgi[order(tgi$NES), ]; tgi$D <- fct_inorder(tgi$D)
  sv(ggplot(tgi, aes(NES, D, fill = NES > 0)) +
       geom_col(width = 0.55, show.legend = FALSE) +
       scale_fill_manual(values = c("TRUE" = "#C0392B",
                                    "FALSE" = "#2980B9")) +
       geom_vline(xintercept = 0, linewidth = 0.2) +
       labs(x = "NES", y = "") +
       tn(7) + theme(axis.text.y = element_text(size = 6)),
     "FigS2a_GSEA_infection", 100, 70)
}

# SFig2b — Microbiome infection (paper SFig 2b)
mi <- tryCatch({
  md2 <- read.csv("results/tables/Microbiome_Infection_DiffAbundance.csv")
  md2$Genus <- sapply(strsplit(as.character(md2$Taxon), ";"),
                      function(x) trimws(tail(x, 2)[1]))
  md2 %>% group_by(Genus) %>% dplyr::slice(1) %>% ungroup() %>%
    dplyr::select(Genus, Log2FC, P_value, P_adjusted) %>% as.data.frame()
}, error = function(e) NULL)
if (!is.null(mi)) {
  tmi <- head(mi %>% arrange(P_value), 8) %>% arrange(Log2FC)
  tmi$Genus <- fct_inorder(tmi$Genus)
  tmi$sg <- ifelse(tmi$P_value < 0.05, "p < 0.05", "p >= 0.05")
  sv(ggplot(tmi, aes(Log2FC, Genus, fill = sg)) +
       geom_col(width = 0.5) +
       scale_fill_manual(values = c("p < 0.05" = "#C0392B",
                                    "p >= 0.05" = "#95A5A6")) +
       geom_vline(xintercept = 0, linewidth = 0.2) +
       labs(x = "log2FC(Infection+ vs \u2212)", y = "", fill = "") +
       tn(7) + theme(axis.text.y = element_text(size = 6, face = "italic"),
                     legend.position = "bottom",
                     legend.text = element_text(size = 6)),
     "FigS2b_microbiome", 90, 70)
}

# SFig2c — Network (paper SFig 2c)
if (file.exists("results/Integration_Results.RData")) {
  load("results/Integration_Results.RData")
  edge_sum <- edges %>%
    mutate(Pair = paste(pmin(Layer1, Layer2), "-", pmax(Layer1, Layer2))) %>%
    group_by(Pair) %>% summarise(N = n(), .groups = "drop") %>%
    arrange(desc(N))
  edge_sum$Pair <- fct_reorder(edge_sum$Pair, edge_sum$N)
  sv(ggplot(edge_sum, aes(N, Pair)) +
       geom_col(fill = "#34495E", width = 0.55) +
       labs(x = "Edges (|rho|>0.5)", y = "") +
       tn(8) + theme(axis.text.y = element_text(size = 6)),
     "FigS2c_network", 100, 80)
}

# SFig2d — GSVA overview (paper SFig 2d)
sv(ggplot(glg, aes(S_d, PW, fill = Sc)) +
     geom_tile() +
     scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                          midpoint = 0, name = "Score") +
     labs(x = "", y = "") +
     tn(6) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 3),
                   axis.text.y = element_text(size = 6)),
   "FigS2d_GSVA_overview", 120, 70)

# SFig2e — AUC summary (paper SFig 2e)
auc_comp <- data.frame(
  Model = c("RA vs Sarcoidosis\n(multi-omics)", "BALF Th17.1\n(infection)",
            "BAL/PB Th17.1\nratio", "GSVA OXPHOS", "BALF MCP3"),
  AUC = c(0.962, 0.870, 0.861, 0.861, 0.866))
auc_comp$Model <- fct_inorder(auc_comp$Model)
sv(ggplot(auc_comp, aes(AUC, Model)) +
     geom_col(fill = "#C0392B", width = 0.5) +
     geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60",
                linewidth = 0.2) +
     labs(x = "AUC", y = "") + xlim(0, 1) +
     tn(8) + theme(axis.text.y = element_text(size = 7)),
   "FigS2e_AUC_summary")

# SFig2f — CT ILD stratified (paper SFig 2f)
ra_ct_all <- ct_unique[ct_unique$Sample_Group == "RA" &
                          !ct_unique$Sample_ID %in% EXCL, ]
ild_status <- ifelse(ra_ct_all$ILD_Score > 0 |
                       grepl("ILD", ra_ct_all$Subgroup, ignore.case = TRUE),
                     "ILD+", "ILD-")
bm_ct_cols <- c("PB_Th17_CCR6pos_CXCR3neg", "PB_Th17.1_CCR6pos_CXCR3pos",
                "BALF_Neutrophil_Percent")
bm_ct_labs <- c("PB Th17", "PB Th17.1", "BAL Neutrophil")
plots_s2f <- list()
for (bi in seq_along(bm_ct_cols)) {
  vals <- as.numeric(ra_ct_all[[bm_ct_cols[bi]]])
  dhl_s <- as.numeric(ra_ct_all$CT_Delta_Healthy_Lung_pct)
  ok <- !is.na(vals) & !is.na(dhl_s) & !is.na(ild_status)
  df_s <- data.frame(BM = vals[ok], DHL = dhl_s[ok], ILD = ild_status[ok])
  sp_ild <- by(df_s, df_s$ILD, function(d)
    cor.test(d$BM, d$DHL, method = "spearman", exact = TRUE))
  plots_s2f[[bi]] <- ggplot(df_s, aes(BM, DHL, color = ILD)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4,
                formula = y ~ x) +
    scale_color_manual(values = c("ILD+" = "#C0392B", "ILD-" = "#3498DB")) +
    labs(x = bm_ct_labs[bi], y = "\u0394 Healthy Lung (%)", color = "") +
    tn(7) + theme(legend.position = "bottom")
}
if (length(plots_s2f) > 0) {
  sv(wrap_plots(plots_s2f, nrow = 1),
     "FigS2f_CT_ILD_stratified", 180, 70)
}

# SFig3 — Healthy serum reference (paper SFig 3)
cat("\n--- Supplementary Fig. 3 ---\n")
col_healthy <- "#27AE60"
cyto_info <- list(
  list(col = "Serum_GCSF", label = "G-CSF"),
  list(col = "Serum_IL6", label = "IL-6"),
  list(col = "Serum_IL13", label = "IL-13"),
  list(col = "Serum_IL17A", label = "IL-17A"),
  list(col = "Serum_IFNa2", label = "IFN-\u03b12"),
  list(col = "Serum_IL22", label = "IL-22"),
  list(col = "Serum_IL4", label = "IL-4"),
  list(col = "Serum_IL2", label = "IL-2"),
  list(col = "Serum_GROa", label = "GRO-\u03b1"),
  list(col = "Serum_FGF2", label = "FGF-2")
)
plots_s3 <- list()
for (ci in seq_along(cyto_info)) {
  info <- cyto_info[[ci]]
  vals <- as.numeric(master_data[[info$col]])
  grp <- ifelse(master_data$Sample_Group == "RA", "RA", "Sarcoidosis")
  ok <- !is.na(vals)
  df <- data.frame(Value = vals[ok],
                   Group = factor(grp[ok], levels = c("RA", "Sarcoidosis")))
  wt <- wilcox.test(Value ~ Group, data = df, exact=TRUE)
  p_lab <- fmt_p(wt$p.value)

  cyto_short <- gsub("^Serum_", "", info$col)
  row_idx <- which(s3_table$Cytokine == cyto_short)
  h_med <- h_lo <- h_hi <- NA
  if (length(row_idx) == 1) {
    hm_str <- s3_table$Healthy_median_95CI[row_idx]
    hm_vals <- as.numeric(strsplit(gsub("[][,]", "", hm_str), " +")[[1]])
    h_med <- hm_vals[1]; h_lo <- hm_vals[2]; h_hi <- hm_vals[3]
  }
  ymax <- max(df$Value, h_hi, na.rm = TRUE) * 1.25

  p <- ggplot(df, aes(Group, Value, fill = Group)) +
    {if (!is.na(h_lo) && !is.na(h_hi))
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = h_lo, ymax = h_hi,
               fill = col_healthy, alpha = 0.15)} +
    {if (!is.na(h_med))
      geom_hline(yintercept = h_med, linetype = "dashed",
                 color = col_healthy, linewidth = 0.4)} +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3,
                 alpha = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = c("RA" = "#C0392B",
                                 "Sarcoidosis" = "#2980B9")) +
    labs(x = "", y = "pg/mL", title = info$label) +
    annotate("text", x = 1.5, y = ymax, label = p_lab, size = 2.2,
             fontface = "italic") +
    coord_cartesian(ylim = c(0, ymax * 1.05), clip = "off") +
    tn(7) + theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5),
                  axis.text.x = element_text(size = 6))
  plots_s3[[ci]] <- p
}
sv(wrap_plots(plots_s3, ncol = 5) +
     plot_annotation(theme = theme(plot.margin = ggplot2::margin(2, 2, 2, 2))),
   "FigS3_healthy_serum_95CI", 250, 110)

# SFig4 — Medication confounding (paper SFig 4)
cat("\n--- Supplementary Fig. 4 ---\n")
ra_samples_m <- rownames(factors)[grepl("^(KY[0-9]|RA[0-9])", rownames(factors))]
ra_samples_m <- setdiff(ra_samples_m, EXCL)
med_cols <- c("steroid.1.0.", "MTX", "NSAIDS")
med_labels <- c("GC", "MTX", "NSAIDs")
factor_names <- c("IDD", "IIV")
plots_s4 <- list(); pi <- 1
for (fi in 1:2) {
  fn <- factor_names[fi]
  fac_vals <- factors[ra_samples_m, fi]
  for (mi in seq_along(med_cols)) {
    mc <- med_cols[mi]; ml <- med_labels[mi]
    med_vals <- ct_unique[[mc]][match(ra_samples_m, ct_unique$Sample_ID)]
    ok <- !is.na(fac_vals) & !is.na(med_vals)
    on_vals <- fac_vals[ok & med_vals == 1]
    off_vals <- fac_vals[ok & med_vals == 0]
    if (length(on_vals) < 2 || length(off_vals) < 2) next
    wt <- wilcox.test(on_vals, off_vals, exact=TRUE)
    p_lab <- fmt_p(wt$p.value)
    df <- data.frame(
      Value = c(on_vals, off_vals),
      Med = factor(c(rep(paste0(ml, " (+)"), length(on_vals)),
                     rep(paste0(ml, " (\u2212)"), length(off_vals))),
                   levels = c(paste0(ml, " (\u2212)"), paste0(ml, " (+)"))))
    yrange <- range(df$Value)
    ymax <- yrange[2] + diff(yrange) * 0.25
    p <- ggplot(df, aes(Med, Value, fill = Med)) +
      geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3,
                   alpha = 0.7) +
      geom_jitter(width = 0.12, size = 0.8, alpha = 0.6,
                  show.legend = FALSE) +
      scale_fill_manual(values = c("#2980B9", "#C0392B")) +
      labs(x = "", y = fn, title = ml) +
      annotate("text", x = 1.5, y = ymax, label = p_lab, size = 2.2,
               fontface = "italic") +
      coord_cartesian(ylim = c(min(yrange[1], 0), ymax * 1.05),
                      clip = "off") +
      tn(7) + theme(legend.position = "none",
                    axis.text.x = element_text(size = 6))
    plots_s4[[pi]] <- p; pi <- pi + 1
  }
}
sv(wrap_plots(plots_s4, ncol = 3) +
     plot_annotation(theme = theme(plot.margin = ggplot2::margin(2, 2, 2, 2))),
   "FigS4_medication", 160, 100)

# SFig5 — Ligand-cell associations (paper SFig 5)
cat("\n--- Supplementary Fig. 5 ---\n")
balf_cyto_cols <- grep("^BALF_", colnames(master_data), value = TRUE)
balf_cyto_cols <- balf_cyto_cols[!grepl(
  "Th17|Treg|Th1|Th2|Macro|Neutro|CD|Lymph|Eosino|Baso|Plasma|Other",
  balf_cyto_cols)]
balf_cyto_cols <- balf_cyto_cols[sapply(balf_cyto_cols,
                                        function(x) is.numeric(master_data[[x]]))]
if (length(balf_cyto_cols) > 15) balf_cyto_cols <- balf_cyto_cols[1:15]
cell_types_lc <- intersect(c("Macrophage", "T_cell", "NK", "B_cell", "Plasma",
                              "Neutrophil", "DC"), colnames(bp_cell))
cor_mat <- matrix(NA, nrow = length(balf_cyto_cols), ncol = length(cell_types_lc))
rownames(cor_mat) <- balf_cyto_cols; colnames(cor_mat) <- cell_types_lc
pval_mat_lc <- cor_mat
common_ids <- intersect(rownames(bp_cell), master_data$Sample_ID)
for (i in seq_along(balf_cyto_cols)) {
  cyto_vals <- as.numeric(master_data[[balf_cyto_cols[i]]][
    match(common_ids, master_data$Sample_ID)])
  for (j in seq_along(cell_types_lc)) {
    cell_vals <- bp_cell[common_ids, cell_types_lc[j]]
    v <- !is.na(cyto_vals) & !is.na(cell_vals)
    if (sum(v) >= 10) {
      sp <- cor.test(cyto_vals[v], cell_vals[v], method = "spearman", exact = TRUE)
      cor_mat[i, j] <- sp$estimate; pval_mat_lc[i, j] <- sp$p.value
    }
  }
}
rownames(cor_mat) <- gsub("^BALF_", "", rownames(cor_mat))
cor_long <- as.data.frame(cor_mat) %>% rownames_to_column("Cytokine") %>%
  pivot_longer(-Cytokine, names_to = "CellType", values_to = "Rho")
pval_long_lc <- as.data.frame(pval_mat_lc) %>%
  rownames_to_column("Cytokine") %>%
  pivot_longer(-Cytokine, names_to = "CellType", values_to = "P")
cor_long$P <- pval_long_lc$P
cor_long$Sig <- ifelse(cor_long$P < 0.05, "*", "")
sv(ggplot(cor_long, aes(CellType, Cytokine, fill = Rho)) +
     geom_tile(color = "white", linewidth = 0.5) +
     geom_text(aes(label = Sig), size = 3, color = "black", vjust = 0.8) +
     scale_fill_gradient2(low = "#2980B9", mid = "white", high = "#C0392B",
                          midpoint = 0, limits = c(-1, 1),
                          name = "Spearman \u03C1") +
     labs(x = "BayesPrism cell type", y = "BAL cytokine") +
     tn(8) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                   axis.text.y = element_text(size = 6)),
   "FigS5_ligand_cell", 110, 100)


# =================================================================
cat("\n========== ALL PANELS COMPLETE ==========\n")
cat(sprintf("Output directory: %s\n", normalizePath(panel_dir)))
# =================================================================
