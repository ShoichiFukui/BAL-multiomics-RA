#!/usr/bin/env Rscript
# =================================================================
# Comprehensive Figure Generation Script
# BAL multi-omics — the accompanying manuscript
# =================================================================
# Generates all figure panels from pre-computed analysis results.
# Output panel names match PPT v21 assembly.
#
# Prerequisites (run in order):
#   01_BayesPrism_Deconvolution.R -> output_v3_BayesPrism/
#   02_PostDeconvolution.R        -> results/RA_ILD_Workspace.RData
#   analysis_modules/Enhanced_Analysis.R -> results/tables/*.csv
#   analysis_modules/CT_Multiomics.R     -> results/CT_Multiomics_Results.RData

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

BASEDIR <- Sys.getenv("PROJECT_DIR", unset = getwd())
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
               "Serum_Cytokine" = "#E74C3C", "BALF_FCM" = "#27AE60",
               "PB_FCM" = "#9B59B6", "Microbiome" = "#7F8C8D")
view_labels <- c("Expression" = "Expression", "BALF_Cytokine" = "BAL Cytokine",
                 "Serum_Cytokine" = "Serum Cytokine", "BALF_FCM" = "BAL FCM",
                 "PB_FCM" = "PB FCM", "Microbiome" = "Microbiome")

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
  mutate(Group = ifelse(grepl("^Sarcoidosis", SID), "Sarcoidosis", "RA"),
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
   "Fig1b_UMAP_dataset", 85, 75)

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
   "Fig1c_UMAP_celltype", 90, 85)

# Fig1d — Dotplot (paper Fig 1d)
dd$M <- factor(dd$M, levels = rev(mks))
sv(ggplot(dd, aes(CT, M, size = Pct, color = Avg)) +
     geom_point() +
     scale_size_continuous(range = c(0.5, 3.5), name = "%Exp") +
     scale_color_gradient(low = "grey90", high = "#C0392B", name = "Avg") +
     labs(x = "", y = "") +
     tn(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                   axis.text.y = element_text(size = 6, face = "italic")),
   "Fig1d_dotplot", 100, 80)

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
   "Fig1e_composition", 110, 75)

rm(ref, ed); gc()

# Fig1f — BayesPrism stacked bar (paper Fig 1f)
bp_l <- bp_cell %>% rownames_to_column("S") %>%
  mutate(G = ifelse(grepl("^Sarcoidosis", S), "Sarcoidosis", "RA"),
         S_d = ifelse(grepl("^Sarcoidosis", S),
                      paste0("Sarcoidosis", gsub("^Sarcoidosis", "", S)),
                      paste0("RA", gsub("^RA", "", S)))) %>%
  pivot_longer(-c(S, G, S_d), names_to = "CT", values_to = "Fr") %>%
  mutate(S_d = fct_reorder(S_d, as.numeric(factor(G))))
sv(ggplot(bp_l, aes(S_d, Fr, fill = CT)) +
     geom_col(width = 0.85) + scale_fill_manual(values = CC) +
     labs(x = "", y = "Proportion", fill = "") +
     tn(7) + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 4),
                   legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = element_text(size = 6)),
   "Fig1f_stacked_bar", 140, 70)


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
   "Fig2a_PCA_subgroup")

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
   "Fig2b_volcano", 115, 90)

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
   "Fig2c_GSEA_GO", 110, 85)

# Fig2d — GSEA KEGG (paper Fig 2d)
tkk <- head(gkk[order(gkk$p.adjust), ], 10) %>%
  mutate(D = str_trunc(Description, 44))
tkk <- tkk[order(tkk$NES), ]; tkk$D <- fct_inorder(tkk$D)
sv(ggplot(tkk, aes(NES, D, fill = NES > 0)) +
     geom_col(width = 0.6, show.legend = FALSE) +
     scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "NES", y = "") + tn(7) + theme(axis.text.y = element_text(size = 8)),
   "Fig2d_GSEA_KEGG", 110, 75)

# Fig2e — GSVA heatmap (paper Fig 2e)
gm <- gsva_scores[, common_samples]
so <- common_samples[order(meta$Sample_Group[match(common_samples, meta$Sample_ID)])]
glg <- as.data.frame(gm) %>% rownames_to_column("PW") %>%
  pivot_longer(-PW, names_to = "S", values_to = "Sc")
glg$S_d <- ifelse(grepl("^Sarcoidosis", glg$S),
                  paste0("Sarcoidosis", gsub("^Sarcoidosis", "", glg$S)),
                  paste0("RA", gsub("^RA", "", glg$S)))
so_d <- ifelse(grepl("^Sarcoidosis", so),
               paste0("Sarcoidosis", gsub("^Sarcoidosis", "", so)),
               paste0("RA", gsub("^RA", "", so)))
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
   "Fig2e_GSVA_heatmap", 120, 85)

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
ga_pl <- ifelse(grepl("^Sarcoidosis", rownames(bp_cell)), "Sarcoidosis", "RA")
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
   "Fig2i_Th_subsets", 100, 75)

# Fig2j-k — BAL vs PB (paper Fig 2j-k)
sv(bvp("BALF_Th17.1_CCR6pos_CXCR3pos", "PB_Th17.1_CCR6pos_CXCR3pos", "Th17.1"),
   "Fig2j_BAL_PB_Th171")
sv(bvp("BALF_Treg_CD127neg_CD25pos", "PB_Treg_CD127neg_CD25pos", "Treg"),
   "Fig2k_BAL_PB_Treg")

# Fig2l — Fully fold-internal nested LOOCV ROC (paper Fig 2l; AUC = 0.780)
# Leakage-free: gene pre-selection (top-50 by TRAINING variance) and supervised
# feature ranking are performed WITHIN each LOOCV fold; no full-cohort prefilter
# (which would inflate AUC to ~0.962). Reproduces 03_Analysis.R.
{
  Gl <- t(expr_matrix)
  sl <- intersect(intersect(rownames(Gl), rownames(cyto_mat)), rownames(fcm_mat))
  gl <- factor(meta$Sample_Group[match(sl, meta$Sample_ID)])
  Gl <- Gl[sl, , drop = FALSE]; Cl <- cyto_mat[sl, , drop = FALSE]; Fl <- fcm_mat[sl, , drop = FALSE]
  nL <- length(sl); NG <- 50
  cV <- function(M) { nr <- nrow(M); cs <- colSums(M); (colSums(M * M) - cs * cs / nr) / (nr - 1) }
  foldL <- vector("list", nL)
  for (i in 1:nL) {
    tr <- setdiff(1:nL, i)
    topg <- names(sort(cV(Gl[tr, , drop = FALSE]), decreasing = TRUE))[1:NG]
    X <- cbind(Gl[, topg, drop = FALSE], Cl, Fl)
    med <- apply(X[tr, , drop = FALSE], 2, function(z) { z[is.infinite(z)] <- NA; median(z, na.rm = TRUE) })
    Xi <- X
    for (j in seq_len(ncol(Xi))) { z <- Xi[, j]; z[is.na(z) | is.infinite(z)] <- med[j]; Xi[, j] <- z }
    tv <- cV(Xi[tr, , drop = FALSE]); keep <- names(tv[!is.na(tv) & tv > 0])
    foldL[[i]] <- list(tr = tr, X = Xi[, keep, drop = FALSE])
  }
  mwp <- function(Xtr, g) {
    lev <- levels(g); g1 <- g == lev[1]; n1 <- sum(g1); n2 <- sum(!g1); N <- n1 + n2
    R <- apply(Xtr, 2, rank); R1 <- colSums(R[g1, , drop = FALSE]); U <- R1 - n1 * (n1 + 1) / 2; mu <- n1 * n2 / 2
    tie <- apply(Xtr, 2, function(x) { tt <- table(x); sum(tt^3 - tt) })
    sig <- sqrt((n1 * n2 / 12) * ((N + 1) - tie / (N * (N - 1)))); z <- (U - mu) / sig; 2 * pnorm(-abs(z))
  }
  rl <- function(labels) {
    prob <- numeric(nL)
    for (i in 1:nL) {
      f <- foldL[[i]]; tr <- f$tr; X <- f$X; ytr <- labels[tr]
      p <- mwp(X[tr, , drop = FALSE], ytr); top <- names(sort(p))[1:min(50, length(p))]
      rf <- randomForest(x = X[tr, top, drop = FALSE], y = ytr, ntree = 500)
      prob[i] <- stats::predict(rf, X[i, top, drop = FALSE], type = "prob")[, "RA"]
    }
    prob
  }
  set.seed(42); prob_obs <- rl(gl)
  roc1 <- roc(gl, prob_obs, levels = c("Control", "RA"), direction = "<", quiet = TRUE)
  set.seed(42); ci1 <- ci.auc(roc1, method = "bootstrap", boot.n = 2000, quiet = TRUE)
  rd1 <- data.frame(S = roc1$sensitivities, Sp = 1 - roc1$specificities); rd1 <- rd1[order(rd1$Sp, rd1$S), ]
  sv(ggplot(rd1, aes(Sp, S)) +
       geom_step(color = "#C0392B", linewidth = 0.8, direction = "vh") +
       geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
       annotate("text", x = 0.55, y = 0.1, size = 2.8, hjust = 0,
                label = sprintf("AUC=%.3f\n95%%CI:%.3f-%.3f", as.numeric(auc(roc1)), ci1[1], ci1[3])) +
       labs(x = "1-Specificity", y = "Sensitivity") +
       coord_equal() + tn(8),
     "Fig2l_ROC_foldinternal")
}


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
   "Fig3a_PCA_infection")

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
   "Fig3b_Cliff_delta", 95, 95)

# Fig3c-e — Infection boxplots (paper Fig 3c-e)
sv(bxi("BALF_Th17.1_CCR6pos_CXCR3pos", "BALF Th17.1 (%)", ra_md),
   "Fig3c_Th171_infection")
sv(bxi("BALF_Treg_CD127neg_CD25pos", "BALF Treg (%)", ra_md),
   "Fig3d_Treg_infection")

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
   "Fig3e_ratio_infection")

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
   "Fig3g_surfactant", 90, 90)

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
   "Fig3h_OXPHOS_infection")


# =================================================================
#  FIGURE 4: CT-quantified lung disease progression
# =================================================================
cat("\n========== FIGURE 4 ==========\n")
# Fig4a: manual CT images (not generated)

# Fig4b — ILD vs GGO (paper Fig 4b)
sv(sct(as.numeric(ra_ext2$ILD_Score),
       as.numeric(ra_ext2$CT_Delta_GGO_pct),
       "ILD score", "Delta GGO (%)"),
   "Fig4b_ILD_vs_GGO")

# Fig4c — PB Th17.1 vs DeltaHL (paper Fig 4c)
sv(sct(as.numeric(ra_ext2$`PB_Th17.1_CCR6pos_CXCR3pos`),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "PB Th17.1 (%)", "Delta Healthy Lung (%)"),
   "Fig4d_PBTh171_vs_DeltaHL")

# Fig4d — PB Th17 vs DeltaHL (paper Fig 4d)
sv(sct(as.numeric(ra_ext2$`PB_Th17_CCR6pos_CXCR3neg`),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "PB Th17 (%)", "Delta Healthy Lung (%)"),
   "Fig4c_PBTh17_vs_DeltaHL")

# Fig4e — Neutrophil vs DeltaHL (paper Fig 4e)
sv(sct(as.numeric(ra_ext2$BALF_Neutrophil_Percent),
       as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),
       "BALF Neutrophil (%)", "Delta Healthy Lung (%)"),
   "Fig4e_Neutrophil_vs_DeltaHL")

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

# Fig5e — RA-only MOFA2 R2 heatmap (paper Fig 5e)
r2m_ra <- ra_r2$r2_per_factor[[1]]
r2l_ra <- as.data.frame(r2m_ra) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to = "View", values_to = "R2")
r2l_ra$View <- recode(r2l_ra$View,
                      "BALF_FCM" = "BAL FCM", "PB_FCM" = "PB FCM",
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
   "Fig5d_RA_MOFA2_R2", 100, 65)

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
   "Fig5e_RA_factor_infection")

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
   "Fig5i_RA_progression_ROC", 80, 75)


# =================================================================
#  SUPPLEMENTARY FIGURES
# =================================================================
cat("\n========== SUPPLEMENTARY FIGURES ==========\n")


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
     "SFig3a_GSEA_infection", 100, 70)
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
     "SFig3b_microbiome", 90, 70)
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
   "SFig2_healthy_serum_95CI", 250, 110)

# SFig4 — Medication confounding (paper SFig 4)
cat("\n--- Supplementary Fig. 4 ---\n")
ra_samples_m <- rownames(factors)[grepl("^RA[0-9]", rownames(factors))]
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
   "SFig4_medication", 160, 100)



# =================================================================
cat("\n========== ALL PANELS COMPLETE ==========\n")
cat(sprintf("Output directory: %s\n", normalizePath(panel_dir)))
# =================================================================


# ============================================================================
# Final figure/output step
# ============================================================================
# =====================================================================
# Supplementary Fig. 1 (deconvolution validation): 3-panel Bland-Altman
# comparing BayesPrism deconvolution vs paired flow cytometry for
# macrophages (a), lymphocytes (b), neutrophils (c). Supersedes the
# earlier 5-panel layout (the scRNA-seq UMAP panels were removed).
# Run after 02_PostDeconvolution.R. Assumes cwd = parent project dir.
# =====================================================================
suppressMessages({library(tidyverse); library(readxl); library(patchwork)})
EXCL <- character(0)  # released dataset is the final n = 35 cohort (no-op filter)
bp_cell <- read.csv("output_v3_BayesPrism/celltype_proportions_cellFraction.csv", row.names = 1)
bp_cell <- bp_cell[!rownames(bp_cell) %in% EXCL, ]
fcm <- read_excel("FCM_integrated_data_transformed.xlsx", sheet = "BALF_analysis")
cf <- intersect(rownames(bp_cell), fcm$Sample_ID)
fi_idx <- match(cf, fcm$Sample_ID)
lc <- intersect(c("T_cell", "NK", "B_cell", "Plasma"), colnames(bp_cell))
panel_dir <- "results/panels"; dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.1), face = "bold", hjust = 0,
                                    margin = ggplot2::margin(0, 0, 3, 0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 3), "mm"))
}
ba_panel <- function(D, F, lb, ccc, tag) {
  v <- !is.na(D) & !is.na(F); D <- D[v]; F <- F[v]
  # Concordance correlation coefficient (Lin) computed from data (overrides hardcoded arg)
  ccc <- 2 * cor(D, F) * sd(D) * sd(F) / (var(D) + var(F) + (mean(D) - mean(F))^2)
  m <- (D + F) / 2; d <- D - F
  bias <- mean(d); s <- sd(d); lo <- bias - 1.96 * s; hi <- bias + 1.96 * s
  df <- data.frame(m = m, d = d, G = ifelse(grepl("^Sarcoidosis", cf[v]), "Sarcoidosis", "RA"))
  ggplot(df, aes(m, d, color = G)) +
    geom_hline(yintercept = bias, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_hline(yintercept = c(lo, hi), linetype = "dotted", color = "#C0392B", linewidth = 0.4) +
    geom_point(size = 1.4, alpha = 0.8) +
    scale_color_manual(values = c("RA" = "#C0392B", "Sarcoidosis" = "#2980B9")) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.8,
             label = sprintf("CCC = %.3f", ccc)) +
    labs(title = paste0(tag, "   ", lb), x = "Mean of methods (%)",
         y = "Deconv. − FCM (%)", color = "") +
    tn(7) + theme(legend.position = "none")
}
pa <- ba_panel(bp_cell[cf, "Macrophage"] * 100, as.numeric(fcm$BALF_Macrophage_Percent[fi_idx]), "Macrophage", 0.041, "a")
pb <- ba_panel(rowSums(bp_cell[cf, lc]) * 100, as.numeric(fcm$BALF_Lymphocyte_Percent[fi_idx]), "Lymphocyte", 0.215, "b")
pc <- ba_panel(bp_cell[cf, "Neutrophil"] * 100, as.numeric(fcm$BALF_Neutrophil_Percent[fi_idx]), "Neutrophil", 0.563, "c")
fig <- pa + pb + pc + plot_layout(nrow = 1)
ggsave(file.path(panel_dir, "SFig1_BlandAltman_abc.png"), fig, width = 180, height = 62, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(panel_dir, "SFig1_BlandAltman_abc.pdf"), fig, width = 180, height = 62, units = "mm", dpi = 300, device = cairo_pdf)
cat(sprintf("Saved Supplementary Fig. 1 (Bland-Altman a-c); n pairs = %d\n", length(cf)))


# ============================================================================
# Final figure/output step
# ============================================================================
# =====================================================================
# Builds Supplementary Data 3: a single comprehensive infection biomarker
# table (all omics + top transcripts). For every biomarker it reports
# median(+/-), N, ROC AUC + bootstrap 95% CI, Cliff's delta, Wilcoxon p,
# BH-adjusted p, and direction. This unifies the former single-biomarker
# ROC table (AUC + 95% CI) with the all-omics effect-size table
# (Cliff's delta + Wilcoxon p), which are mathematically linked
# (Cliff's delta = 2*AUC - 1 for a two-group comparison).
# Run after 03_Analysis.R. Assumes cwd = parent project dir.
# =====================================================================
suppressMessages({library(effsize); library(dplyr); library(pROC)})
load("results/RA_ILD_Workspace.RData")
TBL <- "results/tables"
allo <- read.csv(file.path(TBL, "Infection_AllOmics_Comparison.csv"), check.names = FALSE)
roc  <- read.csv(file.path(TBL, "Infection_SingleBiomarker_ROC.csv"), check.names = FALSE)

ra <- master_data[master_data$Sample_Group == "RA" & !is.na(master_data$respiratory_infection), ]
ra$InfStat <- factor(ifelse(ra$respiratory_infection == 1, "pos", "neg"), levels = c("neg", "pos"))
samp <- ra$Sample_ID
cat(sprintf("RA infection subset: %d (pos=%d, neg=%d)\n",
            length(samp), sum(ra$InfStat == "pos"), sum(ra$InfStat == "neg")))

auc_ci <- function(vals, grp) {
  ok <- !is.na(vals); v <- vals[ok]; g <- grp[ok]
  if (length(unique(g)) < 2) return(c(AUC = NA, lo = NA, hi = NA))
  r <- roc(g, v, levels = c("neg", "pos"), direction = "auto", quiet = TRUE)
  set.seed(42)
  ci <- tryCatch(as.numeric(ci.auc(r, method = "bootstrap", boot.n = 2000, quiet = TRUE)),
                 error = function(e) c(NA, NA, NA))
  c(AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
}
stat_row <- function(vals, variable, category) {
  pos <- vals[ra$InfStat == "pos"]; neg <- vals[ra$InfStat == "neg"]
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (length(pos) < 3 || length(neg) < 3) return(NULL)
  wt <- wilcox.test(pos, neg, exact = TRUE)
  cd <- as.numeric(cliff.delta(pos, neg)$estimate)
  ac <- auc_ci(vals, ra$InfStat)
  data.frame(Variable = variable, Category = category,
             Median_pos = median(pos), Median_neg = median(neg),
             N_pos = length(pos), N_neg = length(neg),
             AUC = unname(ac["AUC"]), CI_lo = unname(ac["lo"]), CI_hi = unname(ac["hi"]),
             Cliff_delta = cd, P_value = wt$p.value,
             Direction = ifelse(cd > 0, "Up in Inf+", "Down in Inf+"), stringsAsFactors = FALSE)
}
get_vals <- function(variable, category) {
  if (category %in% c("BALF_Cytokine", "Serum_Cytokine", "FCM")) {
    if (!(variable %in% colnames(ra))) return(NULL); return(as.numeric(ra[[variable]]))
  }
  if (category == "Deconvolution") {
    if (!exists("deconv_mat") || !(variable %in% colnames(deconv_mat))) return(NULL)
    return(deconv_mat[match(samp, rownames(deconv_mat)), variable])
  }
  if (category == "GSVA_Pathway") {
    if (!exists("gsva_scores") || !(variable %in% rownames(gsva_scores))) return(NULL)
    return(gsva_scores[variable, match(samp, colnames(gsva_scores))])
  }
  if (category == "Microbiome") {
    if (!exists("microbiome_results")) return(NULL)
    div <- microbiome_results$alpha_diversity; return(div[[variable]][match(samp, div$Sample_ID)])
  }
  NULL
}
rows <- list()
for (i in seq_len(nrow(allo))) {
  v <- get_vals(allo$Variable[i], allo$Category[i]); if (is.null(v)) next
  r <- stat_row(v, allo$Variable[i], allo$Category[i]); if (!is.null(r)) rows[[length(rows) + 1]] <- r
}
genes <- sub("^Gene_", "", roc$Biomarker[grepl("^Gene_", roc$Biomarker)])
gsamp <- intersect(samp, colnames(expr_matrix)); gidx <- match(gsamp, samp)
for (g in genes) {
  if (!(g %in% rownames(expr_matrix))) next
  v <- rep(NA_real_, length(samp)); v[gidx] <- as.numeric(expr_matrix[g, gsamp])
  r <- stat_row(v, paste0("Gene_", g), "Gene (RNA-seq)"); if (!is.null(r)) rows[[length(rows) + 1]] <- r
}
merged <- bind_rows(rows)
merged$Adjusted_P <- p.adjust(merged$P_value, method = "BH")
merged <- merged %>% arrange(desc(AUC))
out <- data.frame(
  Variable = merged$Variable, Category = merged$Category,
  `Median Infection (+)` = signif(merged$Median_pos, 4), `Median Infection (-)` = signif(merged$Median_neg, 4),
  `N (+)` = merged$N_pos, `N (-)` = merged$N_neg, AUC = round(merged$AUC, 3),
  `95% CI lower` = round(merged$CI_lo, 3), `95% CI upper` = round(merged$CI_hi, 3),
  `Cliff's delta` = round(merged$Cliff_delta, 3), `P value` = signif(merged$P_value, 3),
  `Adjusted P value` = signif(merged$Adjusted_P, 3), Direction = merged$Direction, check.names = FALSE)
write.csv(out, file.path(TBL, "Merged_Infection_Data3.csv"), row.names = FALSE)
cat(sprintf("Saved Supplementary Data 3: %d biomarkers (genes = %d)\n",
            nrow(out), sum(merged$Category == "Gene (RNA-seq)")))


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Figure 5 "Infection prediction (RA)" ROC panel: the integrated factor curve
# is the RA-only MOFA Factor 1 (ra_factors[,1]). Respiratory infection is an
# outcome defined only within RA, so the integrated axis comes from the RA-only
# MOFA, not the full-cohort model. Single best biomarker vs the RA-only
# Factor 1, sharing theme/legend/size with the CT-progression ROC.
#
# ROC method: roc(outcome, predictor, direction="auto") (median-based; AUC can
# be <0.5). AUCs are computed, not hardcoded.
# Output: results/panels/Fig5f_infection_ROC_RAfactor.{png,pdf} (80x75mm)
# ============================================================================
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


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Supplementary Fig. 2 "infection x ILD-stratified" panel: BAL Th17.1 and
# GSVA OXPHOS by prospective respiratory-infection status, stratified by
# baseline ILD status. The infection biomarkers are ILD-independent
# (ILD-: Th17.1 p=0.014, OXPHOS p=0.007; ILD+: NS).
#
# Output: results/panels/SFig3c_infection_ILD_stratified.{png,pdf}
# ============================================================================
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
EXCL <- character(0)
load("results/RA_ILD_Workspace.RData")
gsva <- get("gsva_scores")
panel_dir <- "results/panels"

md <- master_data[!master_data$Sample_ID %in% EXCL, ]
ra <- md[!grepl("^Sarcoidosis", md$Sample_ID), ]                       # RA only, n=24
ra$ILD <- factor(ifelse(ra$ILD_Score > 0, "ILD+", "ILD−"), levels = c("ILD−", "ILD+"))
ra$Infection <- factor(ifelse(ra$respiratory_infection == 1, "+", "−"), levels = c("−", "+"))
ra$`BAL Th17.1 (%)` <- as.numeric(ra$BALF_Th17.1_CCR6pos_CXCR3pos)
ra$`GSVA OXPHOS`    <- as.numeric(gsva["OXPHOS", match(ra$Sample_ID, colnames(gsva))])

long <- ra %>%
  select(Sample_ID, ILD, Infection, `BAL Th17.1 (%)`, `GSVA OXPHOS`) %>%
  pivot_longer(c(`BAL Th17.1 (%)`, `GSVA OXPHOS`), names_to = "Variable", values_to = "Value") %>%
  filter(!is.na(Value) & !is.na(Infection))
long$Variable <- factor(long$Variable, levels = c("BAL Th17.1 (%)", "GSVA OXPHOS"))

# stratified Wilcoxon p per (Variable, ILD)
pdat <- long %>% group_by(Variable, ILD) %>%
  summarise(p = tryCatch(suppressWarnings(
              wilcox.test(Value ~ Infection, exact = TRUE)$p.value), error = function(e) NA),
            ymax = max(Value), .groups = "drop")
pdat$lab <- sprintf("p = %.3f", pdat$p)
print(pdat)

tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.85), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          strip.text = element_text(size = rel(0.9), face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(2, 3, 2, 3), "mm"))
}

p <- ggplot(long, aes(Infection, Value)) +
  geom_boxplot(aes(fill = Infection), outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.7, color = "grey20") +
  facet_grid(Variable ~ ILD, scales = "free_y", switch = "y") +
  geom_text(data = pdat, aes(x = 1.5, y = ymax, label = lab),
            inherit.aes = FALSE, size = 2.5, vjust = -0.2, fontface = "italic") +
  scale_fill_manual(values = c("−" = "#4C78A8", "+" = "#C0392B")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = "Respiratory infection", y = NULL) +
  tn(8) + theme(strip.placement = "outside")

ggsave(file.path(panel_dir, "SFig3c_infection_ILD_stratified.pdf"), p,
       width = 95, height = 85, units = "mm", dpi = 300, device = cairo_pdf)
ggsave(file.path(panel_dir, "SFig3c_infection_ILD_stratified.png"), p,
       width = 95, height = 85, units = "mm", dpi = 300, bg = "white")
cat("Saved SFig3c_infection_ILD_stratified (95x85mm)\n")


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Figure 5 panels with descriptive "Factor 1 / Factor 2" labels, and Fig 5b
# showing Factor 1 only (RA vs Sarcoidosis).
#   5a R2 heatmap  : Factor labels
#   5b boxplot     : Factor 1 only (RA vs Sarcoidosis)
#   5c disease ROC : Factor 1 vs best serum biomarker
# Outputs: Fig5a_R2_noacronym, Fig5b_Factor1_only, Fig5c_diseaseROC_noacronym
# ============================================================================
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(pROC); library(patchwork) })
load("results/RA_ILD_Workspace.RData")
e <- new.env(); load("results/MOFA2_6views_Results.RData", envir = e)
factors <- get("factors", e); r2 <- get("r2", e)
panel_dir <- "results/panels"

view_labels <- c(Expression="Expression", BALF_Cytokine="BAL Cytokine",
                 Serum_Cytokine="Serum Cytokine", BALF_FCM="BAL FCM",
                 PB_FCM="PB FCM", Microbiome="Microbiome")
view_cols <- c(Expression="#E67E22", BALF_Cytokine="#3498DB", Serum_Cytokine="#E74C3C",
               BALF_FCM="#27AE60", PB_FCM="#9B59B6", Microbiome="#95A5A6")

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
sv <- function(p, name, w, h) {
  ggsave(file.path(panel_dir, paste0(name,".pdf")), p, width=w, height=h, units="mm", dpi=300, device=cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name,".png")), p, width=w, height=h, units="mm", dpi=300, bg="white")
  cat("  saved", name, "\n")
}
group_m <- ifelse(grepl("^Sarcoidosis", rownames(factors)), "Sarcoidosis", "RA")

# ---------- 5a: R2 heatmap with Factor 1..5 (no IDD/IIV) ----------
fnames <- c(Factor1="Factor 1", Factor2="Factor 2", Factor3="Factor 3",
            Factor4="Factor 4", Factor5="Factor 5")
r2m <- r2$r2_per_factor[[1]]
r2l <- as.data.frame(r2m) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to="View", values_to="R2")
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2m)))
r2l$View <- factor(r2l$View, levels = colnames(r2m))
r2l$FLabel <- factor(fnames[as.character(r2l$Factor)],
                     levels = rev(fnames[rownames(r2m)]))
p_hm <- ggplot(r2l, aes(View, FLabel, fill = R2)) +
  geom_tile(color="white", linewidth=0.8) +
  scale_fill_gradient(low="white", high="#2471A3", name="R² (%)") +
  geom_text(aes(label=ifelse(R2>0.5, sprintf("%.1f",R2), ""),
                color=ifelse(R2>25,"white","black")), size=2.5, fontface="bold", show.legend=FALSE) +
  scale_color_identity() + scale_x_discrete(labels=view_labels) +
  labs(x="", y="") +
  tn(8) + theme(axis.text.x=element_text(angle=40, hjust=1, size=7,
                  color=view_cols[levels(r2l$View)]),
                legend.position="left", legend.key.height=unit(10,"mm"))
total_r2 <- data.frame(View=factor(names(r2$r2_total[[1]]), levels=colnames(r2m)),
                       Total=as.numeric(r2$r2_total[[1]]))
p_bar <- ggplot(total_r2, aes(Total, View, fill=View)) +
  geom_col(width=0.65) + scale_fill_manual(values=view_cols) +
  scale_y_discrete(labels=view_labels) + labs(x="Total R² (%)", y="") +
  tn(7) + theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank())
sv(p_hm + p_bar + plot_layout(widths=c(4,1.2)), "Fig5a_R2_noacronym", 130, 75)

# ---------- 5b: Factor 1 only, RA vs Sarcoidosis ----------
pval <- suppressWarnings(wilcox.test(factors[group_m=="RA",1], factors[group_m=="Sarcoidosis",1],
                                     exact=TRUE)$p.value)
plab <- if (pval < 0.0001) "p < 0.0001" else sprintf("p = %.4f", pval)
cat(sprintf("  Factor 1 disease wilcox exact p = %.6f -> %s\n", pval, plab))
df5b <- data.frame(Group=factor(group_m, levels=c("RA","Sarcoidosis")), Value=factors[,1])
sv(ggplot(df5b, aes(Group, Value, fill=Group)) +
     geom_boxplot(outlier.shape=NA, width=0.6, linewidth=0.3) +
     geom_jitter(width=0.15, size=0.9, alpha=0.6, color="grey25") +
     annotate("text", x=1.5, y=max(df5b$Value), label=plab, size=2.8, fontface="italic", vjust=-0.2) +
     scale_fill_manual(values=c(RA="#C0392B", Sarcoidosis="#2980B9")) +
     scale_y_continuous(expand=expansion(mult=c(0.05,0.15))) +
     labs(x="", y="Factor 1 value") + ggtitle("Factor 1") +
     tn(8) + theme(legend.position="none"),
   "Fig5b_Factor1_only", 62, 75)

# ---------- 5c: disease ROC with "Factor 1" (no IDD) ----------
y_disease <- ifelse(group_m=="RA",1,0)
roc_f1 <- roc(y_disease, factors[,1], quiet=TRUE, direction="auto")
serum <- master_data[match(rownames(factors), master_data$Sample_ID), ]
best<-NULL; ba<-0
for(cc in grep("^Serum_", colnames(serum), value=TRUE)){
  v<-as.numeric(serum[[cc]])
  if(sum(!is.na(v))>=30){ r<-tryCatch(roc(y_disease,v,quiet=TRUE,direction="auto"),error=function(x)NULL)
    if(!is.null(r)&&as.numeric(auc(r))>ba){ba<-as.numeric(auc(r));best<-cc}}}
roc_best<-roc(y_disease, as.numeric(serum[[best]]), quiet=TRUE, direction="auto")
cf<-gsub("^Serum_","",best); cf<-gsub("^GCSF$","G-CSF",cf); cf<-gsub("^IL(\\d)","IL-\\1",cf)
roc_df<-rbind(
  data.frame(Sens=roc_f1$sensitivities, Spec=1-roc_f1$specificities,
             Model=sprintf("Factor 1 (AUC = %.3f)", auc(roc_f1))),
  data.frame(Sens=roc_best$sensitivities, Spec=1-roc_best$specificities,
             Model=sprintf("%s (AUC = %.3f)", cf, auc(roc_best))))
cat(sprintf("  Factor 1 disease AUC=%.3f, %s AUC=%.3f\n", auc(roc_f1), cf, auc(roc_best)))
sv(ggplot(roc_df, aes(Spec, Sens, color=Model)) +
     geom_step(linewidth=0.6, direction="vh") +
     geom_abline(slope=1, intercept=0, linetype="dotted", color="grey60", linewidth=0.3) +
     annotate("segment", x=0, xend=1, y=1, yend=1, linetype="dashed", color="#27AE60", linewidth=0.4) +
     annotate("text", x=0.98, y=0.05, size=2.3, hjust=1, color="#27AE60",
              label="Nested LOOCV (AUC = 0.780)", fontface="italic") +
     scale_color_manual(values=c("#E67E22","#C0392B")) +
     labs(x="1 – Specificity", y="Sensitivity", color="") + ggtitle("RA vs Sarcoidosis") +
     tn(8) + theme(legend.position=c(0.62,0.22),
                   legend.background=element_rect(fill=alpha("white",0.95), color=NA),
                   legend.text=element_text(size=6.5), legend.key.height=unit(3,"mm")),
   "Fig5c_diseaseROC_noacronym", 80, 75)


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Regenerate SuppFig3 so the Greek-alpha panel titles (GRO-α, IFN-α2) render in
# the SAME bold style as the other titles. Helvetica-bold lacks α (font
# fallback made GRO-α/IFN-α2 look different); render those two titles as bold
# plotmath expressions so the α is a proper bold symbol.
# Healthy comparison remains descriptive (green 95% CI band); the in-panel
# p-value is the RA-vs-Sarcoidosis Wilcoxon test (the disease comparison).
# Overwrites results/panels/SFig2_healthy_serum_95CI.{png,pdf}
# ============================================================================
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(patchwork) })
load("results/RA_ILD_Workspace.RData")
s3_table <- read.csv("results/tables/S3_Healthy_Serum_Comparison.csv")
panel_dir <- "results/panels"
col_ra <- "#C0392B"; col_sarc <- "#3498DB"; col_healthy <- "#27AE60"

tn <- function(bs = 7) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.85), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.15), face = "bold", hjust = 0.5,
                                    margin = ggplot2::margin(0,0,2,0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3),
          axis.line = element_blank())
}
fmt_p <- function(p) if (p < 0.0001) "p < 0.0001" else sprintf("p = %.4f", p)

# label / title: use a bold plotmath expression for the alpha-containing names
cyto_info <- list(
  list(col="Serum_GCSF",  ttl="G-CSF"),
  list(col="Serum_IL6",   ttl="IL-6"),
  list(col="Serum_IL13",  ttl="IL-13"),
  list(col="Serum_IL17A", ttl="IL-17A"),
  list(col="Serum_IFNa2", ttl=bquote(bold("IFN-" * alpha * "2"))),
  list(col="Serum_IL22",  ttl="IL-22"),
  list(col="Serum_IL4",   ttl="IL-4"),
  list(col="Serum_IL2",   ttl="IL-2"),
  list(col="Serum_GROa",  ttl=bquote(bold("GRO-" * alpha))),
  list(col="Serum_FGF2",  ttl="FGF-2"))

md <- master_data; plots_s3 <- list()
for (ci in seq_along(cyto_info)) {
  info <- cyto_info[[ci]]
  vals <- as.numeric(md[[info$col]]); grp <- ifelse(md$Sample_Group=="RA","RA","Sarcoidosis")
  ok <- !is.na(vals)
  df <- data.frame(Value=vals[ok], Group=factor(grp[ok], levels=c("RA","Sarcoidosis")))
  wt <- suppressWarnings(wilcox.test(Value ~ Group, data=df, exact=TRUE))
  p_lab <- fmt_p(wt$p.value)
  cyto_short <- gsub("^Serum_","",info$col); row_idx <- which(s3_table$Cytokine==cyto_short)
  h_med <- h_lo <- h_hi <- NA
  if (length(row_idx)==1) {
    hm <- as.numeric(strsplit(gsub("[][,]","", s3_table$Healthy_median_95CI[row_idx])," +")[[1]])
    h_med <- hm[1]; h_lo <- hm[2]; h_hi <- hm[3]
  }
  ymax <- max(df$Value, h_hi, na.rm=TRUE)*1.25
  p <- ggplot(df, aes(Group, Value, fill=Group)) +
    {if(!is.na(h_lo)&&!is.na(h_hi)) annotate("rect", xmin=-Inf, xmax=Inf, ymin=h_lo, ymax=h_hi, fill=col_healthy, alpha=0.15)} +
    {if(!is.na(h_med)) geom_hline(yintercept=h_med, linetype="dashed", color=col_healthy, linewidth=0.4)} +
    geom_boxplot(outlier.shape=NA, width=0.55, linewidth=0.3, alpha=0.7) +
    geom_jitter(width=0.12, size=0.8, alpha=0.6, show.legend=FALSE) +
    scale_fill_manual(values=c(RA=col_ra, Sarcoidosis=col_sarc)) +
    labs(x="", y="pg/mL", title=info$ttl) +
    annotate("text", x=1.5, y=ymax, label=p_lab, size=2.2, fontface="italic") +
    coord_cartesian(ylim=c(0, ymax*1.05), clip="off") +
    tn(7) + theme(legend.position="none", axis.text.x=element_text(size=6))
  plots_s3[[ci]] <- p
}
p_s3 <- wrap_plots(plots_s3, ncol=5) + plot_annotation(theme=theme(plot.margin=ggplot2::margin(2,2,2,2)))
ggsave(file.path(panel_dir,"SFig2_healthy_serum_95CI.pdf"), p_s3, width=250, height=110, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(panel_dir,"SFig2_healthy_serum_95CI.png"), p_s3, width=250, height=110, units="mm", dpi=300, bg="white")
cat("Saved SFig2_healthy_serum_95CI (alpha titles as bold plotmath)\n")


# ============================================================================
# Final figure/output step
# ============================================================================
# Remove the green "Nested LOOCV (AUC = 0.962)" dashed line + text from Fig5c
# (redundant with Fig2l; overlapped the legend). CUD colors retained.
setwd(BASEDIR)
suppressPackageStartupMessages({library(tidyverse);library(pROC)})
CUD_ORANGE<-"#FF9900"; CUD_BLUE<-"#0041FF"; panel_dir<-"results/panels"
load("results/RA_ILD_Workspace.RData"); master_all<-master_data
tn<-function(bs=8){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.05),face="bold",hjust=0,margin=ggplot2::margin(0,0,3,0)),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,3,2,3),"mm"))}
add_ep<-function(spec,sens){ord<-order(spec,sens);spec<-spec[ord];sens<-sens[ord]
  if(spec[1]!=0||sens[1]!=0){spec<-c(0,spec);sens<-c(0,sens)}
  if(tail(spec,1)!=1||tail(sens,1)!=1){spec<-c(spec,1);sens<-c(sens,1)};list(spec=spec,sens=sens)}
e6<-new.env(); load("results/MOFA2_6views_Results.RData",envir=e6); factors6<-get("factors",e6)
y<-ifelse(grepl("^Sarcoidosis",rownames(factors6)),0,1)
roc_f1<-roc(y,factors6[,1],quiet=TRUE,direction="auto")
serum<-master_all[match(rownames(factors6),master_all$Sample_ID),]
best<-NULL;ba<-0
for(cc in grep("^Serum_",colnames(serum),value=TRUE)){v<-as.numeric(serum[[cc]])
  if(sum(!is.na(v))>=30){r<-tryCatch(roc(y,v,quiet=TRUE,direction="auto"),error=function(x)NULL)
    if(!is.null(r)&&as.numeric(auc(r))>ba){ba<-as.numeric(auc(r));best<-cc}}}
roc_best<-roc(y,as.numeric(serum[[best]]),quiet=TRUE,direction="auto")
cf<-gsub("^Serum_","",best);cf<-gsub("^GCSF$","G-CSF",cf);cf<-gsub("^IL(\\d)","IL-\\1",cf)
ep1<-add_ep(1-roc_f1$specificities,roc_f1$sensitivities);ep2<-add_ep(1-roc_best$specificities,roc_best$sensitivities)
df<-rbind(data.frame(Spec=ep1$spec,Sens=ep1$sens,Model=sprintf("Factor 1 (AUC = %.3f)",auc(roc_f1))),
          data.frame(Spec=ep2$spec,Sens=ep2$sens,Model=sprintf("%s (AUC = %.3f)",cf,auc(roc_best))))
df$Model<-factor(df$Model,levels=unique(df$Model))
p<-ggplot(df,aes(Spec,Sens,color=Model))+geom_step(linewidth=0.6,direction="vh")+
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey60",linewidth=0.3)+
  scale_color_manual(values=c(CUD_BLUE,CUD_ORANGE))+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="on")+
  labs(x="1 – Specificity",y="Sensitivity",color="")+ggtitle("RA vs Sarcoidosis")+
  tn(8)+theme(legend.position=c(0.62,0.18),legend.background=element_rect(fill=alpha("white",0.95),color=NA),
              legend.text=element_text(size=6.5),legend.key.height=unit(3,"mm"))
ggsave(file.path(panel_dir,"Fig5c_diseaseROC_noacronym.pdf"),p,width=80,height=75,units="mm",dpi=300,device=cairo_pdf)
ggsave(file.path(panel_dir,"Fig5c_diseaseROC_noacronym.png"),p,width=80,height=75,units="mm",dpi=300,bg="white")
cat(sprintf("5c regenerated WITHOUT LOOCV annot: F1=%.3f, %s=%.3f\n",auc(roc_f1),cf,auc(roc_best)))


# ============================================================================
# Final figure/output step
# ============================================================================
# Medication-confounding for the RA-specific outcome predictors used in Fig 5:
# RA-only MOFA Factor 1 (ra_factors[,1]) and BAL Th17.1, vs GC/MTX/NSAID.
# Writes Medication_Confounding_RAfactor.csv and regenerates SFig4_medication
# (2 rows x 3 medications).
suppressPackageStartupMessages({library(tidyverse);library(ggplot2);library(patchwork)})
setwd(BASEDIR); panel_dir<-"results/panels"
EXCL<-character(0)
tn<-function(bs=7){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.0),face="bold",hjust=0.5,margin=ggplot2::margin(0,0,2,0)),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,2,2,2),"mm"))}
sv<-function(p,name,w,h){ggsave(file.path(panel_dir,paste0(name,".pdf")),p,width=w,height=h,units="mm",dpi=300,device=cairo_pdf)
  ggsave(file.path(panel_dir,paste0(name,".png")),p,width=w,height=h,units="mm",dpi=300,bg="white");cat("Saved",name,"\n")}
fmt_p<-function(p){if(p<0.001)return("p < 0.001");sprintf("p = %.2f",p)}
col_ra<-"#C0392B"; col_sarc<-"#2980B9"
load("results/RA_ILD_Workspace.RData")
e<-new.env(); load("results/MOFA2_RA_only_Results.RData",envir=e); ra_factors<-get("ra_factors",e)
ct<-read.csv("results/tables/master_data_with_CT_all.csv"); ct<-ct[!duplicated(ct$Sample_ID),]
ids<-rownames(ra_factors)
md<-master_data[!master_data$Sample_ID %in% EXCL,]
vars<-list(list(nm="RA Factor 1", lab="RA Factor 1", v=as.numeric(ra_factors[,1])),
           list(nm="BAL Th17.1", lab="BAL Th17.1 (%)", v=as.numeric(md$BALF_Th17.1_CCR6pos_CXCR3pos[match(ids,md$Sample_ID)])))
med_cols<-c("steroid.1.0.","MTX","NSAIDS"); med_lab<-c("GC","MTX","NSAIDs")
rows<-list(); plots<-list(); pi<-1
for(vv in vars){ for(k in seq_along(med_cols)){
  mv<-ct[[med_cols[k]]][match(ids,ct$Sample_ID)]; ok<-!is.na(vv$v)&!is.na(mv)
  on<-vv$v[ok&mv==1]; off<-vv$v[ok&mv==0]
  w<-suppressWarnings(wilcox.test(on,off,exact=TRUE))
  rows[[length(rows)+1]]<-data.frame(Medication=med_lab[k],Variable=vv$nm,N_on=length(on),N_off=length(off),
                                     Median_on=round(median(on),2),Median_off=round(median(off),2),P=round(w$p.value,4))
  df<-data.frame(Value=c(on,off),Med=factor(c(rep(paste0(med_lab[k]," (+)"),length(on)),rep(paste0(med_lab[k]," (−)"),length(off))),
                 levels=c(paste0(med_lab[k]," (−)"),paste0(med_lab[k]," (+)"))))
  yr<-range(df$Value); ymax<-yr[2]+diff(yr)*0.25
  plots[[pi]]<-ggplot(df,aes(Med,Value,fill=Med))+geom_boxplot(outlier.shape=NA,width=0.55,linewidth=0.3,alpha=0.7)+
    geom_jitter(width=0.12,size=0.8,alpha=0.6,show.legend=FALSE)+scale_fill_manual(values=c(col_sarc,col_ra))+
    labs(x="",y=vv$lab,title=med_lab[k])+annotate("text",x=1.5,y=ymax,label=fmt_p(w$p.value),size=2.2,fontface="italic")+
    coord_cartesian(ylim=c(min(yr[1],0),ymax*1.05),clip="off")+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=6))
  pi<-pi+1 }}
res<-do.call(rbind,rows); write.csv(res,"results/tables/Medication_Confounding_RAfactor.csv",row.names=FALSE)
print(res)
p_s4<-wrap_plots(plots,ncol=3)+plot_annotation(theme=theme(plot.margin=ggplot2::margin(2,2,2,2)))
sv(p_s4,"SFig4_medication",160,100)


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Build the 9-sheet Supplementary Tables workbook in body first-mention order.
# ST2 (healthy serum) reports only the KW omnibus and the RA-vs-Sarcoidosis
# test; healthy-comparison hypothesis tests are omitted per study protocol.
#   ST1 = clinical RA vs sarcoidosis
#   ST2 = serum cytokine RA/sarcoidosis/healthy
#   ST3 = clinical by infection status
#   ST4 = single biomarker ROC (infection)
#   ST5 = BAL microbiome by infection
#   ST6 = clinical vs multi-omics
#   ST7 = BAL/serum profiles by ILD
#   ST8 = CT progression correlations
#   ST9 = medication confounding
# Output: results/tables/Supplementary_Tables.xlsx (overwrite)
# ============================================================================
setwd(BASEDIR)
suppressPackageStartupMessages(library(openxlsx))
TBL <- "results/tables"
SRC_XLSX <- "results/tables/Supplementary_Tables.xlsx"
OUT <- "results/tables/Supplementary_Tables.xlsx"
MINUS<-"−"; SUB2<-"₂"; RHO<-"ρ"; ENDASH<-"–"
rnd<-function(x,d) ifelse(is.na(x),NA,round(as.numeric(x),d))
sgf<-function(x,d=3) ifelse(is.na(x),NA,signif(as.numeric(x),d))

title_style  <- createStyle(fontSize=11, textDecoration="bold")
header_style <- createStyle(fontSize=10, textDecoration="bold", border="TopBottom",
                            borderStyle="thin", fgFill="#F2F2F2", halign="left",
                            valign="center", wrapText=TRUE)
data_style   <- createStyle(fontSize=10, halign="left", valign="top")
section_style<- createStyle(fontSize=10, textDecoration="bold")
wb <- createWorkbook()
add_sheet <- function(sheet,title,headers,dat){
  addWorksheet(wb,sheet)
  writeData(wb,sheet,title,startRow=1,startCol=1); addStyle(wb,sheet,title_style,rows=1,cols=1)
  colnames(dat)<-headers
  writeData(wb,sheet,dat,startRow=2,startCol=1,headerStyle=header_style)
  if(nrow(dat)>0) addStyle(wb,sheet,data_style,rows=3:(2+nrow(dat)),cols=1:ncol(dat),gridExpand=TRUE,stack=TRUE)
  setColWidths(wb,sheet,cols=1:ncol(dat),widths="auto")
  cat(sprintf("  %s: %d x %d  | %s\n",sheet,nrow(dat),ncol(dat),substr(title,1,45)))
}

## ST1 = clinical RA vs sarcoidosis (copied verbatim from the source workbook)
# readWorkbook can return numeric character refs as literal "&#NNNN;" strings; decode them.
dehtml<-function(x){
  if(is.na(x)) return(x)
  repeat{ g<-regmatches(x,regexpr("&#[0-9]+;",x)); if(length(g)==0) break
    x<-sub("&#[0-9]+;", intToUtf8(as.integer(sub("&#([0-9]+);","\\1",g))), x) }
  x<-gsub("&gt;",">",x,fixed=TRUE); x<-gsub("&lt;","<",x,fixed=TRUE); x<-gsub("&amp;","&",x,fixed=TRUE); x }
wb_old<-loadWorkbook(SRC_XLSX); st1<-readWorkbook(wb_old,sheet="ST1",colNames=FALSE,skipEmptyRows=FALSE)[,1:3]
st1[]<-lapply(st1,function(col) vapply(col,dehtml,character(1)))
st1[]<-lapply(st1,function(col) gsub("†","",col))   # drop dagger footnote marker -> plain N/A
addWorksheet(wb,"ST1"); writeData(wb,"ST1",st1,startRow=1,startCol=1,colNames=FALSE)
addStyle(wb,"ST1",title_style,rows=1,cols=1); addStyle(wb,"ST1",header_style,rows=2,cols=1:3,gridExpand=TRUE)
for(r in 3:nrow(st1)){ c2<-st1[r,2];c3<-st1[r,3]
  sec<-!is.na(st1[r,1])&&(is.na(c2)||c2=="")&&(is.na(c3)||c3=="")
  addStyle(wb,"ST1",if(sec)section_style else data_style,rows=r,cols=1:3,gridExpand=TRUE,stack=TRUE)}
setColWidths(wb,"ST1",cols=1:3,widths="auto"); cat("  ST1: verbatim\n")

## ST2 = serum cytokine; drop RA-vs-Healthy and Sarc-vs-Healthy P
d<-read.csv(file.path(TBL,"S3_Healthy_Serum_Comparison.csv"),check.names=FALSE)
d[["RA_vs_Sarc_p"]]<-sgf(d[["RA_vs_Sarc_p"]],3)
add_sheet("ST2","Supplementary Table 2. Serum cytokine comparison: RA, sarcoidosis, and healthy controls.",
  c("Cytokine","Healthy\nN","Healthy\nmedian [95% CI]","Sarcoidosis\nN","Sarcoidosis\nmedian [IQR]",
    "RA\nN","RA\nmedian [IQR]","RA vs Sarc\nP"),
  d[,c("Cytokine","Healthy_n","Healthy_median_95CI","Sarcoidosis_n","Sarcoidosis_median_IQR","RA_n","RA_median_IQR","RA_vs_Sarc_p")])

## ST3 = clinical by infection status
d<-read.csv(file.path(TBL,"Infection_Clinical_Characteristics.csv"),check.names=FALSE,colClasses="character")
add_sheet("ST3","Supplementary Table 3. Clinical characteristics by prospective infection status (RA, n = 24).",
  c("Variable","Infection (+)\n(n = 6)",paste0("Infection (",MINUS,")\n(n = 18)"),"P value"),
  d[,c("Variable","InfPos","InfNeg","P")])

## ST4 = single biomarker ROC, sorted by descending AUC
d<-read.csv(file.path(TBL,"Infection_SingleBiomarker_ROC.csv"),check.names=FALSE); d<-d[order(-d$AUC),]
d$AUC<-rnd(d$AUC,3); d$CI_lo<-rnd(d$CI_lo,3); d$CI_hi<-rnd(d$CI_hi,3)
add_sheet("ST4","Single-biomarker ROC analysis for infection prediction (source data for Figure 5).",
  c("Biomarker","AUC","95% CI\nlower","95% CI\nupper","Direction","N"),
  d[,c("Biomarker","AUC","CI_lo","CI_hi","Direction","N")])

## NOTE: in the final submission the sheets below are split across files --
## ST5 -> Supplementary Data 4 (microbiome); ST6 -> Supplementary Table 4 (clinical vs multi-omics);
## ST7 -> Supplementary Data 5 (ILD profiles); ST8 -> Supplementary Data 6 (CT correlations);
## ST9 -> Supplementary Table 5 (medication). ST4 is source data for Figure 5 (single-biomarker ROC).
## ST1-ST3 = Supplementary Tables 1-3. Sheet tab names (ST*) are kept for internal reference.
## ST5 = microbiome by infection -> Supplementary Data 4
d<-read.csv(file.path(TBL,"Microbiome_Genus_Infection.csv"),check.names=FALSE); d<-d[order(d$P_value),]
d$Mean_InfPos<-rnd(d$Mean_InfPos,4); d$Mean_InfNeg<-rnd(d$Mean_InfNeg,4); d$Log2FC<-rnd(d$Log2FC,2)
d$P_value<-sgf(d$P_value,3); d$P_adjusted<-sgf(d$P_adjusted,3)
add_sheet("ST5","Supplementary Data 4. BAL microbiome genus-level comparison by prospective infection status (RA, n = 24).",
  c("Genus","Mean (Infection +)",paste0("Mean (Infection ",MINUS,")"),paste0("Log",SUB2,"FC"),"P value","FDR-adjusted P"),
  d[,c("Genus","Mean_InfPos","Mean_InfNeg","Log2FC","P_value","P_adjusted")])

## ST6 = clinical markers vs multi-omics
d<-read.csv(file.path(TBL,"Clinical_vs_Multiomics_Comparison.csv"),check.names=FALSE)
d$Inf_AUC<-rnd(d$Inf_AUC,3); d$Inf_P<-rnd(d$Inf_P,3); d$Prog_AUC<-rnd(d$Prog_AUC,3); d$Prog_P<-rnd(d$Prog_P,3)
add_sheet("ST6","Supplementary Table 4. Clinical markers versus multi-omics biomarkers (exploratory).",
  c("Category","Variable","Infection\nAUC","Infection\nP value","Progression\nAUC","Progression\nP value"),
  d[,c("Category","Variable","Inf_AUC","Inf_P","Prog_AUC","Prog_P")])

## ST7 = profiles by ILD status, sorted by ascending P
d<-read.csv(file.path(TBL,"ST4_ILD_stratification_updated.csv"),check.names=FALSE); d<-d[order(d$P_value),]
d$Mean_ILD_pos<-rnd(d$Mean_ILD_pos,2); d$Mean_ILD_neg<-rnd(d$Mean_ILD_neg,2)
d$P_value<-sgf(d$P_value,3); d$P_adjusted<-sgf(d$P_adjusted,3)
add_sheet("ST7","Supplementary Data 5. BAL cellular and serum cytokine profiles by ILD status (ILD+, n = 11; ILD−, n = 13).",
  c("Variable","Mean\nILD (+)",paste0("Mean\nILD (",MINUS,")"),"N\nILD (+)",paste0("N\nILD (",MINUS,")"),"P value","FDR-adjusted\nP value"),
  d[,c("Variable","Mean_ILD_pos","Mean_ILD_neg","N_ILD_pos","N_ILD_neg","P_value","P_adjusted")])

## ST8 = CT progression correlations, sorted by ascending p
d<-read.csv(file.path(TBL,"CT_Delta_vs_NonGene_RAonly.csv"),check.names=FALSE); d<-d[order(d$p),]
d$rho<-rnd(d$rho,3); d$p<-sgf(d$p,3)
add_sheet("ST8","Supplementary Data 6. CT progression correlations with non-gene biomarkers (RA only, n = 24).",
  c("CT variable","Biomarker",paste0("Spearman ",RHO),"P value","N"),
  d[,c("CT_Var","Biomarker","rho","p","n")])

## ST9 = medication confounding
d<-read.csv(file.path(TBL,"Medication_Confounding_RAfactor.csv"),check.names=FALSE)
d$Median_on<-rnd(d$Median_on,2); d$Median_off<-rnd(d$Median_off,2); d$P<-rnd(d$P,4)
add_sheet("ST9","Supplementary Table 5. Medication confounding analysis for RA-specific outcome predictors (RA Factor 1 and BAL Th17.1).",
  c("Medication","Variable","N (on)","N (off)","Median (on)","Median (off)","P value"),
  d[,c("Medication","Variable","N_on","N_off","Median_on","Median_off","P")])

saveWorkbook(wb,OUT,overwrite=TRUE); cat(sprintf("\nSaved: %s\n",OUT))


# ============================================================================
# Final figure/output step
# ============================================================================
# ============================================================================
# Recolor the three Figure 5 ROC panels (5c disease, 5f infection, 5i CT
# progression) with Color-Universal-Design (CUD) safe colors so the two curves
# in each panel are easy to tell apart for P/D-type colour vision.
#   single biomarker curve    -> orange #FF9900 (warm)
#   integrative MOFA factor   -> blue   #0041FF (cool)
#   5c nested-LOOCV reference -> green  #35A16B (bluish green, dashed)
# Only the colour values differ from the source panels; all data, logic, AUCs,
# endpoints, titles, and positions are unchanged.
# Overwrites the three panel PNG/PDFs in results/panels/.
# ============================================================================
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
group_m <- ifelse(grepl("^Sarcoidosis", rownames(factors6)), "Sarcoidosis", "RA")
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
           label = "Nested LOOCV (AUC = 0.780)", fontface = "italic") +
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
sv(p5i, "Fig5i_RA_progression_ROC", 80, 75)

cat("\nDone: 3 ROC panels recolored (CUD orange/blue).\n")


# ============================================================================
# Final figure/output step
# ============================================================================
# Fig5h: color points by baseline ILD status (CUD: ILD− blue, ILD+ orange) + legend.
# Single overall regression line + rho/p (overall Spearman) retained.
setwd(BASEDIR)
suppressPackageStartupMessages({library(tidyverse);library(ggplot2)})
EXCL<-character(0); panel_dir<-"results/panels"
CUD_BLUE<-"#0041FF"; CUD_ORANGE<-"#FF9900"
load("results/RA_ILD_Workspace.RData")
e<-new.env(); load("results/MOFA2_RA_only_Results.RData",envir=e); ra_factors<-get("ra_factors",e)
ct<-read.csv("results/tables/master_data_with_CT_all.csv"); ct<-ct[!duplicated(ct$Sample_ID),]
ids<-rownames(ra_factors)
dhl<-as.numeric(ct$CT_Delta_Healthy_Lung_pct[match(ids,ct$Sample_ID)])
ild<-master_data$ILD_Group[match(ids,master_data$Sample_ID)]
ild<-ifelse(ild=="ILD_Positive","ILD+",ifelse(ild=="ILD_Negative","ILD−",NA))
f1<-as.numeric(ra_factors[,1]); ok<-!is.na(f1)&!is.na(dhl)
sp<-suppressWarnings(cor.test(f1[ok],dhl[ok],method="spearman",exact=TRUE))
cat(sprintf("rho=%.3f p=%.4f | ILD+ n=%d ILD- n=%d\n",sp$estimate,sp$p.value,sum(ild[ok]=="ILD+",na.rm=T),sum(ild[ok]=="ILD−",na.rm=T)))
df<-data.frame(Fval=f1[ok],DHL=dhl[ok],ILD=factor(ild[ok],levels=c("ILD−","ILD+")))
tn<-function(bs=8){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        legend.text=element_text(size=rel(0.8)),legend.title=element_text(size=rel(0.85),face="bold"),
        legend.key.size=unit(3,"mm"),legend.background=element_rect(fill=alpha("white",0.95),color=NA),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,3,2,3),"mm"))}
p<-ggplot(df,aes(Fval,DHL))+
  geom_smooth(method="lm",se=TRUE,color="grey30",linewidth=0.4,fill="grey90",formula=y~x)+
  geom_hline(yintercept=0,linetype="dotted",color="grey60",linewidth=0.2)+
  geom_point(aes(color=ILD),size=1.8,alpha=0.85)+
  scale_color_manual(values=c("ILD−"=CUD_BLUE,"ILD+"=CUD_ORANGE),name="ILD status")+
  annotate("text",x=-Inf,y=Inf,hjust=-0.05,vjust=1.3,size=3.3,
           label=paste0("ρ=",sprintf("%.3f",sp$estimate),"\np=",sprintf("%.4f",sp$p.value)))+
  labs(x="RA Factor 1",y="Δ Healthy Lung (%)")+
  tn(8)+theme(legend.position=c(0.83,0.22))
ggsave(file.path(panel_dir,"Fig5h_RA_factor_CT.pdf"),p,width=80,height=70,units="mm",dpi=300,device=cairo_pdf)
ggsave(file.path(panel_dir,"Fig5h_RA_factor_CT.png"),p,width=80,height=70,units="mm",dpi=300,bg="white")
cat("Saved Fig5h_RA_factor_CT (ILD-colored, CUD)\n")


# ============================================================================
# Final figure/output step
# ============================================================================
# Regenerate Supplementary Fig 4 (medication confounding) with the MOFA factor
# axis labels as "Factor 1" / "Factor 2".
suppressPackageStartupMessages({library(tidyverse);library(ggplot2);library(patchwork)})
setwd(BASEDIR); panel_dir<-"results/panels"
EXCL<-character(0)
tn<-function(bs=7){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.0),face="bold",hjust=0.5,margin=ggplot2::margin(0,0,2,0)),
        legend.text=element_text(size=rel(0.8)),legend.key.size=unit(3,"mm"),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,2,2,2),"mm"))}
sv<-function(p,name,w,h){ggsave(file.path(panel_dir,paste0(name,".pdf")),p,width=w,height=h,units="mm",dpi=300,device=cairo_pdf)
  ggsave(file.path(panel_dir,paste0(name,".png")),p,width=w,height=h,units="mm",dpi=300,bg="white");cat("Saved",name,"\n")}
fmt_p<-function(p){if(p<0.001)return("p < 0.001");sprintf("p = %.2f",p)}
col_ra<-"#C0392B"; col_sarc<-"#2980B9"
load("results/RA_ILD_Workspace.RData"); load("results/MOFA2_6views_Results.RData")
ct_data<-read.csv("results/tables/master_data_with_CT_all.csv"); ct_unique<-ct_data[!duplicated(ct_data$Sample_ID),]
ra_samples<-setdiff(rownames(factors)[grepl("^RA[0-9]",rownames(factors))],EXCL)
med_cols<-c("steroid.1.0.","MTX","NSAIDS"); med_labels<-c("GC","MTX","NSAIDs")
factor_names<-c("Factor 1","Factor 2")
plots_s4<-list(); pi<-1
for(fi in 1:2){ fn<-factor_names[fi]; fac_vals<-factors[ra_samples,fi]
  for(mi in seq_along(med_cols)){ mc<-med_cols[mi]; ml<-med_labels[mi]
    med_vals<-ct_unique[[mc]][match(ra_samples,ct_unique$Sample_ID)]; ok<-!is.na(fac_vals)&!is.na(med_vals)
    on_vals<-fac_vals[ok&med_vals==1]; off_vals<-fac_vals[ok&med_vals==0]
    if(length(on_vals)<2||length(off_vals)<2) next
    wt<-wilcox.test(on_vals,off_vals,exact=TRUE); p_lab<-fmt_p(wt$p.value)
    df<-data.frame(Value=c(on_vals,off_vals),
      Med=factor(c(rep(paste0(ml," (+)"),length(on_vals)),rep(paste0(ml," (−)"),length(off_vals))),
                 levels=c(paste0(ml," (−)"),paste0(ml," (+)"))))
    yr<-range(df$Value); ymax<-yr[2]+diff(yr)*0.25
    p<-ggplot(df,aes(Med,Value,fill=Med))+geom_boxplot(outlier.shape=NA,width=0.55,linewidth=0.3,alpha=0.7)+
      geom_jitter(width=0.12,size=0.8,alpha=0.6,show.legend=FALSE)+scale_fill_manual(values=c(col_sarc,col_ra))+
      labs(x="",y=fn,title=ml)+annotate("text",x=1.5,y=ymax,label=p_lab,size=2.2,fontface="italic")+
      coord_cartesian(ylim=c(min(yr[1],0),ymax*1.05),clip="off")+tn(7)+
      theme(legend.position="none",axis.text.x=element_text(size=6))
    plots_s4[[pi]]<-p; pi<-pi+1 }}
p_s4<-wrap_plots(plots_s4,ncol=3)+plot_annotation(theme=theme(plot.margin=ggplot2::margin(2,2,2,2)))
sv(p_s4,"SFig4_medication",160,100)
