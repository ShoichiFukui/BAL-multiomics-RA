# panels_fig1_atlas.R -- Figure 1b-e: scRNA-seq reference atlas (public GEO data): UMAP by dataset (b) and by cell type (c),
#   canonical marker dot plot (d) and cell-type composition of each dataset (e).
# Produces: <ART_FIGS>/n31/panels/fig1/Fig1b-e.png/.pdf.
# Inputs: output_BayesPrism/BAL_reference_author_annotated.rds (01a_reference_atlas.R).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(Seurat) })
panel_dir <- file.path(ART_FIGS, "n31", "panels"); dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
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

CC <- c(Macrophage = "#E74C3C", T_cell = "#3498DB", NK = "#2ECC71",
        B_cell = "#9B59B6", Plasma = "#E67E22", Neutrophil = "#F1C40F",
        DC = "#8B4513", Epithelial = "#E91E63", Mast = "#95A5A6",
        Other = "#BDC3C7")
ref <- readRDS("output_BayesPrism/BAL_reference_author_annotated.rds")
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
# Figure 1b -- UMAP by dataset
col_ds <- c("GSE145926" = "#E74C3C", "GSE193782" = "#3498DB", "GSE184735" = "#27AE60")
sv(ggplot(ds, aes(U1, U2, color = DS)) +
     geom_point(size = 0.06, alpha = 0.2) +
     scale_color_manual(values = col_ds,
                        labels = c("GSE145926" = "COVID-19",
                                   "GSE193782" = "Healthy",
                                   "GSE184735" = "Sarcoidosis")) +
     guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
     labs(x = "UMAP 1", y = "UMAP 2", color = "") +
     tn(8) + theme(legend.position = "bottom"),
   "Fig1b", 85, 75)
# Figure 1c -- UMAP by cell type
sv(ggplot(ds, aes(U1, U2, color = CT)) +
     geom_point(size = 0.06, alpha = 0.2) +
     scale_color_manual(values = CC, breaks = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"), labels = function(x) sub("_", " ", x)) +
     guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                 ncol = 5, byrow = TRUE)) +
     labs(x = "UMAP 1", y = "UMAP 2", color = "") +
     tn(8) + theme(legend.position = "bottom",
                   legend.text = element_text(size = 6),
                   legend.key.size = unit(3, "mm")),
   "Fig1c", 90, 85)
# Figure 1d -- canonical marker dot plot
dd$M <- factor(dd$M, levels = rev(mks)); dd$CT <- factor(dd$CT, levels = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"))
sv(ggplot(dd, aes(CT, M, size = Pct, color = Avg)) +
     geom_point() +
     scale_size_continuous(range = c(0.5, 3.5), name = "%Exp") +
     scale_color_gradient(low = "grey90", high = "#C0392B", name = "Avg") +
     scale_x_discrete(labels = function(x) sub("_", " ", x)) + labs(x = "", y = "") +
     tn(7) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                   axis.text.y = element_text(size = 6, face = "italic"),
                   legend.text = element_text(size = 10.5), legend.title = element_text(size = 11), legend.key.size = unit(4, "mm")),
   "Fig1d", 100, 80)
# Figure 1e -- cell-type composition of each reference dataset
cm <- ds %>% mutate(CT = factor(CT, levels = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"))) %>% group_by(DS, CT) %>% summarise(N = n(), .groups = "drop") %>%
  group_by(DS) %>% mutate(P = N / sum(N))
cm$DS <- recode(cm$DS, "GSE145926" = "COVID-19", "GSE193782" = "Healthy",
                "GSE184735" = "Sarcoidosis")
sv(ggplot(cm, aes(DS, P, fill = CT)) +
     geom_col(width = 0.6) + scale_fill_manual(values = CC, breaks = c("Macrophage", "T_cell", "NK", "B_cell", "Plasma", "Neutrophil", "DC", "Epithelial", "Mast", "Other"), labels = function(x) sub("_", " ", x)) +
     labs(x = "", y = "Proportion", fill = "") +
     tn(8) + theme(legend.position = "right",
                   legend.text = element_text(size = 11),
                   legend.key.size = unit(4.2, "mm")),
   "Fig1e", 110, 75)
