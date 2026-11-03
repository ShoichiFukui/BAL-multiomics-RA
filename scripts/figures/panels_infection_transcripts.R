# panels_infection_transcripts.R -- Figure 3g (surfactant / alveolar-epithelial transcripts, infection + vs -) and Supplementary Fig. 3a
#   (GO BP GSEA, infection + vs -), RA n = 24.
# Produces: <ART_FIGS>/n31/panels/infection/Fig3g.png/.pdf and SFig3a.png/.pdf.
# Inputs: DEG_Infection_pos_vs_neg.csv and GSEA_Infection_GO_BP.csv (03_DEG_infection_within_RA.R).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(tidyverse); library(forcats) })
TB <- file.path(ART_TABLES, "n31"); panel_dir <- file.path(ART_FIGS, "n31", "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
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
}

# Figure 3g -- the nine surfactant / alveolar-epithelial transcripts (DESeq2, genome-wide BH)
deg_inf <- read.csv(file.path(TB, "DEG_Infection_pos_vs_neg.csv"))
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
     labs(x = "log2FC (Infection + vs −)", y = "", fill = "") +
     tn(8) + theme(legend.position = "bottom", legend.direction = "horizontal",
                   legend.text = element_text(size = 9),
                   legend.key.size = unit(3.5, "mm"),
                   axis.text.y = element_text(face = "italic")),
   "Fig3g", 90, 90)

# Supplementary Fig. 3a -- 13 GO BP sets share the lowest adjusted P (P at the fgsea floor of 1e-10); the eight with the largest |NES| are shown
ginf <- read.csv(file.path(TB, "GSEA_Infection_GO_BP.csv"))
stopifnot(sum(ginf$p.adjust == min(ginf$p.adjust)) == 13)
tgi <- head(ginf[order(ginf$p.adjust, -abs(ginf$NES)), ], 8) %>%
  mutate(D = ifelse(nchar(Description) > 42, str_wrap(Description, 36), Description))
tgi <- tgi[order(tgi$NES), ]; tgi$D <- fct_inorder(tgi$D)
sv(ggplot(tgi, aes(NES, D, fill = NES > 0)) +
     geom_col(width = 0.55, show.legend = FALSE) +
     scale_fill_manual(values = c("TRUE" = "#C0392B",
                                  "FALSE" = "#2980B9")) +
     geom_vline(xintercept = 0, linewidth = 0.2) +
     labs(x = "NES", y = "") +
     tn(7) + theme(axis.text.y = element_text(size = 6)),
   "SFig3a", 100, 70)
