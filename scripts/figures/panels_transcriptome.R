# panels_transcriptome.R -- Figure 2a-d and 2g panels.
# Produces: panels/Fig2*.png and .pdf.
# Inputs: 02 outputs, results/RA_ILD_Workspace.RData, scripts/figures/panel_map.json (shipped).
suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel); library(jsonlite)
})
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
cfg <- set_cohort("n31")

TB  <- file.path(cfg$tables)
tb  <- function(f) read.csv(file.path(TB, paste0("transcriptome_", f)),
                            check.names = FALSE)
panel_dir <- file.path(ART_FIGS, "n31", "panels"); dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
pmap <- fromJSON("scripts/figures/panel_map.json")   # pixel sizes of the figure panels (shipped with the code)

# ---- theme (verbatim from the figure code of the article) ------------------
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
CS <- c("Sarcoidosis" = "#2980B9", "RA-nonILD" = "#27AE60", "RA-ILD" = "#C0392B")

sv <- function(p, name, size_of = name) {
  px <- pmap[[size_of]]$px
  ggsave(file.path(panel_dir, paste0(name, ".png")), p,
         width = px[1], height = px[2], units = "px", dpi = 300, bg = "white")
  ggsave(file.path(panel_dir, paste0(name, ".pdf")), p,
         width = px[1] / 300 * 25.4, height = px[2] / 300 * 25.4, units = "mm",
         device = cairo_pdf)
  
}

# ---- data ---------------------------------------------------------------------
pcs  <- tb("PCA_scores.csv")
pve  <- tb("PCA_variance.csv")
deg  <- tb("DEG_ShrunkLFC.csv")
ggo  <- tb("GSEA_GO_BP.csv")
gkk  <- tb("GSEA_KEGG.csv")
igt  <- tb("IGH_genes.csv")
met  <- tb("metrics.csv")
check_n(pcs$Sample_ID, cfg, "PCA_scores")

e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
meta <- e$meta; rm(e); invisible(gc())

# ---- GSEA term selection ------------------------------------------------------
sel_orig <- function(g) head(g[order(g$p.adjust), ], 10)
dedup <- function(x) { x <- x[order(x$p.adjust, -abs(x$NES)), ]; x[!duplicated(paste(x$setSize, x$core_enrichment)), ] }
sel_terms <- function(g, k = 5) {
  g <- g[g$p.adjust < 0.05, ]; g <- rbind(dedup(g[g$NES > 0, ]), dedup(g[g$NES < 0, ]))
  pick <- function(x) head(x[order(x$p.adjust, -abs(x$NES)), ], k)
  rbind(pick(g[g$NES > 0, ]), pick(g[g$NES < 0, ]))
}

pd <- pcs %>% transmute(PC1, PC2, SID = Sample_ID,
  Sub = recode(as.character(meta$Subgroup[match(SID, meta$Sample_ID)]),
               "Control" = "Sarcoidosis", "RA_nonILD" = "RA-nonILD", "RA_ILD" = "RA-ILD"))
stopifnot(!any(is.na(pd$Sub)), all((pd$Sub == "Sarcoidosis") == is_sarc(pd$SID)))
ve <- setNames(pve$Pct_variance, pve$PC)
xl <- sprintf("PC1 (%.1f%%)", ve["PC1"])
yl <- sprintf("PC2 (%.1f%%)", ve["PC2"])
p2a <- ggplot(pd %>% filter(!is.na(Sub)), aes(PC1, PC2, color = Sub)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = CS) +
  labs(x = xl, y = yl, color = "") +
  tn(8) + theme(legend.position = "bottom")
sv(p2a, "Fig2a")

deg$sig <- case_when(
  deg$padj < 0.05 & deg$log2FoldChange > 0.585 ~ "Up (padj<0.05, |LFC|>0.585)",
  deg$padj < 0.05 & deg$log2FoldChange < -0.585 ~ "Down (padj<0.05, |LFC|>0.585)",
  deg$padj < 0.05 ~ "Sig (padj<0.05)",
  TRUE ~ "NS (padj>=0.05)")
deg$lab <- ifelse(!is.na(deg$padj) & deg$padj < 0.05, deg$Gene_Symbol, NA)
p2b <- ggplot(deg, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(data = deg %>% filter(grepl("^NS", sig)), size = 1.0, alpha = 0.05) +
  geom_point(data = deg %>% filter(!grepl("^NS", sig)), size = 1.8, alpha = 0.9) +
  geom_text_repel(data = deg %>% filter(!is.na(lab)), aes(label = lab),
                  size = 2.2, fontface = "italic", color = "black",
                  show.legend = FALSE, max.overlaps = 10,
                  segment.size = 0.3, segment.color = "grey40",
                  min.segment.length = 0, box.padding = 0.5,
                  point.padding = 0.3, force = 2, force_pull = 0.5, seed = 42,
                  max.time = 10, max.iter = 1e6) +
  scale_color_manual(values = c(
    "Up (padj<0.05, |LFC|>0.585)" = "#C0392B",
    "Down (padj<0.05, |LFC|>0.585)" = "#2980B9",
    "Sig (padj<0.05)" = "#E67E22",
    "NS (padj>=0.05)" = "grey70"),
    guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "grey50", linewidth = 0.2) +
  xlim(-2, 2) + labs(x = "Shrunken log2FC", y = "-log10(padj)", color = "") +
  tn(8) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 2.5, alpha = 1))) + theme(legend.position = "bottom", legend.location = "plot", legend.justification = "left", legend.margin = margin(-4, 0, 0, 0), plot.margin = unit(c(2, 6, 2, 3), "mm"),
                legend.text = element_text(size = 9.5),
                legend.key.size = unit(3.5, "mm"))
sv(p2b, "Fig2b")

prep <- function(t) { t <- t %>% mutate(D = str_trunc(Description, 44))
  t <- t[order(t$NES), ]; t$D <- fct_inorder(t$D); t }
mk2c <- function(tgo) ggplot(tgo, aes(NES, D, fill = NES > 0)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9"),
                    labels = c("TRUE" = "Up in RA", "FALSE" = "Down in RA")) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(x = "NES", y = "", fill = "") +
  tn(7) + theme(axis.text.y = element_text(size = 8),
                legend.position = "bottom",
                legend.text = element_text(size = 7),
                plot.margin = unit(c(1, 2, 4, 2), "mm"))
tgo <- prep(sel_terms(ggo))   # GO BP: five up / five down terms by adjusted P, duplicates removed
p2c <- mk2c(tgo); sv(p2c, "Fig2c")

mk2d <- function(tkk) ggplot(tkk, aes(NES, D, fill = NES > 0)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2980B9"),
                    labels = c("TRUE" = "Up in RA", "FALSE" = "Down in RA")) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(x = "NES", y = "", fill = "") + tn(7) + theme(axis.text.y = element_text(size = 8), legend.position = "bottom", legend.text = element_text(size = 7), plot.margin = unit(c(1, 2, 4, 2), "mm"))
tkk <- prep(sel_orig(gkk))   # KEGG: the ten terms with the lowest adjusted P
p2d <- mk2d(tkk); sv(p2d, "Fig2d")

ig <- deg[grepl("^IGHG|^IGHM$|^IGKC$|^IGHA[12]$|^JCHAIN$", deg$Gene_Symbol) &
            !is.na(deg$padj), ] %>% arrange(desc(abs(log2FoldChange))) %>% head(7)
ig$Gene_Symbol <- fct_reorder(ig$Gene_Symbol, ig$log2FoldChange)
ig$sg <- ifelse(ig$padj < 0.05, "padj < 0.05", "padj >= 0.05")
p2h <- ggplot(ig, aes(log2FoldChange, Gene_Symbol, fill = sg)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("padj < 0.05" = "#C0392B", "padj >= 0.05" = "#95A5A6")) +
  geom_vline(xintercept = 0, linewidth = 0.2) +
  labs(x = "log2FC (RA vs Sarcoidosis)", y = "", fill = "") +
  tn(8) + theme(legend.position = "bottom",
                legend.text = element_text(size = 6),
                axis.text.y = element_text(face = "italic"))
sv(p2h, "Fig2g")

chk <- list()
ok <- function(name, cond, detail = "") {
  chk[[length(chk) + 1]] <<- data.frame(check = name, pass = isTRUE(cond), detail = detail)
  if (!isTRUE(cond)) stop("Verification failed: ", name, " ", detail)
}
m <- setNames(met$Value, met$Metric)

b <- ggplot_build(p2a)$data[[1]]
ok("2a n points = 31", nrow(b) == 31)
ok("2a points == PCA_scores", isTRUE(all.equal(sort(b$x), sort(pcs$PC1))) &&
     isTRUE(all.equal(sort(b$y), sort(pcs$PC2))))
ok("2a axis labels", xl == sprintf("PC1 (%.1f%%)", as.numeric(m["PCA_PC1_pct"])) &&
     yl == sprintf("PC2 (%.1f%%)", as.numeric(m["PCA_PC2_pct"])), paste(xl, yl))
ok("2a subgroup counts", TRUE, paste(names(table(pd$Sub)), table(pd$Sub), collapse = ", "))

nsig <- table(deg$sig)
ok("2b padj<0.05 = 143 (metrics)", sum(deg$padj < 0.05, na.rm = TRUE) == as.numeric(m["DEG_padj_lt_0.05"]))
ok("2b Up = metrics", nsig[["Up (padj<0.05, |LFC|>0.585)"]] == as.numeric(m["DEG_absLFC_gt_0.585_up"]),
   sprintf("Up %d", nsig[["Up (padj<0.05, |LFC|>0.585)"]]))
ok("2b Down = metrics", nsig[["Down (padj<0.05, |LFC|>0.585)"]] == as.numeric(m["DEG_absLFC_gt_0.585_down"]),
   sprintf("Down %d", nsig[["Down (padj<0.05, |LFC|>0.585)"]]))
bb <- ggplot_build(p2b)$data
ok("2b no significant gene dropped by xlim", nrow(bb[[2]]) == 143)
lab_df <- bb[[3]]
ok("2b labels are padj<0.05 genes", all(lab_df$label %in% deg$Gene_Symbol[deg$padj < 0.05 & !is.na(deg$padj)]),
   sprintf("%d label candidates (ggrepel hides overlapping ones, max.overlaps=10)", nrow(lab_df)))

bc <- ggplot_build(p2c)$data[[1]]
ok("2c NES == GSEA_GO_BP selected", isTRUE(all.equal(sort(bc$x), sort(tgo$NES))))
bd <- ggplot_build(p2d)$data[[1]]
ok("2d NES == GSEA_KEGG selected", isTRUE(all.equal(sort(bd$x), sort(tkk$NES))))
ok("2d = the first 10 rows by stable order of p.adjust",
   identical(as.character(tkk$ID[order(match(tkk$ID, gkk$ID))]), head(gkk$ID[order(gkk$p.adjust)], 10)) &&
   nrow(tkk) == 10, sprintf("%d terms tied at the 10th p.adjust (%.3g)",
   sum(gkk$p.adjust == sort(gkk$p.adjust)[10]), sort(gkk$p.adjust)[10]))
ok("2c NES drawn", isTRUE(all.equal(sort(ggplot_build(p2c)$data[[1]]$x), sort(tgo$NES))))
for (pp in list(list("2c", ggo, tgo))) {
  g <- pp[[2]]; t <- pp[[3]]; gs <- g[g$p.adjust < 0.05, ]
  for (dir in c(1, -1)) {
    gd <- gs[sign(gs$NES) == dir, ]; td <- t[sign(t$NES) == dir, ]
    exp_ids <- head(dedup(gd[gd$p.adjust < 0.05, ])$ID, 5)
    ok(sprintf("%s %s: top5 by p.adjust then |NES|", pp[[1]], ifelse(dir > 0, "up", "down")),
       setequal(td$ID, exp_ids) && nrow(td) == min(5, nrow(gd)) &&
       all(td$p.adjust <= min(c(Inf, dedup(gd)$p.adjust[!dedup(gd)$ID %in% td$ID]))),
       sprintf("%d shown of %d significant", nrow(td), nrow(gd)))
  }
}

bh <- ggplot_build(p2h)$data[[1]]
ok("2g 7 genes", nrow(ig) == 7)
ok("2g LFC == DEG table", isTRUE(all.equal(sort(bh$x), sort(ig$log2FoldChange))))
mm <- merge(igt, ig[, c("Gene_Symbol", "log2FoldChange", "padj")], by = "Gene_Symbol")
ok("2g IGHG3/IGHG1/IGHG4/IGHM == IGH_genes.csv", nrow(mm) == 4 &&
     isTRUE(all.equal(mm$log2FoldChange.x, mm$log2FoldChange.y)) &&
     isTRUE(all.equal(mm$padj.x, mm$padj.y)))
ok("2g IGHM grey (padj>=0.05)", ig$sg[ig$Gene_Symbol == "IGHM"] == "padj >= 0.05",
   sprintf("IGHM padj=%.4f", ig$padj[ig$Gene_Symbol == "IGHM"]))

for (nm in c("Fig2a", "Fig2b", "Fig2c", "Fig2d", "Fig2g")) {
  d <- dim(png::readPNG(file.path(panel_dir, paste0(nm, ".png"))))
  ok(paste(nm, "pixel size"), d[2] == pmap[[nm]]$px[1] && d[1] == pmap[[nm]]$px[2],
     sprintf("%d x %d", d[2], d[1]))
}

chk <- do.call(rbind, chk); print(chk, right = FALSE)

vals <- bind_rows(
  data.frame(panel = "Fig2a", item = c("x label", "y label"), value = c(xl, yl)),
  data.frame(panel = "Fig2c", item = as.character(tgo$Description),
             value = sprintf("NES=%.3f padj=%.3g", tgo$NES, tgo$p.adjust)),
  data.frame(panel = "Fig2d", item = as.character(tkk$Description),
             value = sprintf("NES=%.3f padj=%.3g", tkk$NES, tkk$p.adjust)),
  data.frame(panel = "Fig2g", item = as.character(ig$Gene_Symbol),
             value = sprintf("LFC=%.3f padj=%.4f", ig$log2FoldChange, ig$padj)))
write.csv(vals, file.path(panel_dir, "panels_transcriptome_values.csv"), row.names = FALSE)
write.csv(chk, file.path(panel_dir, "panels_transcriptome_checks.csv"), row.names = FALSE)
print(vals, right = FALSE)
writeLines(capture.output(sessionInfo()), file.path(panel_dir, "panels_transcriptome_sessionInfo.txt"))

