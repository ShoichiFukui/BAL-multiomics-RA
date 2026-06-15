# ============================================================================
# 20260604_regen_FigS3_alpha_consistent.R
# Regenerate SuppFig3 so the Greek-alpha panel titles (GRO-α, IFN-α2) render in
# the SAME bold style as the other titles. Helvetica-bold lacks α (font
# fallback made GRO-α/IFN-α2 look different); render those two titles as bold
# plotmath expressions so the α is a proper bold symbol.
# Healthy comparison remains descriptive (green 95% CI band); the in-panel
# p-value is the RA-vs-Sarcoidosis Wilcoxon test (the disease comparison).
# Overwrites results/panels/FigS3_healthy_serum_95CI.{png,pdf}
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
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
ggsave(file.path(panel_dir,"FigS3_healthy_serum_95CI.pdf"), p_s3, width=250, height=110, units="mm", dpi=300, device=cairo_pdf)
ggsave(file.path(panel_dir,"FigS3_healthy_serum_95CI.png"), p_s3, width=250, height=110, units="mm", dpi=300, bg="white")
cat("Saved FigS3_healthy_serum_95CI (alpha titles as bold plotmath)\n")
