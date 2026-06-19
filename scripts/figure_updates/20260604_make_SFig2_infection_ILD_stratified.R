# ============================================================================
# 20260604_make_SFig2_infection_ILD_stratified.R
# Supplementary Fig. 2 "infection x ILD-stratified" panel: BAL Th17.1 and
# GSVA OXPHOS by prospective respiratory-infection status, stratified by
# baseline ILD status. The infection biomarkers are ILD-independent
# (ILD-: Th17.1 p=0.014, OXPHOS p=0.007; ILD+: NS).
#
# Output: results/panels/FigS2_infection_ILD_stratified.{png,pdf}
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(ggplot2) })
EXCL <- character(0)
load("results/RA_ILD_Workspace.RData")
gsva <- get("gsva_scores")
panel_dir <- "results/panels"

md <- master_data[!master_data$Sample_ID %in% EXCL, ]
ra <- md[!grepl("^(KYC|Sarcoidosis)", md$Sample_ID), ]                       # RA only, n=24
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

ggsave(file.path(panel_dir, "FigS2_infection_ILD_stratified.pdf"), p,
       width = 95, height = 85, units = "mm", dpi = 300, device = cairo_pdf)
ggsave(file.path(panel_dir, "FigS2_infection_ILD_stratified.png"), p,
       width = 95, height = 85, units = "mm", dpi = 300, bg = "white")
cat("Saved FigS2_infection_ILD_stratified (95x85mm)\n")
