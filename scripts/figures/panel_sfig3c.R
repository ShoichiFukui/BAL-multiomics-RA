# panel_sfig3c.R -- Supplementary Fig. 3C: BAL Th17.1 proportion by subsequent respiratory infection, stratified by baseline ILD
#   (ILD- / ILD+; exact Wilcoxon, p = 0.014 / 0.19).
# Produces: <ART_FIGS>/n31/panels/SFig3c.png and .pdf (95 x 85 mm, 1122 x 1003 px at 300 dpi).
# Inputs: results/RA_ILD_Workspace.RData (master_data: BAL Th17.1 %, ILD score), planA_data_RA24.csv (infection outcome; controlled access).
# The jitter is drawn with set.seed(42) immediately before saving. The cohort is read through 00_cohort_config.R.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(tidyverse) })
PD <- file.path(ART_FIGS, "n31", "panels"); dir.create(PD, recursive = TRUE, showWarnings = FALSE)
sv <- function(p, name, w, h, px) {
  ggsave(file.path(PD, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", dpi = 300, device = cairo_pdf)
  ggsave(file.path(PD, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  d <- dim(png::readPNG(file.path(PD, paste0(name, ".png")))); stopifnot(d[2] == px[1], d[1] == px[2]) }
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e); md <- e$master_data; rm(e)
pa <- read.csv(file.path(ART_TABLES, "n31", "planA/planA_data_RA24.csv"))
ra <- md[md$Sample_ID %in% pa$Sample_ID & sample_group(md$Sample_ID) %in% "RA", ]; stopifnot(nrow(ra) == 24)
inf <- pa$respiratory_infection[match(ra$Sample_ID, pa$Sample_ID)]; stopifnot(sum(inf) == 6)
ra$ILD <- factor(ifelse(ra$ILD_Score > 0, "ILD+", "ILD−"), levels = c("ILD−", "ILD+")); ra$Infection <- factor(ifelse(inf == 1, "+", "−"), levels = c("−", "+"))
ra$`BAL Th17.1 (%)` <- as.numeric(ra$BALF_Th17.1_CCR6pos_CXCR3pos)
long <- ra %>% select(Sample_ID, ILD, Infection, `BAL Th17.1 (%)`) %>% pivot_longer(c(`BAL Th17.1 (%)`), names_to = "Variable", values_to = "Value") %>% filter(!is.na(Value) & !is.na(Infection))
long$Variable <- factor(long$Variable, levels = "BAL Th17.1 (%)")
pdat <- long %>% group_by(Variable, ILD) %>% summarise(p = tryCatch(suppressWarnings(wilcox.test(Value ~ Infection, exact = TRUE)$p.value), error = function(e) NA), ymax = max(Value), .groups = "drop")
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))  # figure format = manuscript text
pdat$lab <- pf(pdat$p); print(pdat)
stopifnot(nrow(pdat) == 2, abs(pdat$p[pdat$ILD == "ILD−"] - 0.0141) < 5e-4, abs(pdat$p[pdat$ILD == "ILD+"] - 0.194) < 5e-4)   # p = 0.014 / 0.19
tn2 <- function(bs = 8) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.85), color = "black"), axis.title = element_text(size = rel(0.95), face = "bold"),
        strip.text = element_text(size = rel(0.9), face = "bold"), strip.background = element_rect(fill = "grey92", color = NA),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33), axis.line = element_blank(), legend.position = "none", plot.margin = unit(c(2, 3, 2, 3), "mm"))
p <- ggplot(long, aes(Infection, Value)) + geom_boxplot(aes(fill = Infection), outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.7, color = "grey20") + facet_grid(Variable ~ ILD, scales = "free_y", switch = "y") +
  geom_text(data = pdat, aes(x = 1.5, y = ymax, label = lab), inherit.aes = FALSE, size = 2.5, vjust = -0.2) +
  scale_fill_manual(values = c("−" = "#4C78A8", "+" = "#C0392B")) + scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = "Respiratory infection", y = NULL) + tn2(8) + theme(strip.placement = "outside", strip.background.y = element_blank())   # the row strip is the y-axis title (no grey band)
set.seed(42); sv(p, "SFig3c", 95, 85, c(1122, 1003))
