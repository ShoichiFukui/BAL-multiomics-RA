# 05_RA_only_MOFA2.R -- RA-only MOFA2, five views, expression = VST of the 24 RA samples (Figure 5d-f; Supplementary Fig. 4, whose P values
#   are written to SFig4_medication_n31.csv (source of the figure); Supplementary Data 9 = SD9_weights_n31.csv; infection DeLong test of Figure 5f). RA Factor 1 scores
#   (SourceData_Fig5_n31.csv) feed Figure 5g-i, which are drawn by 07_pulmonary_function_panels.R.
# Produces: mofa_RA/ model, statistics (Results), SFig4 P values, SD9, Source Data table and panels.
# Inputs: results/RA_ILD_Workspace.RData, results/tables/master_data_with_CT_all.csv, planA_data_RA24.csv, medication_at_BAL.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
MOFA_PYTHON <- Sys.getenv("MOFA_PYTHON", "/usr/bin/python3"); Sys.setenv(RETICULATE_PYTHON = MOFA_PYTHON)
suppressPackageStartupMessages({ library(reticulate); use_python(MOFA_PYTHON, required = TRUE)
  library(pROC); library(MOFA2); library(DESeq2); library(ggplot2); library(patchwork) })
pf <- function(p) ifelse(p < 0.001, "p < 0.001", ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p)))  # figure format = manuscript text
rf <- function(r) sprintf("ρ = %.2f", r)
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "mofa_RA"); PD <- file.path(ART_FIGS, "n31", "panels")
dir.create(OD, recursive = TRUE, showWarnings = FALSE); dir.create(PD, recursive = TRUE, showWarnings = FALSE)

load("results/RA_ILD_Workspace.RData")
master_data <- master_data[!master_data$Sample_ID %in% EXCL_ALL, ]
pa <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))
ra_ids <- pa$Sample_ID   # the 24 RA patients, in a fixed sample order
stopifnot(length(ra_ids) == 24, all(sample_group(ra_ids) == "RA"), !any(ra_ids %in% EXCL_ALL))
inf <- pa$respiratory_infection[match(ra_ids, pa$Sample_ID)]; stopifnot(sum(inf == 1) == 6, all(inf %in% 0:1))
dOpa <- -pa$d_aer[match(ra_ids, pa$Sample_ID)]
cm <- count_matrix
vst_on <- function(ids, design) {
  cd <- data.frame(G = factor(ifelse(is_sarc(ids), "Control", "RA")), row.names = ids)
  d <- DESeqDataSetFromMatrix(round(cm[, ids]), cd, if (design == "group") ~G else ~1)
  d <- d[rowSums(counts(d) >= 10) >= 5, ]; d <- DESeq(d, quiet = TRUE); assay(vst(d, blind = FALSE)) }
EXPR <- list(C_vst_RA24 = vst_on(ra_ids, "none"))
mk <- function(cols) { m <- as.matrix(master_data[match(ra_ids, master_data$Sample_ID), cols]); rownames(m) <- ra_ids; m }
# BAL cytokines: names containing a cell or FCM term are dropped; "CD" also drops BALF_sCD40L, so the BAL view has 40 cytokines (Supplementary Methods)
bc <- grep("^BALF_", colnames(master_data), value = TRUE); bc <- bc[!grepl("Th17|Treg|Th1|Th2|Macro|Neutro|CD|Lymph|Eosino|Baso|Plasma|Other|Percent", bc)]
bc <- bc[sapply(bc, function(x) is.numeric(master_data[[x]]))]; balf_view <- mk(bc); balf_view <- balf_view[, colSums(!is.na(balf_view)) >= 18]
sc <- grep("^Serum_", colnames(master_data), value = TRUE); serum_view <- mk(sc); serum_view <- serum_view[, colSums(!is.na(serum_view)) >= 18]
balf_facs <- mk(intersect(c("BALF_Th17.1_CCR6pos_CXCR3pos", "BALF_Th17_CCR6pos_CXCR3neg", "BALF_Treg_CD127neg_CD25pos", "BALF_Th1_CCR6neg_CXCR3pos",
  "BALF_Th2_CCR6neg_CXCR3neg", "BALF_Macrophage_Percent", "BALF_Lymphocyte_Percent", "BALF_Neutrophil_Percent"), colnames(master_data)))
pb_facs <- mk(intersect(c("PB_Th17.1_CCR6pos_CXCR3pos", "PB_Th17_CCR6pos_CXCR3neg", "PB_Treg_CD127neg_CD25pos", "PB_Th1_CCR6neg_CXCR3pos",
  "PB_Th2_CCR6neg_CXCR3neg"), colnames(master_data)))
sp <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = TRUE)); c(unname(r$estimate), r$p.value) }

ex <- EXPR$C_vst_RA24[, ra_ids]; top <- names(sort(apply(ex, 1, var), decreasing = TRUE))[1:30]
dl <- list(Expression = ex[top, ], BALF_Cytokine = t(balf_view), Serum_Cytokine = t(serum_view), BALF_FACS = t(balf_facs), PB_FACS = t(pb_facs))
for (v in names(dl)) { m <- t(dl[[v]]); for (j in 1:ncol(m)) { s <- sd(m[, j], na.rm = TRUE); if (!is.na(s) && s > 0) m[, j] <- scale(m[, j]) }; dl[[v]] <- t(m) }
mo0 <- create_mofa(dl); mo <- get_default_model_options(mo0); mo$num_factors <- 3
to <- get_default_training_options(mo0); to$maxiter <- 1000; to$seed <- 42; to$verbose <- FALSE; to$convergence_mode <- "medium"
mo0 <- prepare_mofa(mo0, model_options = mo, training_options = to)
tr <- run_mofa(mo0, outfile = file.path(OD, "MOFA2_RA_only_RA24vst.hdf5"), use_basilisk = FALSE)
ra_factors <- get_factors(tr)[[1]][ra_ids, ]; colnames(ra_factors) <- paste0("Factor", 1:3)   # MOFA order and signs (seed 42)
r2 <- get_variance_explained(tr)$r2_per_factor[[1]]; rownames(r2) <- paste0("Factor", 1:3)
ra_weights <- lapply(get_weights(tr), function(w) { colnames(w) <- paste0("Factor", 1:3); w })

ra_r2 <- list(r2_per_factor = list(group1 = r2))
save(ra_factors, ra_r2, ra_weights, file = file.path(OD, "MOFA2_RA_only_Results_RA24vst.RData"))

th171 <- as.numeric(master_data$BALF_Th17.1_CCR6pos_CXCR3pos[match(ra_ids, master_data$Sample_ID)])
pbth17 <- as.numeric(master_data$PB_Th17_CCR6pos_CXCR3neg[match(ra_ids, master_data$Sample_ID)])
pbth171 <- as.numeric(master_data$PB_Th17.1_CCR6pos_CXCR3pos[match(ra_ids, master_data$Sample_ID)])
mc <- read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE); mc <- mc[!duplicated(mc$Sample_ID), ]
bneut <- mc$BALF_Neutrophil_Percent[match(ra_ids, mc$Sample_ID)]
roc_th <- roc(inf, th171, quiet = TRUE, direction = "auto")
ST <- list(); f1 <- as.numeric(ra_factors[, 1])
rf1 <- roc(inf, f1, quiet = TRUE, direction = "auto"); dl_p <- roc.test(roc_th, rf1, method = "delong")$p.value
ST$n31 <- c(Inf_p = wilcox.test(f1[inf == 1], f1[inf == 0], exact = TRUE)$p.value, Inf_AUC = as.numeric(auc(rf1)), Th171_AUC = as.numeric(auc(roc_th)),
            DeLong_p = dl_p, rho_dOpa = sp(f1, dOpa)[1], p_dOpa = sp(f1, dOpa)[2])
stats <- data.frame(Statistic = names(ST$n31), n31_model = ST$n31)
write.csv(stats, file.path(OD, "RA_MOFA2_n31_statistics.csv"), row.names = FALSE)
medv <- list(GC = "GC", MTX = "MTX", NSAIDs = "NSAIDS")   # GC and MTX at BAL from medication_at_BAL.csv (med_at_BAL); NSAIDs from the clinical sheet
med <- do.call(rbind, lapply(list(list("RA Factor 1", as.numeric(ra_factors[, 1])), list("BAL Th17.1", th171)), function(vv) do.call(rbind, lapply(names(medv), function(k) {
  mv <- if (k == "NSAIDs") mc[[medv[[k]]]][match(ra_ids, mc$Sample_ID)] else med_at_BAL(ra_ids, c(GC = "GC", MTX = "MTX")[[k]]); ok <- !is.na(vv[[2]]) & !is.na(mv); on <- vv[[2]][ok & mv == 1]; off <- vv[[2]][ok & mv == 0]
  data.frame(Medication = k, Variable = vv[[1]], N_on = length(on), N_off = length(off), Median_on = round(median(on), 2), Median_off = round(median(off), 2),
             P = suppressWarnings(wilcox.test(on, off, exact = TRUE)$p.value)) }))))
med$BH <- p.adjust(med$P, "BH"); print(med, row.names = FALSE, digits = 3)   # BH across the six rows
pm <- function(k, v) med$P[med$Medication == k & med$Variable == v]
stopifnot(med$N_on[med$Medication == "GC"][1] == 7, abs(pm("GC", "RA Factor 1") - 0.576) < 5e-4, abs(pm("GC", "BAL Th17.1") - 0.775) < 5e-4,
          abs(pm("MTX", "RA Factor 1") - 0.452) < 5e-4, abs(pm("MTX", "BAL Th17.1") - 0.245) < 5e-4)   # Supplementary Fig. 4 P values (medication at BAL)
write.csv(med, file.path(OD, "SFig4_medication_n31.csv"), row.names = FALSE)
vmap <- c(Expression = "Expression", BALF_Cytokine = "BALF_Cytokine", Serum_Cytokine = "Serum_Cytokine", BALF_FACS = "BALF_FCM", PB_FACS = "PB_FCM")
sd9 <- function(w) do.call(rbind, lapply(names(vmap), function(v) do.call(rbind, lapply(1:3, function(k)
  data.frame(View = vmap[[v]], Factor = k, Feature = rownames(w[[v]]), Weight = round(w[[v]][, k], 5))))))
write.csv(sd9(ra_weights), file.path(OD, "SD9_weights_n31.csv"), row.names = FALSE)
src_df <- data.frame(Sample_ID = ra_ids, RA_Factor1 = as.numeric(ra_factors[, 1]), Infection = inf, dOpacity = dOpa,
                     PB_Th17 = pbth17, PB_Th17.1 = pbth171, BAL_neutrophil_pct = bneut, BAL_Th17.1 = th171)
write.csv(src_df, file.path(OD, "SourceData_Fig5_n31.csv"), row.names = FALSE)

svp <- function(p, nm, px) { f <- file.path(PD, nm); ggsave(paste0(f, ".png"), p, width = px[1], height = px[2], units = "px", dpi = 300, bg = "white")
  ggsave(paste0(f, ".pdf"), p, width = px[1] / 300 * 25.4, height = px[2] / 300 * 25.4, units = "mm", device = cairo_pdf)
  d <- dim(png::readPNG(paste0(f, ".png"))); stopifnot(d[2] == px[1], d[1] == px[2]) }
f1 <- as.numeric(ra_factors[, 1])
r2d <- r2; rownames(r2d) <- paste("Factor", seq_len(nrow(r2d)))   # row labels "Factor 1" ...
r2l <- do.call(rbind, lapply(rownames(r2d), function(fk) data.frame(Factor = fk, View = colnames(r2d), R2 = r2d[fk, ])))
r2l$View <- c(BALF_FACS = "BAL FCM", PB_FACS = "PB FCM", Expression = "Expression", BALF_Cytokine = "BAL Cytokine", Serum_Cytokine = "Serum Cytokine")[r2l$View]
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2d))); r2l$View <- factor(r2l$View, levels = c("Expression", "BAL Cytokine", "Serum Cytokine", "BAL FCM", "PB FCM"))
tn3 <- function(bs = 8) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"), axis.title = element_text(size = rel(0.95), face = "bold"),
        plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 3, 0)), panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 3), "mm"))
p5d <- ggplot(r2l, aes(View, Factor, fill = R2)) + geom_tile(color = "white", linewidth = 0.5) + scale_fill_gradient(low = "white", high = "#2171B5", name = "R² (%)") +
  geom_text(aes(label = ifelse(R2 >= 1, sprintf("%.1f", R2), ""), color = ifelse(R2 >= 25, "white", "black")), size = 2.5, fontface = "bold", show.legend = FALSE) +
  geom_text(aes(label = ifelse(R2 > 0 & R2 < 1, sprintf("%.1f", R2), "")), size = 2.5, color = "grey30") + scale_color_identity() + labs(x = "", y = "") +
  tn3(8) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 7), legend.position = "right", legend.key.height = unit(12, "mm"), legend.key.width = unit(3, "mm"), legend.title = element_text(size = 7), legend.text = element_text(size = 6))
svp(p5d, "Fig5d", c(1181, 767))
th13 <- function(bs = 13) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(color = "black"), axis.title = element_text(face = "bold"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4), axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 2), "mm"))
d5e <- data.frame(F1 = f1, G = factor(ifelse(inf == 1, "Infection (+)", "Infection (−)"), levels = c("Infection (−)", "Infection (+)")))
p5e <- ggplot(d5e, aes(G, F1, fill = G)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.85) +
  geom_jitter(width = 0.08, size = 1.5, alpha = 0.7, color = "grey15") + scale_fill_manual(values = c("Infection (−)" = "#5DADE2", "Infection (+)" = "#EC7063"), guide = "none") +
  annotate("text", x = 1.5, y = max(f1) * 1.12, label = pf(ST$n31["Inf_p"]), size = 3.3) + labs(x = NULL, y = "RA Factor 1 value") + th13(10.5)
set.seed(42); svp(p5e, "Fig5e", c(825, 741))
CUD_ORANGE <- "#FF9900"; CUD_BLUE <- "#0041FF"
add_ep <- function(spec, sens) { o <- order(spec, sens); spec <- spec[o]; sens <- sens[o]
  if (spec[1] != 0 || sens[1] != 0) { spec <- c(0, spec); sens <- c(0, sens) }; if (tail(spec, 1) != 1 || tail(sens, 1) != 1) { spec <- c(spec, 1); sens <- c(sens, 1) }; list(spec = spec, sens = sens) }
rf1 <- roc(inf, f1, quiet = TRUE, direction = "auto"); e1 <- add_ep(1 - roc_th$specificities, roc_th$sensitivities); e2 <- add_ep(1 - rf1$specificities, rf1$sensitivities)
rd <- rbind(data.frame(Sens = e1$sens, Spec = e1$spec, Model = sprintf("BAL Th17.1 %% (AUC = %.3f)", auc(roc_th))), data.frame(Sens = e2$sens, Spec = e2$spec, Model = sprintf("RA Factor 1 (AUC = %.3f)", auc(rf1))))
rd$Model <- factor(rd$Model, levels = unique(rd$Model))
tn <- function(bs = 8) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"), axis.title = element_text(size = rel(0.95), face = "bold"),
        plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 3, 0)), panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
        axis.line = element_blank(), plot.margin = unit(c(2, 3, 2, 3), "mm"))
p5f <- ggplot(rd, aes(Spec, Sens, color = Model)) + geom_step(linewidth = 0.7, direction = "vh") + geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.3) +
  scale_color_manual(values = c(CUD_ORANGE, CUD_BLUE)) + labs(x = "1 – Specificity", y = "Sensitivity", color = "") + ggtitle("Respiratory infection (RA)") +
  tn(8) + theme(legend.position = c(0.62, 0.22), legend.background = element_rect(fill = alpha("white", 0.95), color = NA), legend.text = element_text(size = 6.5), legend.key.height = unit(3, "mm"))
svp(p5f, "Fig5f", c(944, 885))
RED <- "#C0392B"; BLU <- "#2980B9"
tn7 <- function(bs = 7) theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
  theme(text = element_text(color = "black"), axis.text = element_text(size = rel(0.9), color = "black"), axis.title = element_text(size = rel(0.95), face = "bold"),
        plot.title = element_text(size = rel(1.0), face = "bold", hjust = 0.5, margin = ggplot2::margin(0, 0, 2, 0)), panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
        axis.line = element_blank(), plot.margin = unit(c(2, 2, 2, 2), "mm"))
fmt_p <- function(p) pf(p)
plots <- list()
for (vv in list(list("RA Factor 1", f1), list("BAL Th17.1 (%)", th171))) for (k in names(medv)) {
  mv <- if (k == "NSAIDs") mc[[medv[[k]]]][match(ra_ids, mc$Sample_ID)] else med_at_BAL(ra_ids, c(GC = "GC", MTX = "MTX")[[k]]); ok <- !is.na(vv[[2]]) & !is.na(mv); on <- vv[[2]][ok & mv == 1]; off <- vv[[2]][ok & mv == 0]
  w <- suppressWarnings(wilcox.test(on, off, exact = TRUE))
  df <- data.frame(Value = c(on, off), Med = factor(c(rep(paste0(k, " (+)"), length(on)), rep(paste0(k, " (−)"), length(off))), levels = c(paste0(k, " (−)"), paste0(k, " (+)"))))
  yr <- range(df$Value); ymax <- yr[2] + diff(yr) * 0.25
  plots[[length(plots) + 1]] <- ggplot(df, aes(Med, Value, fill = Med)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, alpha = 0.7) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.6, show.legend = FALSE) + scale_fill_manual(values = c(BLU, RED)) +
    labs(x = "", y = vv[[1]], title = k) + annotate("text", x = 1.5, y = ymax, label = fmt_p(w$p.value), size = 2.2) +
    coord_cartesian(ylim = c(min(yr[1], 0), ymax * 1.05), clip = "off") + tn7(7) + theme(legend.position = "none", axis.text.x = element_text(size = 6)) }
ps4 <- wrap_plots(plots, ncol = 3) + plot_annotation(theme = theme(plot.margin = ggplot2::margin(2, 2, 2, 2)))
svp(ps4, "SFig4", c(1889, 1181))   # the jitter continues the random stream started by set.seed(42) before Figure 5e, so the image is reproducible
print(round(r2, 1))
