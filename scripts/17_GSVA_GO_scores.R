# 17_GSVA_GO_scores.R -- per-sample gene-set scores for oxidative phosphorylation and surfactant homeostasis: GSVA (kcdf Gaussian) of the GO
#   Biological Process sets GOBP_OXIDATIVE_PHOSPHORYLATION (GO:0006119, the same GO term as the pre-ranked GSEA) and
#   GOBP_SURFACTANT_HOMEOSTASIS (msigdbr 25.1.1, collection C5, subcollection GO:BP) on the log-CPM of the 31 analysed samples (genes with
#   >= 10 reads in >= 5 samples). The scores are RNA-seq-derived variables listed outside the FDR family (like the BayesPrism estimates).
#   Used for Figure 3H, the composition-adjusted model (08) and the derived rows of Supplementary Data 3 (06).
#   Statistics: OXPHOS score vs infection (exact Wilcoxon, Cliff's delta, AUC; manuscript P = 0.0094), adjusted for the BAL macrophage /
#   lymphocyte / neutrophil % (lm; P = 0.0279), ILD-stratified (ILD- P = 0.112, ILD+ P = 0.133; descriptive, not shown since
#   Supplementary Fig. 3C was reduced to BAL Th17.1), surfactant score vs infection (P = 0.537).
# Produces: <ART_TABLES>/n31/gsva_GO/gsva_GO_scores_n31.csv, gsva_GO_gene_sets.csv, oxphos_gsva_stats.csv, sessionInfo_gsva.txt.
# Inputs: results/RA_ILD_Workspace.RData (count_matrix), planA_data_RA24.csv, results/tables/master_data_with_CT_all.csv (controlled access).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(msigdbr); library(GSVA); library(pROC) })
R0 <- file.path(ART_TABLES, "n31"); OD <- file.path(R0, "gsva_GO"); dir.create(OD, recursive = TRUE, showWarnings = FALSE)

# ---- GSVA ----
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
ids <- colnames(e$count_matrix); ids <- ids[!ids %in% EXCL_ALL & !is.na(sample_group(ids))]; stopifnot(length(ids) == 31, sum(is_sarc(ids)) == 7)
cnt <- as.matrix(e$count_matrix[, ids]); cnt <- cnt[rowSums(cnt >= 10) >= 5, ]
lcpm <- log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)
GO <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gs <- list(OXPHOS = unique(GO$gene_symbol[GO$gs_name == "GOBP_OXIDATIVE_PHOSPHORYLATION"]),
           Surfactant = unique(GO$gene_symbol[GO$gs_name == "GOBP_SURFACTANT_HOMEOSTASIS"]))

sc <- gsva(gsvaParam(lcpm, gs, kcdf = "Gaussian"), verbose = FALSE)
out <- data.frame(Sample_ID = ids, Group = sample_group(ids), OXPHOS = sc["OXPHOS", ids], Surfactant = sc["Surfactant", ids])
write.csv(out, file.path(OD, "gsva_GO_scores_n31.csv"), row.names = FALSE)
write.csv(data.frame(GeneSet = rep(c("GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_SURFACTANT_HOMEOSTASIS"), lengths(gs)), Gene = unlist(gs, use.names = FALSE),
                     In_data = unlist(gs, use.names = FALSE) %in% rownames(lcpm)), file.path(OD, "gsva_GO_gene_sets.csv"), row.names = FALSE)

# ---- statistics ----
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); ms <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE)); ms <- ms[match(ra$Sample_ID, ms$Sample_ID), ]
num <- function(x) suppressWarnings(as.numeric(x)); inf <- ra$respiratory_infection; stopifnot(sum(inf) == 6)
ox <- out$OXPHOS[match(ra$Sample_ID, out$Sample_ID)]; sf <- out$Surfactant[match(ra$Sample_ID, out$Sample_ID)]
cliff <- function(a, b) mean(outer(a, b, function(x, y) sign(x - y)))
wx <- function(x, g) suppressWarnings(wilcox.test(x[g == 1], x[g == 0], exact = TRUE))$p.value
mac <- num(ms$BALF_Macrophage_Percent); lym <- num(ms$BALF_Lymphocyte_Percent); neu <- num(ms$BALF_Neutrophil_Percent); tot <- num(ms[["BAL_total_cells_per_mL"]])
fit <- summary(lm(ox ~ inf + mac + lym + neu))$coefficients
spr <- function(x, y) { r <- suppressWarnings(cor.test(x, y, method = "spearman", exact = TRUE)); c(rho = unname(r$estimate), p = r$p.value) }
ild <- num(ra$ILD) == 1; th <- num(ra[["BALF_Th17.1_CCR6pos_CXCR3pos"]])
st <- rbind(
  data.frame(Item = "OXPHOS GSVA, infection", Median_pos = median(ox[inf == 1]), Median_neg = median(ox[inf == 0]), P = wx(ox, inf), Cliff = cliff(ox[inf == 1], ox[inf == 0]), AUC = as.numeric(auc(roc(inf, ox, quiet = TRUE)))),
  data.frame(Item = "OXPHOS GSVA, infection adjusted for macrophage/lymphocyte/neutrophil %", Median_pos = NA, Median_neg = NA, P = fit["inf", 4], Cliff = NA, AUC = NA),
  data.frame(Item = "OXPHOS GSVA, infection, ILD-negative", Median_pos = NA, Median_neg = NA, P = wx(ox[!ild], inf[!ild]), Cliff = NA, AUC = NA),
  data.frame(Item = "OXPHOS GSVA, infection, ILD-positive", Median_pos = NA, Median_neg = NA, P = wx(ox[ild], inf[ild]), Cliff = NA, AUC = NA),
  data.frame(Item = "BAL Th17.1, infection, ILD-negative", Median_pos = NA, Median_neg = NA, P = wx(th[!ild], inf[!ild]), Cliff = NA, AUC = NA),
  data.frame(Item = "BAL Th17.1, infection, ILD-positive", Median_pos = NA, Median_neg = NA, P = wx(th[ild], inf[ild]), Cliff = NA, AUC = NA),
  data.frame(Item = "Surfactant GSVA, infection", Median_pos = median(sf[inf == 1]), Median_neg = median(sf[inf == 0]), P = wx(sf, inf), Cliff = cliff(sf[inf == 1], sf[inf == 0]), AUC = NA),
  data.frame(Item = "OXPHOS GSVA, RA vs sarcoidosis", Median_pos = median(out$OXPHOS[out$Group == "RA"]), Median_neg = median(out$OXPHOS[out$Group == "Sarcoidosis"]),
             P = suppressWarnings(wilcox.test(out$OXPHOS[out$Group == "RA"], out$OXPHOS[out$Group == "Sarcoidosis"], exact = TRUE))$p.value, Cliff = NA, AUC = NA))
for (v in list(c("macrophage %", "mac"), c("lymphocyte %", "lym"), c("neutrophil %", "neu"), c("total cells", "tot"), c("BAL Th17.1 %", "th"))) {
  r <- spr(ox, get(v[2])); st <- rbind(st, data.frame(Item = paste("OXPHOS GSVA vs", v[1], "(Spearman rho in Cliff column)"), Median_pos = NA, Median_neg = NA, P = r[["p"]], Cliff = r[["rho"]], AUC = NA)) }

stopifnot(abs(st$P[1] - 0.00942) < 5e-5, abs(st$P[2] - 0.0279) < 5e-4)   # manuscript: P = 0.0094; adjusted P = 0.028
num_cols <- sapply(st, is.numeric); st[num_cols] <- lapply(st[num_cols], function(x) signif(x, 4)); print(st, row.names = FALSE)
write.csv(st, file.path(OD, "oxphos_gsva_stats.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OD, "sessionInfo_gsva.txt"))
