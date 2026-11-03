# 02_transcriptome_RA_vs_sarcoidosis.R -- PCA, DESeq2 and GSEA, RA versus sarcoidosis (Fig. 2A-D, 2G; Suppl. Data 1; results/DEG_RA_vs_sarcoidosis_n31_genome_wide.csv).
# Produces: transcriptome_*.csv in <ART_TABLES>/<cohort>/ (DEG_ShrunkLFC, GSEA_GO_BP, GSEA_KEGG, PCA_scores,
#   PCA_variance, IGH_genes, metrics, sessionInfo); the files released in results/ are listed in results/README.md.
# Inputs: results/RA_ILD_Workspace.RData.
suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(dplyr)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(effsize)
})

PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
cfg <- set_cohort()
out <- function(f) file.path(cfg$tables, paste0("transcriptome_", f))

metrics <- data.frame(Metric = character(), Value = character(), stringsAsFactors = FALSE)
addm <- function(m, v) metrics[nrow(metrics) + 1, ] <<- list(m, as.character(v))

e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
cm_all <- e$count_matrix
gene_annotations <- e$gene_annotations
ids <- setdiff(colnames(cm_all), cfg$EXCL)
ids <- e$common_samples[e$common_samples %in% ids]   # keep original sample order
check_n(ids, cfg, "transcriptome")
grp_lab <- ifelse(is_sarc(ids), "Control", "RA")

addm("n_RA", sum(grp_lab == "RA")); addm("n_Sarc", sum(grp_lab == "Control"))

count_matrix <- as.matrix(cm_all[, ids])
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]

sample_metadata <- data.frame(Sample_Group = factor(grp_lab), row.names = ids)
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_metadata,
                              design = ~Sample_Group)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- DESeq(dds, quiet = TRUE)
vsd <- vst(dds, blind = FALSE); expr_matrix <- assay(vsd)
log_cpm <- cpm(count_matrix[rownames(expr_matrix), ], log = TRUE, prior.count = 1)

addm("genes_after_filter", nrow(expr_matrix))

pca <- prcomp(t(expr_matrix), scale. = TRUE)
ve <- summary(pca)$importance[2, 1:5] * 100

addm("PCA_PC1_pct", sprintf("%.2f", ve[1])); addm("PCA_PC2_pct", sprintf("%.2f", ve[2]))
pc_df <- data.frame(Sample_ID = ids, Group = ifelse(grp_lab == "RA", "RA", "Sarcoidosis"),
                    PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3])
for (pc in c("PC1", "PC2")) {
  w <- wilcox.test(pc_df[[pc]][grp_lab == "RA"], pc_df[[pc]][grp_lab == "Control"], exact = TRUE)
  
  addm(paste0("PCA_", pc, "_RA_vs_Sarc_p"), sprintf("%.4f", w$p.value))
}
write.csv(pc_df, out("PCA_scores.csv"), row.names = FALSE)
write.csv(data.frame(PC = names(ve), Pct_variance = ve), out("PCA_variance.csv"), row.names = FALSE)

deg_meta <- data.frame(Sample_ID = ids, row.names = ids)
deg_meta$Group <- factor(ifelse(is_sarc(ids), "Control", "RA"), levels = c("Control", "RA"))
count_mat_deg <- count_matrix[, ids]
count_mat_deg <- count_mat_deg[rowSums(count_mat_deg >= 10) >= 5, ]
dds2 <- DESeqDataSetFromMatrix(countData = count_mat_deg, colData = deg_meta, design = ~Group)
dds2 <- DESeq(dds2, quiet = TRUE)
res_shrunk <- lfcShrink(dds2, coef = "Group_RA_vs_Control", type = "normal", quiet = TRUE)
rc2 <- assays(dds2)[["replaceCounts"]]
if (!is.null(rc2)) { dds2_r <- dds2; counts(dds2_r) <- rc2
  fx2 <- lfcShrink(dds2_r, coef = "Group_RA_vs_Control", type = "normal", quiet = TRUE)
  res_shrunk$log2FoldChange <- fx2[rownames(res_shrunk), "log2FoldChange"]; res_shrunk$lfcSE <- fx2[rownames(res_shrunk), "lfcSE"] }   # LFC sign agrees with the Wald statistic

res2_df <- as.data.frame(res_shrunk) %>%
  mutate(Gene_Symbol = rownames(res_shrunk)) %>%
  left_join(gene_annotations[, c("Gene_Symbol", "ENTREZID", "GENENAME")], by = "Gene_Symbol") %>%
  arrange(padj)
res2_df$Regulation <- "NS"   # adjusted P < 0.05 with |log2FC| <= 0.585 is labelled sub-threshold, as in Supplementary Data 1
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange > 0] <- "Up (sub-threshold)"
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange < 0] <- "Down (sub-threshold)"
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange > 0.585] <- "Up"
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange < -0.585] <- "Down"
write.csv(res2_df, out(paste0("DEG_ShrunkLFC", ".csv")), row.names = FALSE)

n_sig <- sum(res2_df$padj < 0.05, na.rm = TRUE)
n_up <- sum(res2_df$Regulation == "Up"); n_dn <- sum(res2_df$Regulation == "Down")

addm("DEG_genes_tested", nrow(res2_df)); addm("DEG_padj_lt_0.05", n_sig)
addm("DEG_absLFC_gt_0.585_up", n_up); addm("DEG_absLFC_gt_0.585_down", n_dn)
sig <- res2_df[which(res2_df$padj < 0.05), c("Gene_Symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "Regulation")]
print(sig, row.names = FALSE)
addm("DEG_padj_lt_0.05_genes", paste(sig$Gene_Symbol, collapse = ";"))

igs <- c("IGHG3", "IGHG1", "IGHG4", "IGHM")
ig <- res2_df[match(igs, res2_df$Gene_Symbol), c("Gene_Symbol", "baseMean", "log2FoldChange", "pvalue", "padj")]
ig$rank_by_padj <- match(igs, res2_df$Gene_Symbol)
print(ig, row.names = FALSE)
write.csv(ig, out(paste0("IGH_genes", ".csv")), row.names = FALSE)
for (i in seq_len(nrow(ig))) addm(paste0("IG_", ig$Gene_Symbol[i]),
  sprintf("LFC=%.3f padj=%.4g rank=%d", ig$log2FoldChange[i], ig$padj[i], ig$rank_by_padj[i]))

dds3 <- DESeqDataSetFromMatrix(round(count_matrix[, ids]), deg_meta, ~Group)
dds3 <- dds3[rowSums(counts(dds3)) > 10, ]
dds3 <- DESeq(dds3, quiet = TRUE)
r3 <- as.data.frame(results(dds3, contrast = c("Group", "RA", "Control")))
n3 <- sum(r3$padj < 0.05, na.rm = TRUE)

addm("diag_Sarcoidosis7filter_genes_tested", nrow(r3))
addm("diag_Sarcoidosis7filter_padj_lt_0.05", n3)
addm("diag_Sarcoidosis7filter_sig_genes", paste(rownames(r3)[which(r3$padj < 0.05)], collapse = ";"))
r3u <- as.data.frame(results(dds2))
addm("diag_mainfilter_unshrunk_padj_lt_0.05", sum(r3u$padj < 0.05, na.rm = TRUE))

ranked <- res2_df$log2FoldChange
names(ranked) <- res2_df$ENTREZID
ranked <- ranked[!is.na(names(ranked)) & names(ranked) != ""]
ranked <- ranked[!duplicated(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)

addm("GSEA_ranked_genes", length(ranked))

set.seed(42)
gsea_go <- gseGO(geneList = ranked, OrgDb = org.Hs.eg.db, ont = "BP",
                 pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500,
                 verbose = FALSE)
go_res <- gsea_go@result
write.csv(go_res, out(paste0("GSEA_GO_BP", ".csv")), row.names = FALSE)
n_go <- sum(go_res$p.adjust < 0.05)

addm("GSEA_GO_BP_padj_lt_0.05", n_go)
addm("GSEA_GO_BP_up_in_RA", sum(go_res$p.adjust < 0.05 & go_res$NES > 0))
addm("GSEA_GO_BP_down_in_RA", sum(go_res$p.adjust < 0.05 & go_res$NES < 0))

run_kegg <- function() tryCatch(
  gseKEGG(geneList = ranked, organism = "hsa", pvalueCutoff = 0.05,
          pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, verbose = FALSE),
  error = function(err) NULL)
gsea_kegg <- run_kegg()          # RNG continues from gseGO
if (!is.null(gsea_kegg)) {
  kg_res <- gsea_kegg@result
  write.csv(kg_res, out(paste0("GSEA_KEGG", ".csv")), row.names = FALSE)
  n_kg <- sum(kg_res$p.adjust < 0.05)
  
  addm("GSEA_KEGG_padj_lt_0.05", n_kg)
  addm("GSEA_KEGG_up_in_RA", sum(kg_res$p.adjust < 0.05 & kg_res$NES > 0))
  addm("GSEA_KEGG_down_in_RA", sum(kg_res$p.adjust < 0.05 & kg_res$NES < 0))
}

# ---- 4. Gene-set scores for OXPHOS and surfactant homeostasis are GSVA scores of GO BP sets (17_GSVA_GO_scores.R).

write.csv(metrics, out(paste0("metrics", ".csv")), row.names = FALSE)
print(metrics, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out(paste0("sessionInfo", ".txt")))

