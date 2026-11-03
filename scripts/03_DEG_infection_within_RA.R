# 03_DEG_infection_within_RA.R -- DESeq2 and GSEA, infection (+) versus (-) within RA (Fig. 3g; Suppl. Fig. 3a; results/DEG_infection_within_RA_n24_genome_wide.csv).
# Produces (in <ART_TABLES>/n31/): DEG_Infection_pos_vs_neg.csv, GSEA_Infection_GO_BP.csv, GSEA_Infection_KEGG.csv,
#   DEG_Infection_adjusted_BAL_differential.csv, DEG_Infection_sessionInfo.txt
#   released in results/ under the names listed in results/README.md.
# Inputs: results/RA_ILD_Workspace.RData, planA_data_RA24.csv, results/tables/master_data_with_CT_all.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(DESeq2); library(dplyr); library(clusterProfiler); library(org.Hs.eg.db) })
cfg <- set_cohort("n31"); R0 <- cfg$tables

e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
pa <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))
ra_common <- e$common_samples[sample_group(e$common_samples) %in% "RA"]; ra_common <- ra_common[!ra_common %in% cfg$EXCL]
stopifnot(length(ra_common) == 24, all(ra_common %in% pa$Sample_ID))
inf <- factor(ifelse(pa$respiratory_infection[match(ra_common, pa$Sample_ID)] == 1, "Inf_pos", "Inf_neg"), levels = c("Inf_neg", "Inf_pos"))
stopifnot(sum(inf == "Inf_pos") == 6)
count_inf <- as.matrix(e$count_matrix[, ra_common]); count_inf <- count_inf[rowSums(count_inf >= 10) >= 3, ]

dds_inf <- DESeqDataSetFromMatrix(countData = count_inf, colData = data.frame(row.names = ra_common, Infection = inf), design = ~Infection)
dds_inf <- DESeq(dds_inf, quiet = TRUE)
COEF <- "Infection_Inf_pos_vs_Inf_neg"
res0 <- lfcShrink(dds_inf, coef = COEF, type = "normal", quiet = TRUE)
rc <- assays(dds_inf)[["replaceCounts"]]

if (!is.null(rc)) { dds_r <- dds_inf; counts(dds_r) <- rc
  fx <- lfcShrink(dds_r, coef = COEF, type = "normal", quiet = TRUE)
  res0$log2FoldChange <- fx[rownames(res0), "log2FoldChange"]; res0$lfcSE <- fx[rownames(res0), "lfcSE"] }   # LFC sign agrees with the Wald statistic

deg <- as.data.frame(res0) %>% mutate(Gene_Symbol = rownames(res0)) %>%
  left_join(e$gene_annotations[, c("Gene_Symbol", "ENTREZID", "GENENAME")], by = "Gene_Symbol") %>% arrange(pvalue)
deg$Regulation <- "NS"   # adjusted P < 0.05 with |log2FC| <= 0.585 is labelled sub-threshold, as in Supplementary Data 1
deg$Regulation[deg$padj < 0.05 & deg$log2FoldChange > 0] <- "Up_in_InfPos (sub-threshold)"
deg$Regulation[deg$padj < 0.05 & deg$log2FoldChange < 0] <- "Down_in_InfPos (sub-threshold)"
deg$Regulation[deg$padj < 0.05 & deg$log2FoldChange > 0.585] <- "Up_in_InfPos"
deg$Regulation[deg$padj < 0.05 & deg$log2FoldChange < -0.585] <- "Down_in_InfPos"
write.csv(deg, file.path(R0, "DEG_Infection_pos_vs_neg.csv"), row.names = FALSE)
ok <- is.finite(deg$stat) & is.finite(deg$log2FoldChange)

sf <- c("SFTPB", "ETV5", "LPCAT1", "PGC", "SFTA3", "CLDN18", "SFTPD", "SFTPA1", "SFTPA2", "SFTA2", "SFTPC", "SLC34A2", "ABCA3")   # surfactant genes (Figure 3g)
print(data.frame(Gene = sf, LFC = round(deg$log2FoldChange[match(sf, deg$Gene_Symbol)], 3), stat = round(deg$stat[match(sf, deg$Gene_Symbol)], 2),
                 padj = signif(deg$padj[match(sf, deg$Gene_Symbol)], 3)), row.names = FALSE)

ranked <- deg$log2FoldChange; names(ranked) <- deg$ENTREZID
ranked <- ranked[!is.na(names(ranked)) & names(ranked) != ""]; ranked <- ranked[!duplicated(names(ranked))]; ranked <- sort(ranked, decreasing = TRUE)
set.seed(42)
ggo <- gseGO(geneList = ranked, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, verbose = FALSE)
write.csv(ggo@result, file.path(R0, "GSEA_Infection_GO_BP.csv"), row.names = FALSE)
g <- ggo@result; print(head(g[order(g$p.adjust), c("Description", "NES", "p.adjust")], 8), row.names = FALSE)
kk <- tryCatch(gseKEGG(geneList = ranked, organism = "hsa", pvalueCutoff = 0.1, pAdjustMethod = "BH", minGSSize = 15, maxGSSize = 500, verbose = FALSE),
               error = function(err) NULL)
if (!is.null(kk)) write.csv(kk@result, file.path(R0, "GSEA_Infection_KEGG.csv"), row.names = FALSE)
# Adjustment for the BAL differential count (Results; Supplementary Methods): the same DESeq2 model with the macrophage, lymphocyte and
#   neutrophil percentages (z-scored) as covariates, design ~ mac + lym + neu + Infection; BH-adjusted P of SFTPB, ETV5 and LPCAT1.
ms <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE)); ms <- ms[match(ra_common, ms$Sample_ID), ]
z <- function(v) as.numeric(scale(as.numeric(v)))
cd <- data.frame(row.names = ra_common, Infection = inf, mac = z(ms$BAL_macrophage_pct), lym = z(ms$BAL_lymphocyte_pct), neu = z(ms$BAL_neutrophil_pct))
dda <- DESeq(DESeqDataSetFromMatrix(count_inf, cd, ~ mac + lym + neu + Infection), quiet = TRUE)
ra_adj <- as.data.frame(results(dda, contrast = c("Infection", "Inf_pos", "Inf_neg")))[c("SFTPB", "ETV5", "LPCAT1"), c("log2FoldChange", "pvalue", "padj")]
stopifnot(abs(ra_adj["SFTPB", "padj"] - 0.052) < 5e-4, abs(ra_adj["ETV5", "padj"] - 0.189) < 5e-4)   # values of the article
write.csv(data.frame(Gene = rownames(ra_adj), ra_adj), file.path(R0, "DEG_Infection_adjusted_BAL_differential.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(R0, "DEG_Infection_sessionInfo.txt"))

