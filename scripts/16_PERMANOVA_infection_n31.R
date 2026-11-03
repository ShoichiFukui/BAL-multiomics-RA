# 16_PERMANOVA_infection_n31.R -- Figure 3a PERMANOVA: transcriptome by subsequent respiratory infection within RA (n = 24, 6 events),
#   on the variance-stabilized matrix of the 31 analysed samples (DESeq2 vst, blind = FALSE, genes with >= 10 counts in >= 5 samples,
#   as for Figure 2a), RA columns only; Euclidean distance, vegan::adonis2, 999 permutations, set.seed(42). Manuscript: R2 = 0.091,
#   p = 0.027.
# Produces: <ART_TABLES>/n31/PERMANOVA_infection_RA24.csv (read by figures/panels_fig3a_3b_3h.R) and PERMANOVA_microbiome_infection_RA24.csv
#   (BAL microbiome by infection, Results) and microbiome_alpha_diversity_n31.csv.
# Inputs: results/RA_ILD_Workspace.RData, planA_data_RA24.csv (infection outcome; controlled access), the OTU table (project root).
# Infection status is read from planA_data_RA24.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressPackageStartupMessages({ library(DESeq2); library(vegan) })
cfg <- set_cohort("n31"); R0 <- file.path(ART_TABLES, "n31")
e <- new.env(); load("results/RA_ILD_Workspace.RData", envir = e)
ids <- e$common_samples[e$common_samples %in% setdiff(colnames(e$count_matrix), cfg$EXCL)]; check_n(ids, cfg, "permanova")
cm <- as.matrix(e$count_matrix[, ids]); cm <- cm[rowSums(cm) > 0, ]
grp <- factor(ifelse(is_sarc(ids), "Control", "RA"))
dds <- DESeqDataSetFromMatrix(cm, data.frame(Sample_Group = grp, row.names = ids), ~Sample_Group)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]; dds <- DESeq(dds, quiet = TRUE); vs <- assay(vst(dds, blind = FALSE))
pa <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"))
ra <- ids[grp == "RA"]; inf <- factor(as.numeric(pa$respiratory_infection[match(ra, pa$Sample_ID)])); stopifnot(length(ra) == 24, sum(inf == 1) == 6)
run <- function(m, lab) { set.seed(42); a <- adonis2(dist(t(m[, ra])) ~ inf, permutations = 999); a }
a31 <- run(vs, "31-sample VST (article)")
write.csv(data.frame(Matrix = "VST_31", R2 = a31$R2[1], F = a31$F[1], P = a31$`Pr(>F)`[1]),
          file.path(R0, "PERMANOVA_infection_RA24.csv"), row.names = FALSE)

# BAL microbiome composition by subsequent infection (Results: PERMANOVA R2 = 0.045, p = 0.34): the OTU rule of Supplementary Data 4
#   (18_SD4_microbiome_genus_exact.R): the 24 RA samples, OTUs with >= 5 reads in >= 3 of them, relative abundance per sample over the retained
#   OTUs; Bray-Curtis distances, adonis2, 999 permutations. The permutation P depends on the random-number stream of the installed
#   vegan/permute versions (vegan 2.7.2, permute 0.9.8 here; docs/package_versions.txt); R2 and F do not.
otu <- read.csv("OTU abundance table-2 (merged abundance table).csv", check.names = FALSE)
sc <- grep("Abundance$", colnames(otu)); sc <- sc[colnames(otu)[sc] != "Combined Abundance"]   # one column per sample (not the total column)
om <- as.matrix(otu[, sc]); colnames(om) <- sub("BALF$", "", sub("_S[0-9]+.*Abundance$", "", colnames(om)))
mi <- pa$Sample_ID; stopifnot(length(mi) == 24, all(mi %in% colnames(om)))
m24 <- om[, mi]; m24 <- m24[rowSums(m24 >= 5) >= 3, ]; rel <- sweep(m24, 2, colSums(m24), "/")
minf <- factor(ifelse(pa$respiratory_infection == 1, "Inf_pos", "Inf_neg"), levels = c("Inf_neg", "Inf_pos")); stopifnot(sum(minf == "Inf_pos") == 6)
set.seed(42); am <- adonis2(vegdist(t(rel), method = "bray") ~ minf, permutations = 999)
stopifnot(abs(am$R2[1] - 0.045) < 5e-4)   # value of the article
if (abs(am$`Pr(>F)`[1] - 0.34) > 0.005) warning("PERMANOVA P = ", am$`Pr(>F)`[1], " (article 0.34): the permutation P depends on the vegan/permute versions")
write.csv(data.frame(Matrix = "Microbiome_BrayCurtis_relative", OTUs = nrow(m24), R2 = am$R2[1], F = am$F[1], P = am$`Pr(>F)`[1]),
          file.path(R0, "PERMANOVA_microbiome_infection_RA24.csv"), row.names = FALSE)

# Alpha diversity (Shannon, natural log; Simpson 1 - D) of every sequenced sample from all OTUs: the values stored as diversity_df in
#   RA_ILD_Workspace.RData and used by 04 and 06 (Shannon by infection, p = 0.63, is computed in 06)
dvc <- data.frame(Sample_ID = colnames(om), Shannon = diversity(t(om), "shannon"), Simpson = diversity(t(om), "simpson"))
diversity_df <- e$diversity_df; chk <- diversity_df[diversity_df$Sample_ID %in% dvc$Sample_ID, ]
stopifnot(nrow(chk) >= 31, max(abs(dvc$Shannon[match(chk$Sample_ID, dvc$Sample_ID)] - chk$Shannon)) < 1e-9,
          max(abs(dvc$Simpson[match(chk$Sample_ID, dvc$Sample_ID)] - chk$Simpson)) < 1e-9)
write.csv(dvc[dvc$Sample_ID %in% COHORT_TAB$Sample_ID[COHORT_TAB$status == "analysed"], ], file.path(R0, "microbiome_alpha_diversity_n31.csv"), row.names = FALSE)
