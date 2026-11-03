# 18_SD4_microbiome_genus_exact.R -- Supplementary Data 4 (BAL microbiome genera by subsequent respiratory infection, RA n = 24
#   infection n = 6) and Supplementary Fig. 3b.
# Rule: the 24 RA samples with an infection outcome; an OTU is kept if it has >= 5 reads in >= 3 of these 24 samples; genus = 6th level of
#   the taxonomy (otherwise "Unclassified"); reads summed per genus; relative abundance (%) per sample over the retained OTUs (106 genera)
#   no genus prevalence filter; two-sided Wilcoxon rank-sum test with exact = TRUE (exact where there are no ties, otherwise the normal
#   approximation with continuity correction); BH adjustment across the 106 genera; log2 fold change = log2((mean+ + 0.001)/(mean- + 0.001)).
#   Also writes per genus the number of samples with zero abundance and the group medians (most genera are absent from many samples,
#   which is why Supplementary Data 4 reports means).
# Produces: <ART_TABLES>/n31/microbiome/SD4_microbiome_genus_exact.csv and SD4_zero_median.csv.
# Inputs: "OTU abundance table-2 (merged abundance table).csv" (project root) and
#   planA/planA_data_RA24.csv (Sample_ID, respiratory_infection; controlled access).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
OUT <- file.path(ART_TABLES, "n31", "microbiome"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ra <- read.csv(file.path(ART_TABLES, "n31", "planA/planA_data_RA24.csv"))[, c("Sample_ID", "respiratory_infection")]
stopifnot(nrow(ra) == 24, all(sample_group(ra$Sample_ID) == "RA"), sum(ra$respiratory_infection == 1) == 6)
d <- read.csv("OTU abundance table-2 (merged abundance table).csv", check.names = FALSE)
sid <- sub("BALF$", "", sub("_S[0-9]+.*Abundance$", "", colnames(d)))   # column "<Sample_ID>[BALF]_S<n> Abundance"
k <- match(ra$Sample_ID, sid); stopifnot(!anyNA(k))
m <- as.matrix(d[, k]); colnames(m) <- ra$Sample_ID
keep <- rowSums(m >= 5) >= 3; m <- m[keep, ]
g <- vapply(strsplit(d$Taxonomy[keep], ";\\s*"), function(x) if (length(x) >= 6) trimws(x[6]) else "Unclassified", "")
gm <- rowsum(m, g); rel <- sweep(gm, 2, colSums(gm), "/") * 100
pos <- ra$respiratory_infection == 1
r <- data.frame(Genus = rownames(rel), Mean_InfPos = rowMeans(rel[, pos]), Mean_InfNeg = rowMeans(rel[, !pos]))
r$Log2FC <- log2((r$Mean_InfPos + 0.001) / (r$Mean_InfNeg + 0.001))
r$P_value <- apply(rel, 1, function(v) suppressWarnings(wilcox.test(v[pos], v[!pos], exact = TRUE))$p.value)
r$P_adjusted <- p.adjust(r$P_value, "BH")
r <- r[order(r$P_value), ]; stopifnot(nrow(r) == 106)
write.csv(r, file.path(OUT, "SD4_microbiome_genus_exact.csv"), row.names = FALSE)
z <- data.frame(Genus = rownames(rel), Zeros_of_24 = rowSums(rel == 0), Median_InfPos = apply(rel[, pos], 1, median), Median_InfNeg = apply(rel[, !pos], 1, median))
write.csv(z[match(r$Genus, z$Genus), ], file.path(OUT, "SD4_zero_median.csv"), row.names = FALSE)
n <- z[z$Genus != "Unclassified", ]

