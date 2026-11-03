# 08c_SD2_effect_sizes.R -- Supplementary Data 2 effect sizes (116 variables, including BAL lymphocyte % of the cell differential,
#   BH over the 110 testable). Run after 08e.
# Produces: SD2_effect_sizes_117.csv, supp_tables/SD2_effect_sizes_116.csv.
# Inputs: results/tables/master_data_with_CT_all.csv, planA_data_RA24.csv, SD2_effect_sizes_RA_vs_Sarc.csv (08e).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
R0 <- file.path(ART_TABLES, "n31")
d <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE))
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"), check.names = FALSE)$Sample_ID
sarc <- COHORT_TAB$Sample_ID[COHORT_TAB$group == "Sarcoidosis" & COHORT_TAB$status == "analysed"]   # the seven analysed sarcoidosis patients
d <- d[d$Sample_ID %in% c(ra, sarc), ]; g <- ifelse(d$Sample_ID %in% ra, "RA", "Sarc"); stopifnot(table(g)[["RA"]] == 24, table(g)[["Sarc"]] == 7)
es <- read.csv(file.path(R0, "SD2_effect_sizes_RA_vs_Sarc.csv"), stringsAsFactors = FALSE); stopifnot(nrow(es) == 115, !any(grepl("sCD40L", es$Variable)))
for (v in c("Serum_IL6", "BALF_IL6")) { x <- as.numeric(d[[v]]); p <- wilcox.test(x[g == "RA"], x[g != "RA"], exact = TRUE)$p.value; stopifnot(abs(p - es$p_value[es$Variable == v]) < 1e-6) }
mag <- function(a) { a <- abs(a); if (a < 0.147) "negligible" else if (a < 0.33) "small" else if (a < 0.474) "medium" else "large" }
new <- do.call(rbind, lapply(c("Serum_sCD40L", "BALF_sCD40L"), function(v) {
  x <- as.numeric(d[[v]]); a <- x[g == "RA"]; b <- x[g != "RA"]
  p <- suppressWarnings(wilcox.test(a, b, exact = TRUE)$p.value)
  cd <- mean(outer(a, b, ">")) - mean(outer(a, b, "<"))
  data.frame(Variable = v, n1 = length(a), n2 = length(b), median1 = median(a), median2 = median(b), p_value = p, cliff_delta = cd,
             magnitude = mag(cd), Category = ifelse(grepl("^Serum", v), "Serum_Cytokine", "BALF_Cytokine"), p_adjusted = NA)
}))
print(new)
es2 <- rbind(es, new); es2$p_adjusted <- p.adjust(es2$p_value, method = "BH")

write.csv(es2, file.path(R0, "SD2_effect_sizes_117.csv"), row.names = FALSE)

DROP <- c("BALF_Neutrophil_CD66bpos_CD16pos", "PB_Neutrophil_CD66bpos_CD16pos"); stopifnot(all(DROP %in% es2$Variable))
OD7 <- file.path(R0, "supp_tables"); dir.create(OD7, recursive = TRUE, showWarnings = FALSE)
# BAL lymphocyte % (cell differential; BALF_Lymphocyte_Percent, the same source as Table 1), computed as the other rows
x <- as.numeric(d$BALF_Lymphocyte_Percent); a <- x[g == "RA" & !is.na(x)]; b <- x[g != "RA" & !is.na(x)]
cdl <- mean(sign(outer(a, b, "-")))
lym <- data.frame(Variable = "BALF_Lymphocyte_Percent", n1 = length(a), n2 = length(b), median1 = median(a), median2 = median(b),
                  p_value = suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value, cliff_delta = cdl, magnitude = mag(cdl),
                  Category = "BAL_differential_count", p_adjusted = NA)
stopifnot(!"BALF_Lymphocyte_Percent" %in% es2$Variable, lym$n1 == 24, lym$n2 == 7, abs(lym$p_value - 0.478) < 0.0015)   # = Table 1
es7 <- rbind(es2[!es2$Variable %in% DROP, ], lym); es7$p_adjusted <- p.adjust(es7$p_value, method = "BH")
stopifnot(nrow(es7) == 116, sum(!is.na(es7$p_value)) == 110,
          sum(es7$p_adjusted < 0.05, na.rm = TRUE) == 32,   # Results: 32 of the 116 variables differ after FDR correction (Fig. 2e)
          abs(es7$p_adjusted[es7$Category == "Deconvolution" & es7$Variable == "T_cell"] - 0.0845) < 5e-4)   # Supplementary Data 2, BayesPrism T cell

write.csv(es7, file.path(OD7, "SD2_effect_sizes_116.csv"), row.names = FALSE)
write.csv(es7[, names(es)], file.path(OD7, "SD2_effect_sizes_RA_vs_Sarc.csv"), row.names = FALSE)
# reference row of Supplementary Data 2, outside the screen and its BH family: BAL CD4/CD8 ratio of the clinical lymphocyte subsets
x <- as.numeric(d$CD4_CD8_Ratio); a <- x[g == "RA" & !is.na(x)]; b <- x[g != "RA" & !is.na(x)]; cdr <- mean(sign(outer(a, b, "-")))
ref <- data.frame(Variable = "BAL_CD4_CD8_Ratio", n1 = length(a), n2 = length(b), median1 = median(a), median2 = median(b),
                  p_value = suppressWarnings(wilcox.test(a, b, exact = TRUE))$p.value, cliff_delta = cdr, magnitude = mag(cdr))
print(ref); write.csv(ref, file.path(OD7, "SD2_reference_CD4_CD8_ratio.csv"), row.names = FALSE)
