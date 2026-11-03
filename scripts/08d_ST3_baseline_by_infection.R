# 08d_ST3_baseline_by_infection.R -- Supplementary Table 3: clinical characteristics at BAL by subsequent respiratory infection (RA, n = 24;
#   medication at BAL from medication_at_BAL.csv). Median [IQR] with Wilcoxon rank-sum P (exact unless tied); n/N (%) with Fisher's exact P.
# Produces: ST3_by_infection.tsv.
# Inputs: planA_data_RA24.csv, medication_at_BAL.csv.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
R0 <- file.path(ART_TABLES, "n31")
d <- check_english_columns(read.csv(file.path(R0, "planA/planA_data_RA24.csv"), check.names = FALSE))
stopifnot(nrow(d) == 24, sum(d$respiratory_infection) == 6)
inf <- d$respiratory_infection; num <- function(x) suppressWarnings(as.numeric(x))
med <- function(v) med_at_BAL(d$Sample_ID, v)
cont <- function(nm, x) { q <- function(v) sprintf("%.1f [%.1f–%.1f]", median(v), quantile(v, .25), quantile(v, .75))
  sprintf("%s\t%s\t%s\t%.3f", nm, q(x[inf == 1]), q(x[inf == 0]), suppressWarnings(wilcox.test(x[inf == 1], x[inf == 0], exact = TRUE))$p.value) }
bin <- function(nm, x) { a <- sum(x[inf == 1]); b <- sum(x[inf == 0])
  p <- if (length(unique(x)) > 1) sprintf("%.3f", fisher.test(table(factor(x, 0:1), factor(inf, 0:1)))$p.value) else "—"
  sprintf("%s\t%d/6 (%.0f%%)\t%d/18 (%.0f%%)\t%s", nm, a, 100 * a / 6, b, 100 * b / 18, p) }
out <- c(cont("Age (years)", num(d$Age)), bin("Sex (male)", num(d$gender_Male.1._0.)), bin("Ever smoker", num(d$smoking.1.0.)),
  cont("Body mass index (kg/m²)", num(d$BMI)), cont("RA duration (years)", num(d$Years_since_diagnosis_of_RA)),
  cont("DAS28-CRP", num(d$DAS28.CRP)), cont("CDAI", num(d$CDAI)), cont("CRP (mg/L)", 10 * num(d$CRP)),   # CRP recorded in mg/dL
  bin("MTX use", med("MTX")), bin("Glucocorticoid use", med("GC")), bin("Immunosuppressant or csDMARD other than MTX", med("IS_csDMARD")),
  bin("b/tsDMARD use", med("bts")), bin("NSAIDs use", num(d$NSAIDS)), bin("Pneumococcal vaccination", num(d$pneumococcal.vaccination.1.0.)),
  bin("ILD at baseline (modified Ichikado score > 0)", as.integer(d$ILD_Group == "ILD_Positive")),
  bin("Rheumatoid factor positive", d$RF_positive), bin("Anti-CCP antibody positive", d$anti_CCP_positive))
writeLines(c("Variable\tInfection (+) (n = 6)\tInfection (−) (n = 18)\tP value", out), file.path(R0, "ST3_by_infection.tsv"))
cat(out, sep = "\n")
