# 08_sensitivity_AUC_stability_Firth.R -- infection analyses of the BAL Th17.1 proportion and the OXPHOS score (RA, n = 24):
#   AUC stability (bootstrap, jackknife, optimism, permutation; Results, methods in the Supplementary Methods), OXPHOS GSVA score adjusted
#   for the BAL differential count (Results), Firth adjustment one covariate at a time (Supplementary Table 5; medication at BAL from
#   medication_at_BAL.csv), the exact and permutation P of the BAL Th17.1 AUC (Results) and the sample size of a validation study.
# Produces (in <ART_TABLES>/n31/sensitivity/): AUC_stability.csv, OXPHOS_adjusted_for_BAL_differential.csv, ST5_Firth_adjusted.csv,
#   AUC_permutation_two_sided.csv, validation_sample_size.csv (Discussion).
# Inputs: results/tables/master_data_with_CT_all.csv, planA_data_RA24.csv, gsva_GO_scores_n31.csv (17).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")

suppressPackageStartupMessages({
  library(pROC); library(logistf); library(ppcor)
})

OUTDIR <- file.path(ART_TABLES, "n31", "sensitivity")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- load data --

master <- check_english_columns(read.csv("results/tables/master_data_with_CT_all.csv", check.names = FALSE))   # per-patient table, 35 rows

stopifnot(nrow(master) == 35)
gv <- read.csv(file.path(ART_TABLES, "n31", "gsva_GO/gsva_GO_scores_n31.csv"))   # 17_GSVA_GO_scores.R: GSVA of GO BP OXPHOS / surfactant
master$OXPHOS_score     <- gv$OXPHOS[match(master$Sample_ID, gv$Sample_ID)]       # (GSVA GO score)
master$Surfactant_score <- gv$Surfactant[match(master$Sample_ID, gv$Sample_ID)]

# ------------------------------------------------------------ analysis sets --
RA <- master[master$Sample_Group == "RA", ]

pa <- read.csv(file.path(ART_TABLES, "n31", "planA/planA_data_RA24.csv"))   # infection outcome (extended follow-up)
stopifnot(nrow(RA) == 24, all(RA$Sample_ID %in% pa$Sample_ID))
RA$infection <- as.integer(pa$respiratory_infection[match(RA$Sample_ID, pa$Sample_ID)]); stopifnot(sum(RA$infection) == 6)
for (v in c("MTX", "GC", "IS_csDMARD", "bts")) RA[[paste0(v, "_atBAL")]] <- med_at_BAL(RA$Sample_ID, v)

RA$BAL_Th17_1 <- as.numeric(RA[["BALF_Th17.1_CCR6pos_CXCR3pos"]])
RA$PB_Th17    <- as.numeric(RA[["PB_Th17_CCR6pos_CXCR3neg"]])

auc_stability <- function(pred, outcome, label, B = 2000) {
  ok <- !is.na(pred) & !is.na(outcome)
  x <- pred[ok]; y <- outcome[ok]

  dir0 <- if (mean(x[y == 1]) < mean(x[y == 0])) ">" else "<"
  r <- roc(y, x, quiet = TRUE, direction = dir0)
  apparent <- as.numeric(auc(r))

  set.seed(42)
  ci_b <- as.numeric(ci.auc(r, method = "bootstrap", boot.n = B))

  jk <- vapply(seq_along(y), function(i) {
    if (length(unique(y[-i])) < 2) return(NA_real_)
    as.numeric(auc(roc(y[-i], x[-i], quiet = TRUE, direction = dir0)))
  }, numeric(1))

  set.seed(42)
  opt <- replicate(B, {
    idx <- sample(seq_along(y), replace = TRUE)
    if (length(unique(y[idx])) < 2) return(NA_real_)
    fit <- tryCatch(suppressWarnings(glm(y[idx] ~ x[idx], family = binomial)),
                    error = function(e) NULL)
    if (is.null(fit) || any(is.na(coef(fit)))) return(NA_real_)
    b      <- coef(fit)
    lp_boot <- b[1] + b[2] * x[idx]
    lp_orig <- b[1] + b[2] * x
    a_boot <- as.numeric(auc(roc(y[idx], as.numeric(lp_boot), quiet = TRUE, direction = "<")))
    a_orig <- as.numeric(auc(roc(y,      as.numeric(lp_orig), quiet = TRUE, direction = "<")))
    a_boot - a_orig
  })
  optimism  <- mean(opt, na.rm = TRUE)
  corrected <- apparent - optimism

  set.seed(42)
  perm <- replicate(B, as.numeric(auc(roc(sample(y), x, quiet = TRUE, direction = dir0))))
  p_perm <- (sum(perm >= apparent) + 1) / (B + 1)

  data.frame(
    Marker           = label,
    N                = length(y),
    Events           = sum(y == 1),
    AUC_apparent     = round(apparent, 3),
    Boot_CI_low      = round(ci_b[1], 3),
    Boot_CI_high     = round(ci_b[3], 3),
    Boot_CI_width    = round(ci_b[3] - ci_b[1], 3),
    Jackknife_min    = round(min(jk, na.rm = TRUE), 3),
    Jackknife_max    = round(max(jk, na.rm = TRUE), 3),
    Jackknife_range  = round(max(jk, na.rm = TRUE) - min(jk, na.rm = TRUE), 3),
    Optimism         = round(optimism, 3),
    AUC_corrected    = round(corrected, 3),
    Permutation_P    = round(p_perm, 4),
    stringsAsFactors = FALSE
  )
}

auc_tab <- auc_stability(RA$BAL_Th17_1, RA$infection,  "BAL Th17.1 -> respiratory infection")   # AUC stability, quoted in the Results
print(auc_tab, row.names = FALSE)
write.csv(auc_tab, file.path(OUTDIR, "AUC_stability.csv"), row.names = FALSE)

adj_composition <- function(score, label) {
  d <- data.frame(
    y    = RA[[score]],
    inf  = RA$infection,
    mac  = as.numeric(RA[["BAL_macrophage_pct"]]),
    lym  = as.numeric(RA[["BAL_lymphocyte_pct"]]),
    neu  = as.numeric(RA[["BAL_neutrophil_pct"]])
  )
  d <- d[complete.cases(d), ]
  w  <- suppressWarnings(wilcox.test(y ~ inf, data = d, exact = TRUE))
  f  <- lm(y ~ inf + mac + lym + neu, data = d)
  s  <- summary(f)$coefficients["inf", ]
  data.frame(Score = label, N = nrow(d),
             Unadj_Wilcoxon_P = round(w$p.value, 4),
             Adj_beta_infection = round(s[["Estimate"]], 3),
             Adj_SE  = round(s[["Std. Error"]], 3),
             Adj_P   = round(s[["Pr(>|t|)"]], 4),
             stringsAsFactors = FALSE)
}
res_adj <- rbind(
  adj_composition("OXPHOS_score",     "OXPHOS"),
  adj_composition("Surfactant_score", "Surfactant")
)
stopifnot(abs(res_adj$Adj_P[res_adj$Score == "OXPHOS"] - 0.0279) < 1e-3, abs(res_adj$Unadj_Wilcoxon_P[res_adj$Score == "OXPHOS"] - 0.0094) < 1e-3)   # values of the article
print(res_adj, row.names = FALSE)
write.csv(res_adj, file.path(OUTDIR, "OXPHOS_adjusted_for_BAL_differential.csv"), row.names = FALSE)

getcol <- function(nm) RA[[nm]]

confounders <- list(
  Age                = "Age",
  Male               = "gender_Male.1._0.",
  Smoking_cat        = "non._ex._current_smoking.0.1.2.",   # never / former / current (0/1/2)
  Steroid            = "GC_atBAL",            # medication at BAL (medication_at_BAL.csv)
  IS_csDMARD         = "IS_csDMARD_atBAL",    # immunosuppressant or csDMARD other than methotrexate
  Biologics          = "bts_atBAL",           # all zero at BAL -> skipped
  MTX                = "MTX_atBAL",
  Baseline_ILD       = "ILD_Group",
  DAS28_CRP          = "DAS28.CRP",
  BMI                = "BMI",
  Pneumococcal_vacc  = "pneumococcal.vaccination.1.0."
)

firth_one <- function(cname, cvar) {
  cv <- getcol(cvar)
  if (is.null(cv)) stop("covariate column not found: ", cvar)
  if (is.character(cv) || is.factor(cv)) cv <- as.integer(factor(cv)) - 1L
  cv <- suppressWarnings(as.numeric(cv))
  stopifnot(!all(is.na(cv)))
  d <- data.frame(y = RA$infection, x = RA$BAL_Th17_1, z = as.numeric(cv))
  d <- d[complete.cases(d), ]
  if (length(unique(d$z)) < 2) return(NULL)
  fit <- tryCatch(logistf(y ~ x + z, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  i <- which(names(coef(fit)) == "x")
  data.frame(Confounder = cname, N = nrow(d),
             OR_per_1pct_Th17.1 = exp(coef(fit)[i]),
             CI_low  = exp(fit$ci.lower[i]),
             CI_high = exp(fit$ci.upper[i]),
             P       = fit$prob[i],
             stringsAsFactors = FALSE)
}
d0  <- data.frame(y = RA$infection, x = RA$BAL_Th17_1); d0 <- d0[complete.cases(d0), ]
f0  <- logistf(y ~ x, data = d0)
res_firth <- rbind(
  data.frame(Confounder = "(unadjusted)", N = nrow(d0),
             OR_per_1pct_Th17.1 = exp(coef(f0)[2]),
             CI_low = exp(f0$ci.lower[2]),
             CI_high = exp(f0$ci.upper[2]),
             P = f0$prob[2], stringsAsFactors = FALSE),
  do.call(rbind, lapply(names(confounders),
                        function(k) firth_one(k, confounders[[k]])))
)
stopifnot(abs(res_firth$P[res_firth$Confounder == "MTX"] - 0.0634) < 5e-4, abs(res_firth$P[res_firth$Confounder == "Steroid"] - 0.0272) < 5e-4, !"Biologics" %in% res_firth$Confounder, nrow(res_firth) == 11)   # medication at BAL; 11 rows of Supplementary Table 5
print(res_firth, row.names = FALSE)
write.csv(res_firth, file.path(OUTDIR, "ST5_Firth_adjusted.csv"), row.names = FALSE)

suppressPackageStartupMessages(library(pROC))
R0 <- file.path(ART_TABLES, "n31")
ra <- read.csv(file.path(R0, "planA/planA_data_RA24.csv")); x <- ra[["BALF_Th17.1_CCR6pos_CXCR3pos"]]; y <- ra$respiratory_infection; stopifnot(sum(y) == 6, length(y) == 24)
dir0 <- ">"
apparent <- as.numeric(auc(roc(y, x, quiet = TRUE, direction = dir0))); B <- 2000
set.seed(42); perm <- replicate(B, as.numeric(auc(roc(sample(y), x, quiet = TRUE, direction = dir0))))
p1 <- (sum(perm >= apparent) + 1) / (B + 1); p2 <- (sum(abs(perm - 0.5) >= abs(apparent - 0.5) - 1e-12) + 1) / (B + 1)
r <- rank(x); cb <- combn(24, 6); W <- colSums(matrix(r[cb], nrow = 6)); mu <- 6 * 25 / 2; obs <- sum(r[y == 1])
pe2 <- mean(abs(W - mu) >= abs(obs - mu) - 1e-9); pe1 <- mean(W <= obs + 1e-9)
st4 <- data.frame(AUC = apparent, perm_one_sided_2000 = p1, perm_two_sided_2000 = p2, exact_one_sided = pe1, exact_two_sided = pe2, wilcox_exact = suppressWarnings(wilcox.test(x[y == 1], x[y == 0], exact = TRUE))$p.value)
print(st4, digits = 4); write.csv(st4, file.path(OUTDIR, "AUC_permutation_two_sided.csv"), row.names = FALSE)
stopifnot(abs(p1 - 0.004) < 0.0015)

# Sample size for a validation study (Discussion): events needed to detect AUC 0.70 with 80% power (two-sided alpha 0.05) at the
#   control-to-case ratio of this cohort (18:6 = 3), and the number of patients at 25% incidence
pw <- pROC::power.roc.test(auc = 0.70, power = 0.80, sig.level = 0.05, kappa = 3)
ss <- data.frame(AUC = 0.70, power = 0.80, kappa = 3, events = ceiling(pw$ncases), patients = ceiling(pw$ncases + pw$ncontrols))   # kappa 3 = 25% incidence
print(ss); write.csv(ss, file.path(OUTDIR, "validation_sample_size.csv"), row.names = FALSE)
stopifnot(ss$events == 21, ss$patients == 83)   # values of the article
