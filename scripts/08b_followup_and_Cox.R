# 08b_followup_and_Cox.R -- follow-up for infection and the Cox model for BAL Th17.1 (Results: median follow-up, incidence,
#   Cox hazard ratio per 10 points of BAL Th17.1 with profile-likelihood CI and likelihood-ratio P; Supplementary Methods).
# Produces: followup_cutoff.txt.
# Inputs: planA_data_RA24.csv, followup_by_sample.csv, infection_events.csv (days from BAL).
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
source("scripts/00_cohort_config.R")
suppressMessages(library(survival))
R0 <- file.path(ART_TABLES, "n31"); TXT <- file.path(R0, "followup_cutoff.txt")
d <- read.csv(file.path(R0, "planA/planA_data_RA24.csv"), check.names = FALSE)
fu <- read.csv(file.path(R0, "followup_by_sample.csv"), stringsAsFactors = FALSE)
evt <- read.csv(file.path(R0, "infection_events.csv"), stringsAsFactors = FALSE); ev <- setNames(evt$days_to_infection, evt$Sample_ID); stopifnot(length(ev) == 6)
d <- merge(d, fu, by = "Sample_ID"); d$inf <- as.numeric(ev[d$Sample_ID])   # days from BAL to infection (NA = no infection)
th <- d[["BALF_Th17.1_CCR6pos_CXCR3pos"]]
sink(TXT)
run <- function(label, end, status) {
  t_end <- ifelse(is.na(d$inf), end, d$inf); e <- as.numeric(!is.na(d$inf))   # end = days from BAL to the last visit
  yr <- t_end / 365.25
  cat("\n==", label, "==\n")
  cat("events", sum(e), "| censored at death", sum(status == "death"), "| censored at last visit before cut-off (lost/transferred)", sum(status == "lost"), "| followed to cut-off", sum(status == "continuing"), "\n")
  cat("time to infection or censoring, all 24: median", round(median(yr), 2), "y, range", round(min(yr), 2), "-", round(max(yr), 2), "\n")
  cat("  censored patients only (n =", sum(e == 0), "): median", round(median(yr[e == 0]), 2), "range", round(min(yr[e == 0]), 2), "-", round(max(yr[e == 0]), 2), "\n")
  cat("person-years", round(sum(yr), 1), "| incidence per 100 PY", round(sum(e) / sum(yr) * 100, 1), "\n")
  grp <- th <= median(th)
  cat("events BAL Th17.1 <= median vs >:", sum(e[grp]), "/", sum(grp), "vs", sum(e[!grp]), "/", sum(!grp), "| log-rank p", signif(1 - pchisq(survdiff(Surv(t_end, e) ~ grp)$chisq, 1), 3), "\n")
  s <- summary(coxph(Surv(t_end, e) ~ I(th / 10)))
  cat("Cox HR per 10-point BAL Th17.1", round(s$conf.int[1], 2), "(", round(s$conf.int[3], 2), "-", round(s$conf.int[4], 2), ") Wald p", signif(s$coefficients[5], 3), "LR p", signif(s$logtest[3], 3), "C", round(s$concordance[1], 3), "\n")
  km <- summary(survfit(Surv(yr, e) ~ 1), times = c(1, 3, 5), extend = TRUE); cat("KM infection-free 1/3/5 y:", round(km$surv, 3), "\n")
  cat("longest censoring time (days):", max(end[is.na(d$inf)]), "| follow-up CT after infection censoring:", sum(as.numeric(d$CT2d) > end & is.na(d$inf), na.rm = TRUE), "patients\n")
}
cat("Existing (manuscript, days from planA): median", round(median(d$days_to_follow_up / 365.25), 2), "range", round(range(d$days_to_follow_up / 365.25), 2), "\n")
run("A: first census", d$days_to_last_visit_A, d$status_A)
run("B: last census", d$days_to_last_visit_B, d$status_B)
cat("\nPer-patient (days from BAL: A end, B end, status):\n"); print(d[order(d$Sample_ID), c("Sample_ID", "Infection", "days_to_last_visit_A", "status_A", "days_to_last_visit_B", "status_B")], row.names = FALSE)
sink()

sink(TXT, append = TRUE)
cat("\n== Reverse Kaplan-Meier median follow-up (censoring = event) ==\n")
for (sc in list(list("A: first census", d$days_to_last_visit_A), list("B: last census", d$days_to_last_visit_B))) {
  end <- sc[[2]]; e <- as.numeric(!is.na(d$inf)); t <- ifelse(e == 1, d$inf, end) / 365.25
  rk <- survfit(Surv(t, 1 - e) ~ 1); cat(sc[[1]], ": median", round(summary(rk)$table["median"], 2), "y (range of censoring times", round(min(t[e == 0]), 2), "-", round(max(t[e == 0]), 2), ")\n")
}
sink()

sink(TXT, append = TRUE)
cat("\n== Cox HR per 10 points of BAL Th17.1, scenario B: profile-likelihood CI ==\n")
d$x <- d[["BALF_Th17.1_CCR6pos_CXCR3pos"]] / 10
d$t <- ifelse(is.na(d$inf), d$days_to_last_visit_B, d$inf); d$e <- as.numeric(!is.na(d$inf))
f <- coxph(Surv(t, e) ~ x, d); ll0 <- f$loglik[2]
pl <- function(b) coxph(Surv(t, e) ~ offset(b * x), d)$loglik
g <- function(b) 2 * (ll0 - pl(b)) - qchisq(.95, 1)
lo <- uniroot(g, c(coef(f) - 6, coef(f)))$root; hi <- uniroot(g, c(coef(f), coef(f) + 3))$root
cat("HR", round(exp(coef(f)), 2), "Wald CI", round(exp(confint(f)), 2), "| profile CI", round(exp(c(lo, hi)), 2), "| LR p", signif(summary(f)$logtest[3], 3), "Wald p", signif(summary(f)$coef[5], 3), "\n")
if (requireNamespace("coxphf", quietly = TRUE)) { ff <- coxphf::coxphf(Surv(t, e) ~ x, d); cat("Firth HR", round(exp(ff$coef), 2), "CI", round(ff$ci.lower, 2), round(ff$ci.upper, 2), "p", signif(ff$prob, 3), "\n") } else cat("coxphf not installed\n")
sink()
