# 00_cohort_config.R -- cohort table, output folders and column renames; sourced by every script.
# Inputs: article/cohort_definition.csv (Sample_ID, group, status, public_id).
ART_DIR    <- Sys.getenv("ART_DIR", "article")
ART_TABLES <- file.path(ART_DIR, "tables")
ART_FIGS   <- file.path(ART_DIR, "figures")

COHORT_FILE <- file.path(ART_DIR, "cohort_definition.csv")
if (!file.exists(COHORT_FILE)) stop("controlled-access cohort table not found: ", COHORT_FILE, " (see data_examples/README.md)")
COHORT_TAB <- read.csv(COHORT_FILE, stringsAsFactors = FALSE)
stopifnot(all(c("Sample_ID", "group", "status", "public_id") %in% names(COHORT_TAB)), all(COHORT_TAB$group %in% c("RA", "Sarcoidosis")),
          all(COHORT_TAB$status %in% c("analysed", "excluded_dx", "not_in_cohort")))
EXCL_NOT_IN_COHORT <- COHORT_TAB$Sample_ID[COHORT_TAB$status == "not_in_cohort"]   # sequenced, but not part of the study cohort (never analysed)
EXCL_DX <- COHORT_TAB$Sample_ID[COHORT_TAB$status == "excluded_dx"]           # enrolled as sarcoidosis, not analysed (other diagnosis)
EXCL_ALL  <- c(EXCL_NOT_IN_COHORT, EXCL_DX)
sample_group <- function(ids) COHORT_TAB$group[match(ids, COHORT_TAB$Sample_ID)]          # "RA" / "Sarcoidosis"
is_sarc      <- function(ids) sample_group(ids) == "Sarcoidosis"
public_label <- function(ids) COHORT_TAB$public_id[match(ids, COHORT_TAB$Sample_ID)]      # de-identified label (RA1, Sarcoidosis1, ...)

set_cohort <- function(cohort = "n31") {   # the 31 analysed patients
  excl <- EXCL_ALL; n_exp <- c(RA = 24, Sarc = 7)
  dir.create(file.path(ART_TABLES, cohort), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ART_FIGS, cohort), recursive = TRUE, showWarnings = FALSE)
  list(name = cohort, EXCL = excl, n = n_exp,
       tables = file.path(ART_TABLES, cohort), figs = file.path(ART_FIGS, cohort))
}

check_n <- function(ids, cfg, what = "") {
  ra <- sum(sample_group(ids) == "RA", na.rm = TRUE); sa <- sum(is_sarc(ids), na.rm = TRUE)
  if (ra != cfg$n[["RA"]] || sa != cfg$n[["Sarc"]])
    stop(sprintf("%s: RA=%d Sarc=%d, expected %d/%d", what, ra, sa, cfg$n[["RA"]], cfg$n[["Sarc"]]))
  if (any(ids %in% cfg$EXCL)) stop(what, ": excluded sample present")
  invisible(TRUE)
}

# Medication at the time of BAL (controlled-access file <ART_DIR>/tables/n31/medication_at_BAL.csv; RA patients only): a drug is in use
#   at BAL when its start date is on or before the BAL date (chart notes; MTX, TAC, SSZ). GC (glucocorticoid) is the original baseline flag of
#   the clinical records (different collection basis). Columns: Sample_ID, MTX (methotrexate), GC (glucocorticoid),
#   TAC (tacrolimus), SSZ (sulfasalazine), bts (biologic or targeted synthetic DMARD), IS_csDMARD (= TAC or SSZ: immunosuppressant or
#   csDMARD other than methotrexate; the article's definition also lists bucillamine and iguratimod, which no patient was taking at BAL).
#   0/1, no dates.
med_at_BAL <- function(ids, var) {
  m <- read.csv(file.path(ART_DIR, "tables", "n31", "medication_at_BAL.csv"), stringsAsFactors = FALSE)
  stopifnot(all(c("Sample_ID", "MTX", "GC", "TAC", "SSZ", "bts", "IS_csDMARD") %in% names(m)), nrow(m) == 24,
            sum(m$MTX) == 8, sum(m$GC) == 7, sum(m$IS_csDMARD) == 4, all(m$IS_csDMARD == as.integer(m$TAC == 1 | m$SSZ == 1)),
            all(ids %in% m$Sample_ID), var %in% names(m))
  as.numeric(m[[var]][match(ids, m$Sample_ID)])
}
# The controlled-access tables are shared with English column names (data_examples/README.md, "Column names"). A table that still
#   carries a non-English column name is rejected.
check_english_columns <- function(df) {
  bad <- names(df)[grepl("[^ -~]", names(df))]
  if (length(bad)) stop(length(bad), " non-English column name(s); use the English names listed in data_examples/README.md")
  df
}
