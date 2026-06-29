# ============================================================================
# 00_run_all.R
# End-to-end runner for the manuscript analysis pipeline.
# Sources each step in order. Heavy step 01 (BayesPrism) can be skipped if its
# outputs already exist.
#
# Usage:  Rscript scripts/00_run_all.R     (run from anywhere)
#
# All analysis scripts use paths relative to the PROJECT ROOT (the parent of
# scripts/), e.g. "output_v3_BayesPrism/..." and "results/...". This runner
# therefore sets the working directory to the project root and sources each
# step with its "scripts/" path, so the paths resolve identically whether a
# script is run individually or through this runner.
# ============================================================================

# Resolve scripts/ dir, then move to the project root (its parent).
this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
BASEDIR <- normalizePath(file.path(this_dir, ".."))
Sys.setenv(PROJECT_DIR = BASEDIR)   # consumed by 06_Comprehensive_Figures.R
setwd(BASEDIR)
cat(sprintf("Project root (working directory): %s\n", BASEDIR))

step <- function(label, path) {
  if (!file.exists(path)) stop(sprintf("Script not found: %s (run from the repository, expecting %s)", path, file.path(BASEDIR, path)))
  cat(sprintf("\n========== %s ==========\n", label))
  source(path, local = FALSE, echo = FALSE)
}

# --- Step 1: scRNA-seq atlas + BayesPrism (heavy; skip if already built) ---
if (!file.exists("output_v3_BayesPrism/celltype_proportions_cellFraction.csv")) {
  step("01 BayesPrism deconvolution", "scripts/01_BayesPrism_Deconvolution.R")
} else {
  cat("\nSkipping 01 (BayesPrism outputs already exist)\n")
}

# --- Step 2: master workspace ---
step("02 PostDeconvolution", "scripts/02_PostDeconvolution.R")

# --- Step 3: primary analyses (DEG, GSEA, GSVA, infection, CT) ---
step("03 Analysis (incl. nested LOOCV for Fig 2l)", "scripts/03_Analysis.R")

# --- Step 4: MOFA2 (full + RA-only) ---
step("05 MOFA2 RA only",     "scripts/05_MOFA2_RA_only.R")

# --- Step 5: primary figure panels ---
step("06 Figures (all main + supplementary panels and Supplementary Data 3)", "scripts/06_Comprehensive_Figures.R")


cat("\nAll steps completed.\n")
cat("Outputs:\n")
cat("  results/RA_ILD_Workspace.RData\n")
cat("  results/MOFA2_*_Results.RData\n")
cat("  results/tables/*.csv\n")
cat("  results/panels/*.png, *.pdf\n")
