# ============================================================================
# 00_run_all.R
# End-to-end runner for the manuscript analysis pipeline.
# Sources each step in order. Heavy step 01 (BayesPrism) can be skipped if its
# outputs already exist.
#
# Usage:  Rscript 00_run_all.R
# ============================================================================

# Resolve the script directory so this runner works from any CWD.
this_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(this_dir)

step <- function(label, path) {
  cat(sprintf("\n========== %s ==========\n", label))
  source(path, local = FALSE, echo = FALSE)
}

# --- Step 1: scRNA-seq atlas + BayesPrism (heavy; skip if already built) ---
if (!file.exists("../output_v3_BayesPrism/celltype_proportions_cellFraction.csv")) {
  step("01 BayesPrism deconvolution", "01_BayesPrism_Deconvolution.R")
} else {
  cat("\nSkipping 01 (BayesPrism outputs already exist)\n")
}

# --- Step 2: master workspace ---
step("02 PostDeconvolution", "02_PostDeconvolution.R")

# --- Step 3: primary analyses (DEG, GSEA, GSVA, infection, CT) ---
step("03 Analysis", "03_Analysis.R")
step("03b Nested LOOCV RA vs sarcoidosis (Fig 2l)", "03b_NestedLOOCV_RAvsSarc.R")

# --- Step 4: MOFA2 (full + RA-only) ---
step("04 MOFA2 full cohort", "04_MOFA2_FullCohort.R")
step("05 MOFA2 RA only",     "05_MOFA2_RA_only.R")

# --- Step 5: primary figure panels ---
step("06 Comprehensive figures", "06_Comprehensive_Figures.R")

# --- Step 6: final-version figure updates ---
update_scripts <- sort(list.files("figure_updates", pattern = "\\.R$", full.names = TRUE))
for (s in update_scripts) {
  step(sprintf("Update: %s", basename(s)), s)
}

cat("\nAll steps completed.\n")
cat("Outputs:\n")
cat("  results/RA_ILD_Workspace.RData\n")
cat("  results/MOFA2_*_Results.RData\n")
cat("  results/tables/*.csv\n")
cat("  results/panels/*.png, *.pdf\n")
