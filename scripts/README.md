# Hosogaya BAL Multi-omics — Final Analysis Scripts

This directory contains the curated, reproducible R / Python pipeline used in the
Nature Communications manuscript on rheumatoid arthritis (RA) BAL multi-omics
profiling versus sarcoidosis controls (n = 35; RA = 24, sarcoidosis = 11).

## Directory layout

```
scripts/
├── README.md                                # this file
├── 00_run_all.R                             # end-to-end runner (sources steps in order)
├── 01_BayesPrism_Deconvolution.R            # scRNA-seq atlas + BayesPrism deconvolution
├── 02_PostDeconvolution.R                   # builds RA_ILD_Workspace.RData
├── 03_Analysis.R                            # master analysis script (sources analysis_modules/)
├── 04_MOFA2_FullCohort.R                    # full-cohort MOFA2 (n = 35, 6 views)
├── 05_MOFA2_RA_only.R                       # RA-specific MOFA2 (n = 24, 5 views)
├── 06_Comprehensive_Figures.R               # generates all primary figure panels
├── analysis_modules/                        # sourced by 03_Analysis.R
│   ├── Enhanced_Analysis.R                  # DEG, GSEA, GSVA, effect sizes
│   ├── Infection_Prediction.R              # BAL Th17.1 / OXPHOS infection prediction
│   ├── CT_Multiomics.R                     # CT-quantified ILD progression analysis
│   ├── Integration.R                       # cross-layer correlation analysis
│   └── Final7.R                            # final integration outputs
│   # (the 6-view MOFA2 model is built directly in 03_Analysis.R)
├── figure_updates/                          # final panel revisions
│   ├── 20260601_make_Fig5_infection_ROC_RAfactor.R
│   ├── 20260604_regen_Fig5_panels_descriptive.R
│   ├── 20260605_Fig5c_remove_loocv_annot.R     # Fig 5c without nested-LOOCV annotation
│   ├── 20260605_recolor_Fig5_ROC_CUD.R         # CUD palette for Fig 5c/5f/5i
│   ├── 20260605_regen_Fig5h_ILDcolor_CUD.R     # Fig 5h ILD status coloured (CUD)
│   ├── 20260604_make_SFig2_infection_ILD_stratified.R  # Suppl Fig 3 panel (b)
│   ├── 20260604_regen_FigS3_alpha_consistent.R         # Suppl Fig 2 (serum, α-titles)
│   ├── 20260605_regen_FigS4_factorlabels.R             # Suppl Fig 4 axis labels
│   ├── 20260605_medication_confound_RAfactor.R         # Suppl Fig 4 + ST9 (RA Factor 1 / BAL Th17.1)
│   └── 20260605_rebuild_supp_tables_renumbered.R       # Suppl Tables 1–9 rebuild
└── CT_GMM/
    ├── lung_gmm_batch_v3_FINAL.py                       # production batch GMM pipeline (used for the manuscript)
    └── lung_gmm_interactive.ipynb                       # 3D Slicer Jupyter notebook (optional prototype)
```

## Required software

- **R ≥ 4.3** (DESeq2, Seurat ≥ 4, Harmony, BayesPrism, fgsea, GSVA, MOFA2,
  pROC, randomForest, tidyverse, patchwork, openxlsx)
- **Python ≥ 3.10** (numpy, scikit-learn, SimpleITK, pydicom, lungmask) — for the
  CT GMM pipeline (`CT_GMM/lung_gmm_batch_v3_FINAL.py`); also python-pptx /
  openpyxl for figure packaging utilities
- **3D Slicer v5 + DensityLungCT GMM Segmentation extension** — optional, only
  required if you want to re-run the interactive prototype notebook
  `CT_GMM/lung_gmm_interactive.ipynb`. The manuscript analyses used
  `lung_gmm_batch_v3_FINAL.py` (pure Python), not Slicer.
- Sufficient memory (~32 GB) for `01_BayesPrism_Deconvolution.R`

## Execution order

```bash
# 1. Reference atlas + deconvolution (very heavy; one-time)
Rscript 01_BayesPrism_Deconvolution.R

# 2. Build the master workspace
Rscript 02_PostDeconvolution.R

# 3. Main analyses (DEG, GSEA, GSVA, effect sizes, infection, CT progression)
Rscript 03_Analysis.R

# 4. MOFA2 (full-cohort then RA-only)
Rscript 04_MOFA2_FullCohort.R
Rscript 05_MOFA2_RA_only.R

# 5. Generate all primary panels (Fig 1–5, Suppl Fig 1–4)
Rscript 06_Comprehensive_Figures.R

# 6. Apply final-version panel updates (Fig 5 CUD palette, SFig labels, etc.)
for f in figure_updates/*.R; do Rscript "$f"; done

# 7. CT-quantified ILD progression (production pipeline used for the manuscript)
python CT_GMM/lung_gmm_batch_v3_FINAL.py
#    Or, optionally, run the prototype notebook in 3D Slicer:
#    CT_GMM/lung_gmm_interactive.ipynb
```

## Manuscript-cited statistical conventions

- All Wilcoxon and Spearman tests use `exact = TRUE`.
- `set.seed(42)` is set immediately before any stochastic step
  (random forest, bootstrap CI, permutation testing).
- Multiple testing correction: Benjamini-Hochberg within each high-dimensional
  omics analysis (DEG, GSEA, GSVA, microbiome). Bonferroni correction is
  applied to the targeted 41-cytokine multi-comparison screening for CT
  progression.
- Cohort: n = 35 (RA = 24, sarcoidosis = 11). Three additional samples failed QC
  (one diagnostic reclassification; two with incomplete multi-omics data) and are
  absent from the released dataset.

## Cohort definitions

- **RA**: 2010 ACR/EULAR criteria, n = 24
- **Sarcoidosis controls (KYC prefix)**: granulomatous lung disease comparator, n = 11
- **Healthy serum reference (Saza cohort)**: n = 101; 42 men, 59 women; mean age 58 (SD 9.7)

## Output paths (relative to the parent project directory)

- `results/RA_ILD_Workspace.RData` — master workspace
- `results/MOFA2_6views_Results.RData` — full-cohort MOFA2
- `results/MOFA2_RA_only_Results.RData` — RA-only MOFA2
- `results/tables/*.csv` — all numeric results tables
- `results/panels/*.png|.pdf` — all figure panels

## Notes

- The R scripts assume the working directory is the parent project directory
  (e.g. `~/Hosogaya_BAL/`).
- Raw FASTQ and 16S sequences are deposited at GEO (GSE329884) and SRA
  (BioProject PRJNA1462027) respectively; release schedule per repository policy.
