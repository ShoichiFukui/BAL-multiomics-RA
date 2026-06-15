# Pipeline overview

This document summarizes the end-to-end analysis pipeline, the data flow
between scripts, and the figures / tables that each step produces.

```
┌──────────────────────────────────────────────────────────────────────┐
│  Raw inputs (restricted)                                             │
│  ─────────────────────                                               │
│  • BAL RNA-seq FASTQ          ──► GEO GSE329884                      │
│  • BAL 16S rRNA FASTQ         ──► SRA PRJNA1462027                   │
│  • BAL / serum cytokines      ──► Luminex 41-plex                    │
│  • BAL / PB flow cytometry    ──► FlowJo gated proportions           │
│  • DICOM chest CT (×2 dates)  ──► CT_GMM_score/*/gmm_results.json    │
│  • Clinical metadata          ──► RA_ILD_Workspace.RData             │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 1. scRNA-seq atlas + BayesPrism deconvolution                  │
│  ─────────────────────────────────────────────────                   │
│  scripts/01_BayesPrism_Deconvolution.R                               │
│    • Loads three public BAL scRNA-seq datasets                       │
│      (GSE145926 COVID-19, GSE193782 Healthy, GSE184735 Sarcoidosis)  │
│    • Integrates with Harmony (Korsunsky 2019)                        │
│    • k-NN transfers cell-type labels for GSE184735                   │
│    • Runs BayesPrism (Chu 2022) on the bulk BAL RNA-seq              │
│  Outputs:                                                            │
│    output_v3_BayesPrism/celltype_proportions_cellFraction.csv        │
│    output_v3_BayesPrism/BAL_reference_author_annotated.rds           │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 2. Build the master workspace                                  │
│  ─────────────────────────────────                                   │
│  scripts/02_PostDeconvolution.R                                      │
│    • Merges deconvolution + cytokines + FCM + clinical               │
│  Outputs:                                                            │
│    results/RA_ILD_Workspace.RData (master_data, EXCL, etc.)          │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 3. Primary analyses                                            │
│  ─────────────────────                                               │
│  scripts/03_Analysis.R   (master, sources analysis_modules/)         │
│    ├─ Enhanced_Analysis.R        ── DEG, GSEA, GSVA, effect sizes    │
│    ├─ Infection_Prediction.R     ── BAL Th17.1 / OXPHOS              │
│    ├─ CT_Multiomics.R            ── CT progression analyses          │
│    ├─ Integration.R              ── cross-layer correlations          │
│    ├─ MOFA2.R                    ── MOFA2 wrapper                     │
│    └─ Final7.R                   ── final integration outputs        │
│  Outputs:                                                            │
│    results/tables/*.csv                                       │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 4. Multi-omics factor analysis                                 │
│  ───────────────────────────────────                                 │
│  scripts/04_MOFA2_FullCohort.R      ── full-cohort, n=35, 6 views    │
│  scripts/05_MOFA2_RA_only.R         ── RA-only, n=24, 5 views        │
│  Outputs:                                                            │
│    results/MOFA2_6views_Results.RData                        │
│    results/MOFA2_RA_only_Results.RData                       │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 5. Primary figure panels                                       │
│  ─────────────────────────────                                       │
│  scripts/06_Comprehensive_Figures.R                            │
│  Outputs:                                                            │
│    results/panels/Fig{1..5}*.png|.pdf                        │
│    results/panels/FigS{1..4}*.png|.pdf                       │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Step 6. Final panel updates                                         │
│  ──────────────────────────                                          │
│  scripts/figure_updates/                                             │
│    • Fig 5 CUD-palette ROCs and ILD-coloured scatter                 │
│    • SFig 2 (alpha-consistent serum), SFig 3 (ILD stratification),   │
│      SFig 4 (RA-specific medication confounding)                     │
│    • Supplementary Tables 1–9 (rebuilt with current CSVs)            │
└──────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌──────────────────────────────────────────────────────────────────────┐
│  Side track. CT GMM tissue classification                            │
│  ───────────────────────────────────────                             │
│  scripts/CT_GMM/lung_gmm_batch_v3_FINAL.py                           │
│    • Auto-selects thin-slice non-contrast lung-kernel series         │
│    • Lungmask (U-net R231; Hofmanninger 2020) for lung mask          │
│    • scikit-learn GaussianMixture (5 components, full covariance)    │
│      following Zaffino et al. 2021                                   │
│  Outputs:                                                            │
│    CT_GMM_score/<patient>_<date>/gmm_results.json                    │
│      → consumed by 02_PostDeconvolution.R as CT_* columns            │
└──────────────────────────────────────────────────────────────────────┘
```

## Reproducibility notes

- **Random seeds.** `set.seed(42)` is set immediately before each stochastic
  step (random forest, bootstrap CI, permutation testing). Python uses
  `random_state=42` for the GMM.
- **Cohort.** n = 35 (RA = 24, sarcoidosis = 11). Three additional samples failed
  QC (one diagnostic reclassification; two with incomplete multi-omics data) and
  are absent from the released dataset.
- **Statistical tests.** Wilcoxon and Spearman tests use `exact = TRUE`.
- **Multiple testing correction.** Benjamini–Hochberg within each
  high-dimensional omics analysis; Bonferroni for the targeted 41-cytokine
  multi-comparison screening of CT progression.
- **Software versions.** R ≥ 4.3, Python ≥ 3.10, 3D Slicer ≥ 5 (optional).
  Critical pin: `scikit-learn==0.23.2` for the GMM, matching the version used
  by Zaffino et al. 2021.
- **Online annotation databases (KEGG).** KEGG GSEA fetches pathway annotations
  live from the KEGG REST API at run time, so the exact set of enriched
  pathways can drift slightly as KEGG is updated (the manuscript reported 81
  KEGG pathways at the analysis date; a later run may differ by a few). GO terms
  use a fixed local MSigDB/org.Hs.eg.db and are stable. For exact reproduction,
  record the KEGG release date or cache the KEGG annotation used.
- **`predict()` namespace.** Random-forest probability predictions are called as
  `stats::predict(rf, …, type = "prob")` so the `randomForest` S3 method is
  dispatched even when packages that define their own `predict` (e.g. MOFA2,
  Seurat) are attached. This is behaviour-preserving — it only avoids a masking
  error under newer R/package versions.

## File-to-figure mapping

| Figure | Generating script |
|---|---|
| Fig 1 | `06_Comprehensive_Figures.R` |
| Fig 2 | `06_Comprehensive_Figures.R` |
| Fig 3 | `06_Comprehensive_Figures.R` |
| Fig 4 | `06_Comprehensive_Figures.R` |
| Fig 5 | `06_Comprehensive_Figures.R` + several `figure_updates/*` |
| Suppl Fig 1 | `06_Comprehensive_Figures.R` |
| Suppl Fig 2 | `figure_updates/20260604_regen_FigS3_alpha_consistent.R` |
| Suppl Fig 3 | `figure_updates/20260604_make_SFig2_infection_ILD_stratified.R` |
| Suppl Fig 4 | `figure_updates/20260605_medication_confound_RAfactor.R` |
| Suppl Tables 1–9 | `figure_updates/20260605_rebuild_supp_tables_renumbered.R` |
