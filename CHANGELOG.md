# Changelog

## v1.1.2 (2026-07-13)
Clarified reproducibility of the full-cohort 6-view MOFA2 model.

- `03_Analysis.R`: documented that the MOFA input `results/Integration_7layers_Results.RData`
  (object `layers`) is a patient-level processed intermediate provided on request with the
  controlled-access data package (not committed to this repository), and that its feature
  selection is specified in the Supplementary Methods (Expression, top 20 genes; BAL and
  serum cytokines, top 15 each by log2-scale variance; flow cytometry and genus-level
  microbiome retained after removing zero-variance features; 6-view, 92 features). The load
  now fails with an explicit message pointing to `docs/data_availability.md` instead of being
  reported as a Python/MOFA2 environment error.
- `README.md`, `scripts/README.md`: corrected the statement that the 6-view MOFA2 is "built
  directly in 03_Analysis.R" to reflect that it is *fit* from the provided processed input.
- `data_examples/README.md`: added `Integration_7layers_Results.RData` to the controlled-access
  processed-data manifest. No patient-level data are added to the repository.

## v1.1.1 (2026-07-10)
Clarified multiple-testing (FDR) families and made the release journal-neutral.

- **Multiple-testing families.** The BH-adjusted P values in Supplementary Table 4 are
  corrected only within that 17-variable exploratory comparison table and are NOT
  study-wide FDR values (its multi-omics markers were selected from the top of the full
  screens). Study-wide FDR for the two headline findings is reported in Supplementary
  Data 3 (infection screen, 143 markers; BAL Th17.1 BH P = 0.112) and Supplementary
  Data 7 (CT progression screen, 117 markers; PB Th17 BH P = 0.451). These are the
  values quoted in the manuscript Abstract.
- `06_Comprehensive_Figures.R`: the Supplementary Table 4 sheet now also emits the
  within-table BH-adjusted P columns, so the released code reproduces the submitted table.
- `scripts/README.md`: correction-family documentation expanded accordingly.
- Documentation, comments and output filenames made journal-neutral.

## v1.1.0 (2026-06-30)
Synchronized with the manuscript submission; the released pipeline now
contains only manuscript-reported analyses.

- Removed orphan analyses not in the manuscript/SI: GO-MF GSEA, elastic-net CT
  prediction, random-forest feature importance, leaky semi-nested LOOCV (and its
  leave-one-out AUC sensitivity), Cook's distance, Jonckheere-Terpstra
  dose-response, inert 6th MOFA microbiome view, cross-layer correlation network.
- Fig 2l: fully fold-internal nested LOOCV (leakage-free AUC = 0.780), replacing
  the semi-nested estimate (AUC ~0.962).
- Figure-panel output names aligned to manuscript figure numbers.
- Full-cohort 6-view MOFA2 consolidated into 03_Analysis.R.
- Manuscript title updated to "...immune signatures associated with...".

## v1.0.0
Public release of the analysis code accompanying the manuscript
"Bronchoalveolar lavage multi-omics identifies immune signatures predicting
infection and lung disease progression in rheumatoid arthritis"
(Fukui, Hosogaya et al.).

- Reproducible R/Python pipeline: scRNA-reference BayesPrism deconvolution,
  differential expression, GSEA/GSVA, prospective infection prediction,
  CT GMM tissue quantification, MOFA2 multi-omics integration, and figure generation.
- Patient-level data are not distributed; see docs/data_availability.md for access.
