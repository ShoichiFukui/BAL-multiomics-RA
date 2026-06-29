# Changelog

## v1.1.0 (2026-06-30)
Synchronized with the Nature Communications submission; the released pipeline now
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
