# Bronchoalveolar lavage multi-omics identifies immune signatures associated with infection and lung disease progression in rheumatoid arthritis

Analysis code accompanying the article by **Fukui, Hosogaya et al.** (Arthritis & Rheumatology).

## Study summary

Prospective observational cohort study of bronchoalveolar lavage (BAL) fluid from
**24 patients with rheumatoid arthritis (RA) and 7 with sarcoidosis** (n = 31),
integrating six data modalities: BAL RNA-seq, BAL flow cytometry (FCM), peripheral
blood (PB) FCM, BAL cytokines, serum cytokines (41-plex) and the BAL 16S microbiome.
The study addresses:

1. the features distinguishing the RA lung from sarcoidosis (descriptive comparison),
2. subsequent respiratory infection within RA (6 of 24 patients; median follow-up
   5.1 years; Cox model with profile-likelihood CI),
3. decline in pulmonary function, the main progression measure (change in %FVC
   between the test nearest to BAL and the test nearest to the follow-up CT,
   n = 22). The CT measure (change in total opacity = GMM ground-glass opacity +
   consolidation + dense tissue) is exploratory.

Key findings (Results of the Abstract): Compared with sarcoidosis, the RA lung was
distinguished mainly by the BAL CD4+ T cell composition, with a smaller Th17.1
(CCR6+CXCR3+) proportion. Patients who later developed respiratory infection had a
lower BAL Th17.1 proportion, and lower expression of oxidative phosphorylation genes
and of the surfactant-related genes *SFTPB* and *ETV5* (BAL Th17.1 area under the
curve 0.87, 95% confidence interval 0.59–1.00). A lower BAL Th17.1 count was
associated with a greater decline in %FVC (ρ = 0.60), whereas blood Th17.1 cells
were not. No baseline marker was associated with the CT-quantified change. Neither
BAL Th17.1 association survived false discovery rate correction.

## Repository contents

```
BAL-multiomics-RA/
├── README.md                     # this file
├── LICENSE                       # MIT License (code)
├── CITATION.cff                  # citation metadata
├── CHANGELOG.md                  # version history
├── .gitignore
├── requirements.txt              # Python dependencies (CT GMM and MOFA2 environments)
├── R_packages.txt                # R package list
├── scripts/                      # analysis code of the article (n = 31)
│   ├── README.md                 # execution order, inputs, article item → script map
│   ├── 00_cohort_config.R … 18_SD4_microbiome_genus_exact.R   # figures and tables
│   ├── figures/                  # figure panels
│   └── CT_GMM/                   # CT GMM tissue classification (production script + prototype notebook)
├── results/                      # genome-wide DESeq2 and GSEA tables (gene / term level, no patient data)
├── docs/
│   ├── data_availability.md      # repositories and access to the individual-level data
│   └── package_versions.txt      # R/Python package versions of the environment used to check the scripts
└── data_examples/                # README only: where the inputs go (no patient-level data)
```

## Analysis scripts

The scripts in `scripts/` produce the figures and tables of the article on the 31
analysed patients, who are defined in the controlled-access table
`cohort_definition.csv`. They start from the processed individual-level data
(per-patient table, cytokine and flow cytometry matrices, bulk RNA-seq counts, 16S
OTU table), which are available on request; `01a_reference_atlas.R` builds the
scRNA-seq reference from public GEO data. [`scripts/README.md`](scripts/README.md) gives the execution
order, the input manifest and a **table mapping every item of the article
(Figures 1–5, Supplementary Figures 1–4, Supplementary Tables 1–8, Supplementary
Data 1–9) to the script and output file that produce it**.

`results/` holds the genome-wide DESeq2 tables (RA vs sarcoidosis, n = 31;
infection within RA, n = 24) and the four GSEA tables. Shrunken log2 fold changes
(`DESeq2::lfcShrink`, `type = "normal"`) are computed on the Cook's-outlier-replaced
counts so that their sign agrees with the Wald statistic, and P and adjusted P are
those of the default `DESeq()` fit. See [`results/README.md`](results/README.md).

## Data availability

Individual-level data are **not included**. RNA-seq: GEO GSE329884; 16S: SRA PRJNA1462027; per-figure
values: Source Data, Supplementary Tables and Supplementary Data of the article; genome-wide DESeq2 and GSEA
tables: `results/`. Individual-level clinical, cytokine, flow cytometry, CT, pulmonary-function and outcome
data are available from the corresponding author on reasonable request under a data sharing agreement (no
costs). The files the scripts expect are listed in [`data_examples/README.md`](data_examples/README.md); see
[`docs/data_availability.md`](docs/data_availability.md) for details.

## Software requirements

### R

Tested with **R 4.5.0** (macOS). Major packages: `DESeq2`, `edgeR`, `Seurat` (v5),
`harmony`, `BayesPrism`, `clusterProfiler`, `fgsea`, `GSVA`, `msigdbr`, `MOFA2`, `pROC`, `survival`,
`logistf`, `ppcor`, `DescTools`, `effsize`, `vegan`, `tidyverse`, `patchwork`,
`readxl`. Full list: [`R_packages.txt`](R_packages.txt); installer:
`scripts/install_dependencies.R`.

### Python

Two separate environments (see [`requirements.txt`](requirements.txt)):

- **CT GMM** (`scripts/CT_GMM/lung_gmm_batch.py`): Python 3.12 with scikit-learn 1.8.0,
  numpy 1.26.4, SimpleITK, pydicom, lungmask 0.2.20 and PyTorch 2.2.2 (the environment in which the
  article's CT measures were computed; exact versions in `requirements.txt`). The pipeline writes one row per
  scan (`<case key>_T<k>`: tissue percentages and volumes) to `batch_summary.csv`. The per-patient table
  `planA_data_RA24.csv` holds these values for the baseline and follow-up scans (`CT_T1_*`, `CT_T2_*`,
  `CT_Delta_* = T2 − T1`); aerated lung `T1_aer`/`T2_aer` = Air + Healthy Lung (%), total opacity =
  100 − aerated lung, and `d_aer` = `T2_aer` − `T1_aer`.
- **MOFA2** backend: `mofapy2`, used from R through `reticulate`; the interpreter is
  set by the environment variable `MOFA_PYTHON`.

`scripts/CT_GMM/lung_gmm_interactive.ipynb` is an optional 3D Slicer prototype; the
article used `lung_gmm_batch.py`.

## Quick reproduction

```bash
python  scripts/CT_GMM/lung_gmm_batch.py        # CT tissue classification (DICOM + lungmask)
Rscript scripts/01a_reference_atlas.R           # scRNA-seq reference from public GEO data (~32 GB RAM)
Rscript scripts/01b_BayesPrism_n31.R            # deconvolution of the 31 samples
Rscript scripts/02_transcriptome_RA_vs_sarcoidosis.R
# ... remaining scripts in the order given in scripts/README.md
```

## Cohort

- **RA** (n = 24): 2010 ACR/EULAR classification criteria; 11 with and 13 without ILD at baseline.
- **Sarcoidosis** (n = 7).
- The analysed patients are defined in the controlled-access table
  `cohort_definition.csv`, read by `scripts/00_cohort_config.R`; no
  patient identifier appears in the code.
- **Healthy serum reference** (n = 101): Saza cohort; 42 men, 59 women; mean age 58 (SD 9.7).

## Statistical conventions

- Wilcoxon rank-sum and Spearman tests are called with `exact = TRUE`; R uses the
  exact distribution unless there are ties, in which case it uses the normal
  approximation (Methods: "exact unless tied"). All tests are two-sided.
- Random seeds, as run: `set.seed(42)` before each bootstrap CI, permutation
  test, GO GSEA and MOFA2 training (`seed = 42`). The KEGG GSEA in 02 and 03
  follows the GO GSEA after a single `set.seed(42)`.
- Multiple testing: **Benjamini–Hochberg (FDR) only**, within each analysis family
  (DEG, GSEA, microbiome, the RA-vs-sarcoidosis screen of
  Supplementary Data 2, and the biomarker screens). The biomarker screens use **one
  common panel of 114 directly measured baseline markers** for every outcome (infection, change in
  %FVC, change in total opacity), adjusted within each outcome. The family is stated
  in each Supplementary Table/Data legend. Secondary small-sample analyses
  (Supplementary Tables 5–7) are reported with nominal P values.
- Gene-set scores: per-sample GSVA scores (kcdf Gaussian, log-CPM of the 31
  samples) of the GO Biological Process sets GOBP_OXIDATIVE_PHOSPHORYLATION and
  GOBP_SURFACTANT_HOMEOSTASIS (msigdbr 25.1.1); they are derived variables
  outside the FDR family, like the BayesPrism estimates (`17_GSVA_GO_scores.R`).
- DESeq2 log2 fold changes are shrunk with `type = "normal"`.

## License

Code: MIT License (see [`LICENSE`](LICENSE)). The license covers the code only.
Patient-level data and protected health information are not covered by it; they
are governed by the terms of the repositories (GEO, SRA) and, for individual-level
data, by the data sharing agreement.

## Citation

If you use this code in academic work, please cite the article and this software.
See [`CITATION.cff`](CITATION.cff).

Archived on Zenodo: [10.5281/zenodo.20699675](https://doi.org/10.5281/zenodo.20699675)
(concept DOI; resolves to the latest version). Source:
https://github.com/ShoichiFukui/BAL-multiomics-RA

## Contact

- Corresponding author: Naoki Hosogaya, M.D., Ph.D., Clinical Research Center, Nagasaki University Hospital
- Ethics approval: Nagasaki University Hospital Institutional Review Board, approval numbers 19021801 and 20051819
