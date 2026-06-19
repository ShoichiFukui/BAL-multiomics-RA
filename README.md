# Multi-omics characterization of bronchoalveolar lavage fluid in rheumatoid arthritis

Analysis code accompanying the manuscript by **Fukui, Hosogaya et al.**
(*Nature Communications*, under review).

## Study summary

Prospective observational cohort study integrating six independent BAL data
modalities (host RNA-seq, BAL cytokines, serum cytokines, BAL flow cytometry,
peripheral blood flow cytometry, and BAL 16S microbiome) from **24 rheumatoid
arthritis (RA) patients and 11 sarcoidosis controls** (n = 35 total) to:

1. identify molecular features distinguishing the RA lung mucosal compartment
   from sarcoidosis,
2. prospectively predict future respiratory infection within RA, and
3. correlate baseline BAL parameters with longitudinal CT-quantified lung
   disease progression.

Key findings: a compartmentalized Th17-lineage immune axis (local BAL Th17.1
reduction predicting infection; systemic PB Th17 expansion linked to ILD
progression), in addition to RA-specific B cell / plasma cell activation in the
lung mucosa.

## Repository contents

```
github_release/
├── README.md                     # this file
├── LICENSE                       # MIT License for code
├── CITATION.cff                  # citation metadata
├── .gitignore
├── requirements.txt              # Python dependencies
├── R_packages.txt                # R package list
├── scripts/                      # full reproducible analysis pipeline
│   ├── README.md                 # pipeline-level documentation
│   ├── 01_BayesPrism_Deconvolution.R
│   ├── 02_PostDeconvolution.R
│   ├── 03_Analysis.R
│   ├── 04_MOFA2_FullCohort.R
│   ├── 05_MOFA2_RA_only.R
│   ├── 06_Comprehensive_Figures.R
│   ├── analysis_modules/         # sourced by 03_Analysis.R
│   ├── figure_updates/           # final-version panel updates
│   └── CT_GMM/                   # CT GMM tissue classification
│       ├── lung_gmm_batch_v3_FINAL.py    # production GMM pipeline
│       └── lung_gmm_interactive.ipynb    # interactive Slicer notebook (prototype)
├── docs/
│   ├── data_availability.md      # public data repository links
│   └── pipeline_overview.md      # step-by-step pipeline guide
└── data_examples/                # template/dummy CSVs (no patient-level data)
```

## Data availability

Patient-level raw and processed data are **not included** in this repository to
protect participant privacy. Public sequence repositories provide access under
the schedules below.

| Modality | Repository | Accession | Status |
|---|---|---|---|
| BAL RNA-seq (FASTQ + counts) | GEO | GSE329884 | Approved; **public release 2027-05-02** |
| 16S rRNA microbiome (FASTQ) | NCBI BioProject / SRA | PRJNA1462027 | Submitted; release per repository policy |
| Healthy serum cytokine reference (n = 101, Saza cohort) | Reported in Ref. 21 (Saza cohort) | See manuscript Methods § Cytokine profiling | Restricted; available on request |
| Public scRNA-seq BAL reference atlas | GEO | GSE145926, GSE193782, GSE184735 | Public |

The code runs on the **processed, de-identified n = 35 dataset** (the three
QC-excluded samples are already removed); the raw 37-sample data is not
distributed. Only the minimal processed data needed to reproduce the reported
results is made available, in two tiers:

- **Open:** sequence data (GEO GSE329884; SRA PRJNA1462027) and this code.
- **Controlled access on request:** individual-level clinical metadata, BAL/serum
  cytokine concentrations, BAL/PB flow cytometry, CT-derived tissue fractions, and
  infection outcomes — from the corresponding author to qualified investigators
  under an IRB-approved data sharing agreement (Nagasaki University Hospital IRB
  19021801 and 2005819).

This managed-access arrangement for sensitive patient data is consistent with
*Nature Communications* policy: the data-availability statement specifies the
access conditions, and the underlying data are provided to the editors and
referees during peer review on request. A per-figure **Source Data** file
accompanies the manuscript.

See [`docs/data_availability.md`](docs/data_availability.md) for full details
and access URLs.

## Software requirements

### R (≥ 4.3)

Major packages: `DESeq2`, `Seurat (≥ 4)`, `Harmony`, `BayesPrism`, `fgsea`,
`GSVA`, `MOFA2`, `pROC`, `randomForest`, `tidyverse`, `vegan`, `patchwork`,
`openxlsx`. See [`R_packages.txt`](R_packages.txt) for the full list.

### Python (≥ 3.10)

Required for the CT GMM tissue-classification pipeline. Major packages:
`numpy`, `scikit-learn`, `SimpleITK`, `pydicom`, `lungmask`. See
[`requirements.txt`](requirements.txt).

### Optional: 3D Slicer (for the prototype notebook)

`CT_GMM/lung_gmm_interactive.ipynb` is the Slicer-based prototype that informed
the development of the production pipeline. The **production analysis uses
`lung_gmm_batch_v3_FINAL.py`** (pure Python).

## Quick reproduction

```bash
# 1. Heavy: build scRNA-seq reference and BayesPrism deconvolution (~32 GB RAM)
Rscript scripts/01_BayesPrism_Deconvolution.R

# 2. Assemble the master workspace
Rscript scripts/02_PostDeconvolution.R

# 3. Run primary analyses (DEG, GSEA, GSVA, infection, CT progression)
Rscript scripts/03_Analysis.R

# 4. MOFA2 integration
Rscript scripts/04_MOFA2_FullCohort.R
Rscript scripts/05_MOFA2_RA_only.R

# 5. Generate all primary figure panels
Rscript scripts/06_Comprehensive_Figures.R

# 6. Apply final-version panel updates
for f in scripts/figure_updates/*.R; do Rscript "$f"; done

# 7. CT tissue classification (requires DICOM scans + lungmask AI model)
python scripts/CT_GMM/lung_gmm_batch_v3_FINAL.py
```

## Cohort definitions

- **RA** (n = 24): meets 2010 ACR/EULAR criteria.
- **Sarcoidosis controls** (n = 11, KYC prefix): granulomatous lung disease comparator.
- **Excluded samples**: three samples failed QC (one diagnostic reclassification; two with incomplete multi-omics data) and are absent from the released n = 35 dataset.
- **Healthy serum reference** (n = 101): Saza cohort; 42 men, 59 women; mean age 58 (SD 9.7).

## Statistical conventions

- Wilcoxon and Spearman tests use `exact = TRUE`.
- `set.seed(42)` immediately before any stochastic step.
- Multiple testing correction: **Benjamini–Hochberg** within each
  high-dimensional omics analysis (DEG, GSEA, GSVA, microbiome differential
  abundance). **Bonferroni** correction is applied to the targeted 41-cytokine
  multi-comparison screening for CT progression.

## License

Code: MIT License (see [`LICENSE`](LICENSE)).
Data: governed by repository-specific terms (GEO, NCBI), patient-level data
restricted by IRB protocol.

## Citation

If you use this code, please cite the manuscript (citation will be added upon
publication). See [`CITATION.cff`](CITATION.cff) for machine-readable metadata.

A permanent archive of the code with a citable DOI is hosted on Zenodo:
[10.5281/zenodo.20699675](https://doi.org/10.5281/zenodo.20699675) (concept DOI;
always resolves to the latest version).

## Contact

- Corresponding author: Naoki Hosogaya, M.D., Ph.D. (Nagasaki University Hospital)
- Repository maintainer: [contact via corresponding author]
- IRB approval: Nagasaki University Hospital IRB, approval numbers 19021801 and 2005819
