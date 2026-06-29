# Data examples

This directory is intentionally **empty** in the public repository. It exists
to document where each input data type should be placed for the scripts to
run.

```
data_examples/
├── master_data_template.csv        # column headers expected by the scripts
├── master_data_with_CT_all_template.csv  # CT columns + clinical
└── (place your processed inputs here)
```

## Recommended layout (relative to the project root)

```
<project_root>/
├── data/raw/                       # raw inputs (Luminex, FCM, clinical)
│   └── 41plex_healthy_control.xlsx
├── data/processed/                 # cleaned tables (master_data, etc.)
├── output_v3_BayesPrism/           # BayesPrism intermediate + final outputs
│   ├── BAL_reference_author_annotated.rds
│   └── celltype_proportions_cellFraction.csv
├── results/
│   ├── RA_ILD_Workspace.RData      # master workspace built by 02_*
│   └── NatComm/
│       ├── tables/                 # CSV result tables
│       └── panels/                 # PNG/PDF figure panels
└── CT_GMM_score/                   # per-patient GMM JSON outputs
    └── <patient_id>_<study_date>/
        └── gmm_results.json
```

## Reproducing the figures from the provided processed-data package

Because the upstream steps require a single-cell reference and a managed Python
environment (BayesPrism deconvolution in `01_*`; MOFA2 in `03_*`), the
controlled-access data package provided on request (see
[`../docs/data_availability.md`](../docs/data_availability.md)) includes the
**processed intermediates** so that the figures can be regenerated **without
re-running BayesPrism or MOFA2**. Place the package contents at the project root:

| Provided file (de-identified, IDs `RA#` / `Sarcoidosis#`) | Used by | Path |
|---|---|---|
| `clinical_metadata.xlsx`, `cytokine_multiplex.xlsx`, `FCM_integrated_data_transformed.xlsx`, `RA_CT_transformed.xlsx`, `ct_lung_analysis_merged.xlsx` (PHI removed), `sample_id_mapping.xlsx`, `RNAseq/RawCount_heatmap_table.csv`, `OTU…csv`, `data/raw/41plex_healthy_control.xlsx` | `02_*`, `03_*` | project root |
| `output_v3_BayesPrism/celltype_proportions_*.csv` and `BAL_reference_author_annotated.rds` | deconvolution + Fig 1 UMAPs | `output_v3_BayesPrism/` |
| `results/MOFA2_6views_Results.RData`, `results/MOFA2_RA_only_Results.RData` | MOFA factors (Fig 2l, Fig 5) | `results/` |
| `results/CT_Multiomics_Results.RData`, `results/RA_ILD_Workspace.RData` | CT + master data | `results/` |
| `results/tables/*.csv` (processed result tables) and `results/tables/Supplementary_Tables.xlsx` | figure panels + supplementary tables | `results/tables/` |

Then regenerate the figure panels:

```r
Rscript scripts/06_Comprehensive_Figures.R          # main + supplementary panels
# Fig 1–5 assembled PDFs are produced by analysis_modules/Integration_Outputs.R,
# sourced from 03_Analysis.R, or run it directly after the workspace is in place.
```

`master_data_with_CT_all.csv` (consumed by `06_*`) is written automatically by
`analysis_modules/CT_Multiomics.R`; it is also included in the package.
The analysis numbers themselves (e.g. DEG count, classifier AUCs) are
regenerated end-to-end by `02_*`/`03_*` from the de-identified inputs above.

## How to obtain the underlying data

See [`../docs/data_availability.md`](../docs/data_availability.md) for full
access information (GEO / SRA accessions for sequence data; data sharing
agreement procedure for clinical and cytokine data).

## Why are these files not committed?

Per the Nagasaki University Hospital IRB protocol (approval numbers 19021801
and 2005819), patient-level data may not be redistributed outside of approved
data sharing agreements. Aggregating templates and column-header descriptions
in this directory keeps the analysis pipeline reproducible without disclosing
restricted information.
