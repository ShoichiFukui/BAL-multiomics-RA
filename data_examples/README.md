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
