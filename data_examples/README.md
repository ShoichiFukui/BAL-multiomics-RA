# Data inputs

This directory contains no data. It documents which input files the scripts expect and
where to place them. Individual-level data are available from the corresponding author on
reasonable request to qualified investigators under a data sharing agreement (no costs are
associated with accessing these data); see
[`../docs/data_availability.md`](../docs/data_availability.md).

## De-identification

The individual-level files contain no names, initials, dates of birth or hospital record
numbers. Samples carry study identifiers that are mapped to the de-identified public labels
(`RA#`, `Sarcoidosis#`, as in GEO and the Source Data) only in `cohort_definition.csv`. All
scripts assign groups, exclusions and labels through this table
(`scripts/00_cohort_config.R`: `sample_group()`, `is_sarc()`, `public_label()`) and never
through identifier patterns.

**All files shared on request carry times only as days from the BAL date; no file contains a
calendar date.** This covers the CT scans, pulmonary function tests, infection events and the
last visit. The CT GMM pipeline names its outputs by case key and time point
(`<case key>_T<k>`) and writes no study date, age or DICOM UID.

## Column names

All files shared on request use English column names. The BAL cell counts of
`results/tables/master_data_with_CT_all.csv` are `BAL_total_cells_per_mL`,
`BAL_lymphocytes_per_mL`, `BAL_macrophage_pct`, `BAL_lymphocyte_pct`, `BAL_neutrophil_pct`,
`BAL_Th171_count_sheet` (BAL Th17.1 count per mL of the clinical sheet), `BAL_Treg_count_sheet`
(BAL Treg count per mL) and `BAL_neutrophil_CD66b_CD206neg_pct_sheet` (BAL CD66b+CD206− neutrophils, %, of the clinical sheet). `check_english_columns()` in
`scripts/00_cohort_config.R` stops if a table carries a non-English column name.

## Layout (relative to the project root)

```
<project_root>/
├── results/
│   ├── RA_ILD_Workspace.RData                 # processed data (per-patient table, FCM, cytokine and microbiome matrices,
│   │                                          #   bulk counts, gene annotation)
│   └── tables/master_data_with_CT_all.csv     # per-patient table
├── RNAseq/RawCount_heatmap_table.csv          # bulk RNA-seq counts (also in GEO GSE329884)
├── OTU abundance table-2 (merged abundance table).csv   # 16S OTU table (SRA PRJNA1462027)
├── FCM_integrated_data_transformed.xlsx       # BAL / PB flow cytometry and BAL differential counts
├── output_BayesPrism/BAL_reference_author_annotated.rds  # written by scripts/01a_reference_atlas.R (public data)
├── GSE145926_RAW/, reference/, analysis/      # public scRNA-seq data for 01a (GEO; see below)
└── article/                                   # ART_DIR
    ├── cohort_definition.csv
    ├── 41plex_healthy_control.xlsx
    └── tables/n31/
        ├── planA/planA_data_RA24.csv
        ├── followup_by_sample.csv
        ├── infection_events.csv
        ├── pft_series_by_sample.tsv
        ├── sarc_BAL_WBC.csv
        └── medication_at_BAL.csv
```

## Public data for `scripts/01a_reference_atlas.R`

| File | Content |
|---|---|
| `GSE145926_RAW/*.h5` | 10x count matrices of the COVID-19 / healthy BAL scRNA-seq dataset GSE145926 (Liao et al., Nat Med 2020) |
| `reference/GSE145926_all_cell_annotation.txt`, `reference/GSE145926_myeloid_annotation.txt` | the authors' cell annotations (`all.cell.annotation.meta.txt`, `myeloid.cell.annotation.meta.txt`, https://github.com/zhangzlab/covid_balf) |
| `analysis/GSE193782_SeuratObject.rds` | healthy BAL scRNA-seq Seurat object of GSE193782 (Li et al., Am J Respir Cell Mol Biol 2022) |
| `analysis/GSE184735/extracted/` | 10x matrices of the sarcoidosis BAL scRNA-seq dataset GSE184735 (Liao et al., Eur Respir J 2021), one folder per sample |

## Individual-level files (on request)

| File | Content | Used by |
|---|---|---|
| `results/RA_ILD_Workspace.RData` | processed data of all sequenced samples; the scripts keep the 31 analysed patients | 02–06, 08e, 15–17, figures |
| `results/tables/master_data_with_CT_all.csv` | one row per patient: clinical variables, BAL cell counts, cytokines, flow cytometry, CT GMM fractions | 05–10, 12, 14, 17, figures |
| `FCM_integrated_data_transformed.xlsx`, sheet `BALF_analysis` | flow cytometry (gated proportions) and BAL differential counts | 11, figures |
| `article/cohort_definition.csv` | `Sample_ID`, `group` (RA / Sarcoidosis), `status` (`analysed` for the 31 patients; `excluded_dx` and `not_in_cohort` mark sequenced samples that are not analysed), `public_id` | all scripts |
| `article/tables/n31/planA/planA_data_RA24.csv` | one row per RA patient (n = 24): baseline clinical variables, RF and anti-CCP positivity, infection outcome, `days_to_follow_up`, pulmonary function at the test nearest to BAL and to the follow-up CT, CT tissue fractions | most scripts |
| `article/tables/n31/followup_by_sample.csv` | days from BAL to the last visit and follow-up status (two censuses) | 08b |
| `article/tables/n31/infection_events.csv` | days from BAL to each of the six infection events | 08b |
| `article/tables/n31/pft_series_by_sample.tsv` | all pulmonary-function tests (FVC, %FVC, FEV1, DLCO, %DLCO) with days from BAL | 10 |
| `article/tables/n31/sarc_BAL_WBC.csv` | BAL total cell count of the sarcoidosis patients | 10 |
| `article/tables/n31/medication_at_BAL.csv` | medication in use at BAL, one row per RA patient, 0/1: `Sample_ID`, `MTX` (methotrexate), `GC` (glucocorticoid), `TAC` (tacrolimus), `SSZ` (sulfasalazine), `bts` (biologic or targeted synthetic DMARD), `IS_csDMARD` (immunosuppressant or csDMARD other than methotrexate) | 05, 08, 08d, 10, 12 |
| `article/41plex_healthy_control.xlsx` | serum cytokines of the healthy reference cohort (n = 101) | 08e, figures |
| CT DICOM images | baseline and follow-up chest CT (not shared; the GMM fractions are in the per-patient tables) | `CT_GMM/lung_gmm_batch.py` |
