# BAL multi-omics in rheumatoid arthritis: analysis scripts

The scripts produce every figure and table of the article on the 31 analysed patients
(RA 24, sarcoidosis 7). The table **Article item → script → output** at the end lists,
for every figure and table, the code that produces it.

## Directory layout

```
scripts/
├── README.md                          # this file
├── install_dependencies.R             # CRAN / Bioconductor / GitHub packages
├── 00_cohort_config.R                 # cohort table, output folders, helpers (sourced by every script)
├── 01a_reference_atlas.R              # scRNA-seq reference atlas from public GEO data
├── 01b_BayesPrism_n31.R               # BayesPrism deconvolution of the 31 samples
├── 02_… 18_*.R                        # analyses: tables and values of the article
├── figures/                           # figure panels
│   └── panel_map.json                 # pixel size of each panel
└── CT_GMM/
    ├── lung_gmm_batch.py              # CT GMM tissue classification (used for the article)
    └── lung_gmm_interactive.ipynb     # optional 3D Slicer prototype
```

## Required software

- **R 4.5.0** (tested; macOS). Packages: `../R_packages.txt`, installer `install_dependencies.R`.
- **Python**, two separate environments (`../requirements.txt`): the CT GMM pipeline
  (Python 3.12, scikit-learn 1.8.0, numpy 1.26.4, SimpleITK, pydicom, lungmask/PyTorch; versions in `../requirements.txt`)
  and the MOFA2 backend `mofapy2` (interpreter given by the environment variable `MOFA_PYTHON`).
- About 32 GB RAM for `01a_reference_atlas.R` and `01b_BayesPrism_n31.R`.

The package versions of the
environment in which every script was run and checked against the article are in
`../docs/package_versions.txt`; the checks (`stopifnot` against the values of the article) can stop
under other versions. Results that can depend on the environment:
- KEGG GSEA (02, 03) queries the KEGG REST API at run time (the reported tables were retrieved in
  October–November 2026), so a later run may differ by a few terms.
- BayesPrism (01b) uses `detectCores() - 1` cores; its Gibbs sampling, and therefore the cell
  fractions, can differ slightly between machines. The fractions used in the article are the
  `bayesprism_n31/` outputs.
- PERMANOVA P values (16) depend on the random-number stream of the installed vegan/permute versions;
  R² and F do not.

The main scripts write `sessionInfo()`.

## Conventions

- Run from the **project root** (working directory, or `PROJECT_ROOT`). Every script sources
  `00_cohort_config.R`, which reads the cohort table `cohort_definition.csv`, defines the output
  folders and checks that the controlled-access tables carry English column names. Group
  membership, exclusions and the de-identified labels (`RA#`, `Sarcoidosis#`) are taken only from
  this table (`sample_group()`, `is_sarc()`, `public_label()`).
- Outputs go to `<ART_DIR>/tables/n31/` and `<ART_DIR>/figures/n31/panels/` (`ART_DIR` defaults to
  `article/` under the project root).

## Inputs

Individual-level inputs are available on request (`../data_examples/README.md`).

| File | Content | Used by |
|---|---|---|
| `results/RA_ILD_Workspace.RData` | processed data: per-patient table, flow cytometry, cytokine and microbiome matrices, bulk RNA-seq counts and gene annotation | 02–06, 08e, 15–17, figures |
| `results/tables/master_data_with_CT_all.csv` | per-patient table (English column names) | 05–10, 12, 14, 17, figures |
| `RNAseq/RawCount_heatmap_table.csv` | bulk BAL RNA-seq counts (GEO GSE329884) | 01b |
| `OTU abundance table-2 (merged abundance table).csv` | 16S OTU table (SRA PRJNA1462027) | 18 |
| `FCM_integrated_data_transformed.xlsx` | BAL / PB flow cytometry and BAL differential counts | 11, `figures/panels_deconv_FCM_serum.R` |
| `output_BayesPrism/BAL_reference_author_annotated.rds` | scRNA-seq reference (written by 01a from public data) | 01b, `figures/panels_fig1_atlas.R` |
| `<ART_DIR>/cohort_definition.csv` | Sample_ID, group, status (`analysed` for the 31 patients), public_id | all |
| `<ART_DIR>/tables/n31/planA/planA_data_RA24.csv` | one row per RA patient: baseline variables, infection outcome, pulmonary function nearest to BAL and to the follow-up CT, CT tissue fractions; times as days from BAL | most scripts |
| `<ART_DIR>/tables/n31/followup_by_sample.csv`, `infection_events.csv` | days from BAL to the last visit and to each infection | 08b |
| `<ART_DIR>/tables/n31/pft_series_by_sample.tsv` | pulmonary-function test series, days from BAL | 10 |
| `<ART_DIR>/tables/n31/sarc_BAL_WBC.csv` | BAL total cell count of the sarcoidosis patients | 10 |
| `<ART_DIR>/tables/n31/medication_at_BAL.csv` | medication in use at BAL (0/1), RA patients | 05, 08, 08d, 10, 12 |
| `<ART_DIR>/41plex_healthy_control.xlsx` | healthy serum reference (n = 101) | 08e, `figures/panels_deconv_FCM_serum.R` |

## Execution order

```bash
python  scripts/CT_GMM/lung_gmm_batch.py                 # CT tissue classification (DICOM + lungmask)
Rscript scripts/01a_reference_atlas.R                     # heavy, public data
Rscript scripts/01b_BayesPrism_n31.R                      # heavy
Rscript scripts/11_deconvolution_FCM_validation_n31.R
Rscript scripts/02_transcriptome_RA_vs_sarcoidosis.R
Rscript scripts/03_DEG_infection_within_RA.R
Rscript scripts/04_disease_MOFA2.R
Rscript scripts/05_RA_only_MOFA2.R
Rscript scripts/09_CT_total_opacity.R
Rscript scripts/17_GSVA_GO_scores.R                       # input of 06, 08, 10 and figures
Rscript scripts/06_common_panel_screens.R
Rscript scripts/06b_Fig3b_all_variables.R                  # Source Data of Figure 3b (all variables)
Rscript scripts/07_pulmonary_function_panels.R
Rscript scripts/08_sensitivity_AUC_stability_Firth.R
Rscript scripts/08b_followup_and_Cox.R
Rscript scripts/08e_RA_vs_sarcoidosis_screen_serum.R      # input of 08c
Rscript scripts/08c_SD2_effect_sizes.R
Rscript scripts/08d_ST3_baseline_by_infection.R
Rscript scripts/10_supplementary_table_values.R
Rscript scripts/12_MOFA_AUC_DeLong.R
Rscript scripts/13_ST8_CT_vs_PFT.R
Rscript scripts/14_SD5_profiles_by_ILD.R
Rscript scripts/15_Th_subsets_BAL_PB.R
Rscript scripts/16_PERMANOVA_infection_n31.R
Rscript scripts/18_SD4_microbiome_genus_exact.R
for f in panels_fig1_atlas panels_deconv_FCM_serum panels_transcriptome panels_screen_MOFA panels_fig3a_3b_3h \
         panels_infection_markers panels_infection_transcripts panel_sfig3c; do Rscript scripts/figures/$f.R; done
```

## Article item → script → output

Tables are written to `<ART_DIR>/tables/n31/`, panels to `<ART_DIR>/figures/n31/panels/` as
`<panel>.png` and `.pdf` (for example `Fig2e.png`).

| Item | Content | Script | Output |
|---|---|---|---|
| Table 1 | participant characteristics | not released (note 1) | — |
| Fig. 1a | study design | illustration | — |
| Fig. 1b–e | scRNA-seq reference atlas | `01a`, `figures/panels_fig1_atlas.R` | panels |
| Fig. 1f–i | BayesPrism fractions; deconvolution vs BAL differential counts | `01b`, `11`, `figures/panels_deconv_FCM_serum.R` | `bayesprism_n31/Fig1ghi_SFig1_deconv_vs_differential.csv` |
| Fig. 2a–d, 2g | PCA, DESeq2, GSEA, Ig genes | `02`, `figures/panels_transcriptome.R` | `transcriptome_*.csv`; `../results/` |
| Fig. 2e | RA vs sarcoidosis screen | `08e`, `08c`, `figures/panels_screen_MOFA.R` | `supp_tables/SD2_effect_sizes_116.csv` |
| Fig. 2f | BayesPrism plasma cell fraction | `11`, `figures/panels_deconv_FCM_serum.R` | `bayesprism_n31/Fig2f_bp_cell_plasma_Bcell.csv` |
| Fig. 2h–j | BAL CD4+ subsets; BAL vs PB | `15`, `figures/panels_screen_MOFA.R` | `Fig2h_BAL_Th_subsets.csv`, `Fig2ij_BAL_vs_PB*.csv` |
| Fig. 3a | PCA by infection; PERMANOVA | `16`, `figures/panels_fig3a_3b_3h.R` | `PERMANOVA_infection_RA24.csv` |
| Fig. 3b | infection screen, 114 markers (12 nominal); all compared variables | `06`, `06b`, `figures/panels_fig3a_3b_3h.R` | `common_panel/infection_panel.csv`, `common_panel/Fig3b_all_variables.csv` |
| Fig. 3c–f | BAL Th17.1, Treg, BAL/PB ratio, ROC | `figures/panels_infection_markers.R` | `panels/infection_markers_values.csv` |
| Fig. 3g | surfactant / epithelial transcripts | `03`, `figures/panels_infection_transcripts.R` | `DEG_Infection_pos_vs_neg.csv` |
| Fig. 3h | OXPHOS GSVA score | `17`, `figures/panels_fig3a_3b_3h.R` | `gsva_GO/gsva_GO_scores_n31.csv`, `gsva_GO/oxphos_gsva_stats.csv` |
| Fig. 4a–e, 4g–j | %FVC, %DLCO, INBUILD ROC, %FVC screen, CT vs pulmonary function | `07` (screen from `06`, opacity from `09`) | `panels/panel_values_pft.csv`, `common_panel/dFVC_panel.csv` |
| Fig. 4f | GMM segmentation images | `CT_GMM/lung_gmm_batch.py` | verification images |
| Fig. 5a–c | disease MOFA2 | `04`, `12`, `figures/panels_screen_MOFA.R` | `diseaseMOFA2_*.csv` |
| Fig. 5d–f | RA-only MOFA2, infection | `05` | `mofa_RA/*` |
| Fig. 5g–i | RA Factor 1 vs %FVC and total opacity; INBUILD ROC | `07`, `12` | `panels/panel_values_pft.csv`, `MOFA_AUC_DeLong.csv` |
| Suppl. Fig. 1 | Bland–Altman | `11`, `figures/panels_deconv_FCM_serum.R` | `bayesprism_n31/Fig1ghi_SFig1_deconv_vs_differential.csv` |
| Suppl. Fig. 2 | serum cytokines vs healthy reference | `08e`, `figures/panels_deconv_FCM_serum.R` | `ST1_serum_vs_healthy.csv` |
| Suppl. Fig. 3a, b | infection GSEA; microbiome genera | `03`, `18`, `figures/panels_infection_transcripts.R`, `figures/panels_infection_markers.R` | `GSEA_Infection_GO_BP.csv`, `microbiome/SD4_microbiome_genus_exact.csv` |
| Suppl. Fig. 3c | BAL Th17.1 by infection within ILD− / ILD+ | `figures/panel_sfig3c.R` | panel |
| Suppl. Fig. 4 | medication | `05` | `mofa_RA/SFig4_medication_n31.csv` |
| Suppl. Table 1 | 41 serum cytokines: RA, sarcoidosis, healthy reference | `08e` | `ST1_serum_vs_healthy.csv` |
| Suppl. Table 2 | infection events | chart review; no analysis code | — |
| Suppl. Table 3 | clinical characteristics by infection | `08d` | `ST3_by_infection.tsv` |
| Suppl. Table 4 | BAL CD4+ subset counts | `10` | `supp_tables/ST4_BAL_CD4_counts.csv` |
| Suppl. Table 5 | Firth adjustment, one covariate at a time | `08` | `sensitivity/ST5_Firth_adjusted.csv` |
| Suppl. Table 6 | clinical variables vs BAL and blood markers | `10` | `supp_tables/ST6_clinical_vs_multiomics.csv` |
| Suppl. Table 7 | partial Spearman correlations of the BAL Th17.1 count | `10` | `supp_tables/ST7_partial_correlations.csv` |
| Suppl. Table 8 | CT measures and pulmonary function | `13` | `ST8/ST8_CT_vs_PFT.csv` |
| Results (no table) | AUC stability and permutation P of the BAL Th17.1 proportion; follow-up and Cox model | `08`, `08b` | `sensitivity/AUC_stability.csv`, `sensitivity/AUC_permutation_two_sided.csv` |
| Results (no table) | SFTPB/ETV5 adjusted for the BAL differential count (DESeq2) | `03` | `DEG_Infection_adjusted_BAL_differential.csv` |
| Results (no table) | OXPHOS GSVA score adjusted for the BAL differential count | `08` | `sensitivity/OXPHOS_adjusted_for_BAL_differential.csv` |
| Results (no table) | BAL microbiome composition by infection (Bray–Curtis PERMANOVA, OTU rule of Supplementary Data 4); Shannon index | `16`; `06` | `PERMANOVA_microbiome_infection_RA24.csv`; `common_panel/infection_panel.csv` |
| Results (no table) | Shannon and Simpson indices computed from the OTU table | `16` | `microbiome_alpha_diversity_n31.csv` |
| Results (no table) | BAL Th17.1 count for an absolute %FVC decline of 5 points or more (AUC) | `07` | `panels/panel_values_pft.csv` (Decline5) |
| Discussion | sample size of a validation study (AUC 0.70, 80% power) | `08` | `sensitivity/validation_sample_size.csv` |
| Suppl. Data 1 | DEG RA vs sarcoidosis | `02` | `../results/DEG_RA_vs_sarcoidosis_n31_genome_wide.csv` |
| Suppl. Data 2 | RA vs sarcoidosis effect sizes | `08e` → `08c` | `supp_tables/SD2_effect_sizes_116.csv`, `supp_tables/SD2_reference_CD4_CD8_ratio.csv` |
| Suppl. Data 3 | infection screen (114 markers) and derived rows | `17` → `06`; transcripts from `03` | `common_panel/infection_panel.csv`, `infection_extra.csv` |
| Suppl. Data 4 | microbiome genera by infection | `18` (note 2) | `microbiome/SD4_microbiome_genus_exact.csv` |
| Suppl. Data 5 | profiles by ILD status | `14` | `supp_tables/SD5_profiles_by_ILD.csv` |
| Suppl. Data 6 | %FVC screen | `06` | `common_panel/dFVC_panel.csv` |
| Suppl. Data 7 | total opacity screen | `06` | `common_panel/CT_panel.csv` |
| Suppl. Data 8 | disease MOFA2 weights | `04` | `diseaseMOFA2_weights_all.csv` |
| Suppl. Data 9 | RA-only MOFA2 weights | `05` | `mofa_RA/SD9_weights_n31.csv` |

Notes:

1. Table 1 was assembled from the clinical records; its values can be recomputed from the
   individual-level data available on request.
2. Supplementary Data 4: OTUs with ≥ 5 reads in ≥ 3 of the 24 RA samples, genus = 6th taxonomy
   level, exact Wilcoxon P, BH across 106 genera.
3. Figure panels are drawn at the pixel size of the article (`figures/panel_map.json`). Figure
   1b, 2h, 3c–f and Supplementary Fig. 3b of the article were drawn by an earlier version of the
   same code; the scripts draw the same data with the same labels, and point jitter and margins
   differ slightly.

## Statistical conventions

See the Statistical conventions of the top-level README.md.

## Notes

- Raw sequences: GEO GSE329884 (RNA-seq) and SRA PRJNA1462027 (16S).
- No patient identifier, name, date of birth, record number or calendar date appears in the
  code; samples are referred to through `cohort_definition.csv`.
