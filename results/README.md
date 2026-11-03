# Genome-wide result tables

Gene-level and gene-set-level tables of the article. They contain no
patient-level data. The shrunken log2 fold changes (`DESeq2::lfcShrink`, `type = "normal"`) are taken
from a fit on the Cook's-outlier-replaced counts so that their sign agrees with
the Wald statistic, while `baseMean`, `stat`, `pvalue` and `padj` are those of
the default `DESeq()` fit.

| File | Contrast | Script | Rows |
|---|---|---|---|
| `DEG_RA_vs_sarcoidosis_n31_genome_wide.csv` | RA (n = 24) vs sarcoidosis (n = 7); design `~Group`, genes with ≥ 10 counts in ≥ 5 samples | `scripts/02_transcriptome_RA_vs_sarcoidosis.R` | 24,353 genes; 143 at padj < 0.05 |
| `DEG_infection_within_RA_n24_genome_wide.csv` | infection (+, n = 6) vs (−, n = 18) within RA; design `~Infection`, genes with ≥ 10 counts in ≥ 3 samples | `scripts/03_DEG_infection_within_RA.R` | 25,061 genes; 2,024 at padj < 0.05 |
| `GSEA_RA_vs_sarcoidosis_GO_BP.csv` | GSEA (GO biological process) on the shrunken LFC, RA vs sarcoidosis; `pvalueCutoff = 0.05` | `02_transcriptome_RA_vs_sarcoidosis.R` | 547 terms |
| `GSEA_RA_vs_sarcoidosis_KEGG.csv` | GSEA (KEGG) on the shrunken LFC, RA vs sarcoidosis; `pvalueCutoff = 0.05` | `02_transcriptome_RA_vs_sarcoidosis.R` | 36 terms |
| `GSEA_infection_within_RA_GO_BP.csv` | GSEA (GO biological process), infection (+) vs (−); `pvalueCutoff = 0.1` | `03_DEG_infection_within_RA.R` | 543 terms |
| `GSEA_infection_within_RA_KEGG.csv` | GSEA (KEGG), infection (+) vs (−); `pvalueCutoff = 0.1` | `03_DEG_infection_within_RA.R` | 98 terms |

## Script output → released file

The released files are byte-identical copies of these outputs of the scripts in `scripts/` (`<ART_TABLES>/n31/`):

| Script output | Released file |
|---|---|
| `transcriptome_DEG_ShrunkLFC.csv` (02) | `DEG_RA_vs_sarcoidosis_n31_genome_wide.csv` |
| `transcriptome_GSEA_GO_BP.csv` (02) | `GSEA_RA_vs_sarcoidosis_GO_BP.csv` |
| `transcriptome_GSEA_KEGG.csv` (02) | `GSEA_RA_vs_sarcoidosis_KEGG.csv` |
| `DEG_Infection_pos_vs_neg.csv` (03) | `DEG_infection_within_RA_n24_genome_wide.csv` |
| `GSEA_Infection_GO_BP.csv` (03) | `GSEA_infection_within_RA_GO_BP.csv` |
| `GSEA_Infection_KEGG.csv` (03) | `GSEA_infection_within_RA_KEGG.csv` |

## Columns

DESeq2 tables: `baseMean`, `log2FoldChange` (shrunken; positive = higher in RA,
or higher in the infection (+) group), `lfcSE`, `stat` (Wald statistic),
`pvalue`, `padj` (Benjamini–Hochberg), `Gene_Symbol`, `ENTREZID`, `GENENAME`,
`Regulation` (`Up`/`Down` or `Up_in_InfPos`/`Down_in_InfPos` at |LFC| > 0.585
and padj < 0.05, otherwise `NS`).

GSEA tables (`clusterProfiler::gseGO` / `gseKEGG`, ranking by shrunken LFC,
BH, gene-set size 15–500; `set.seed(42)` is called before `gseGO`, and `gseKEGG` continues the same random-number stream): `ID`, `Description`, `setSize`,
`enrichmentScore`, `NES`, `pvalue`, `p.adjust`, `qvalue`, `rank`,
`leading_edge`, `core_enrichment` (ENTREZ IDs). KEGG annotations are fetched
from the KEGG REST API at run time, so a later run may differ by a few terms.
