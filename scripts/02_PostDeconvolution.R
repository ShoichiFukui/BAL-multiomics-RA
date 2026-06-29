# ============================================================================
# RA-ILD Multi-omics Integrated Analysis Pipeline
# Primary Files Only - For manuscript submission
# ============================================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║    RA-ILD Multi-omics Analysis - Complete Pipeline        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Environment Setup
suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "readxl", "DESeq2", "edgeR", "clusterProfiler", 
            "enrichplot", "org.Hs.eg.db", "GSVA", "vegan",
            "pheatmap", "randomForest", "caret", "pROC")
  for(p in pkgs) {
    if(!requireNamespace(p, quietly=TRUE)) {
      if(p %in% c("DESeq2","edgeR","clusterProfiler","enrichplot","org.Hs.eg.db","GSVA")) {
        BiocManager::install(p, ask=FALSE, update=FALSE)
      } else install.packages(p, repos="https://cloud.r-project.org/")
    }
    library(p, character.only=TRUE)
  }
})

set.seed(42)
options(stringsAsFactors=FALSE)
for(d in c("results","results/figures","results/tables","results/consensus_clustering")) 
  dir.create(d, showWarnings=FALSE)

color_disease <- c("Control"="#4DAF4A", "RA"="#E41A1C")
color_infection <- c("Infection_Negative"="#377EB8", "Infection_Positive"="#FF7F00")

safe_write <- function(data, file, ...) {
  tryCatch({write.csv(data, file, ...); cat(sprintf("  ✓ %s\n", basename(file)))},
           error=function(e) cat(sprintf("  ✗ %s\n", basename(file))))
}

safe_pdf <- function(file, expr, ...) {
  tryCatch({pdf(file, ...); eval(expr); dev.off(); cat(sprintf("  ✓ %s\n", basename(file)))},
           error=function(e) {if(dev.cur()>1) dev.off(); cat(sprintf("  ✗ %s\n", basename(file)))})
}

# Data Loading
cat("\n=== Data Loading ===\n")
data_dir <- "."

rnaseq_raw <- read.csv(file.path(data_dir, "RNAseq/RawCount_heatmap_table.csv"), 
                       row.names=1, check.names=FALSE)
if("Cluster_ID" %in% colnames(rnaseq_raw)) rnaseq_raw <- rnaseq_raw[, !colnames(rnaseq_raw) %in% "Cluster_ID"]

microbiome_raw <- read.csv(file.path(data_dir, "OTU abundance table-2 (merged abundance table).csv"),
                           check.names=FALSE)
sample_cols <- grep("^(RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)
new_names <- gsub("^((?:RA|Sarcoidosis)\\d+).*", "\\1", sample_cols)
colnames(microbiome_raw)[match(sample_cols, colnames(microbiome_raw))] <- new_names
if(any(duplicated(new_names))) {
  microbiome_agg <- microbiome_raw[, c("ID","Name","Taxonomy")]
  for(s in unique(new_names)) {
    cols <- which(colnames(microbiome_raw)==s)
    microbiome_agg[[s]] <- if(length(cols)==1) microbiome_raw[[cols]] else rowMeans(microbiome_raw[,cols])
  }
  microbiome_raw <- microbiome_agg
}

# NOTE: rename the source clinical/microbiome data file to this ASCII filename to match.
clinical_data <- read_excel(file.path(data_dir, "clinical_metadata.xlsx"), sheet="Analysis")
colnames(clinical_data)[1] <- "Sample_ID"
colnames(clinical_data) <- gsub(" ", "_", colnames(clinical_data))
clinical_colnames <- c("Sample_ID","respiratory_infection","Age","Age_over_65","ra_diagnosis",
                       "follow_up_date","gender","male","smoking_binary","smoking_status","tobacco_per_day","smoking_years",
                       "passive_smoking","pneumococcal_vaccine","antibiotic_use","steroid_use","immunosuppressant",
                       "steroid_or_immunosuppressant","biologics","ct_uip","ct_emphysema","ct_bronchodilation",
                       "mycobacterium_positive","aspergillus_positive","performance_status","height","weight","bmi",
                       "bmi_original","ra_duration_years","nsaids","sasp","dmards","methotrexate","hematologic_malignancy",
                       "solid_malignancy","stroke_history","anti_ccp_antibody","crp","crp_under_1","esr_1hr","esr_under_40",
                       "das28_crp","das28_crp_category","das28_esr","das28_esr_category","sdai","sdai_category","cdai",
                       "cdai_category","das28_crp_binary","das28_esr_binary","sdai_binary","cdai_binary","anti_ccp_over_100",
                       "ra_over_3years","ra_3years_or_more","bal_wbc","bal_total_cells","bal_macrophage_percent",
                       "bal_macrophage_count","bal_lymphocyte_percent","bal_lymphocyte_count","bal_neutrophil_percent",
                       "bal_neutrophil_count","bal_eosinophil_percent","bal_basophil_percent","bal_cd25","bal_treg_percent",
                       "bal_treg_count","bal_th_cells","bal_th17_percent","bal_th17_count","bal_th2_percent","bal_th1_percent",
                       "bal_macrophage_subtype","bal_activated_macrophage","bal_neutrophil_cd66b_cd16_positive",
                       "bal_neutrophil_cd66b_cd206_negative","bal_treg_th17_ratio","blood_treg","blood_th17","blood_treg_th17_ratio")
if(ncol(clinical_data)>=length(clinical_colnames)) 
  colnames(clinical_data)[1:length(clinical_colnames)] <- clinical_colnames

ct_data <- read_excel(file.path(data_dir, "RA_CT_transformed.xlsx"), sheet="analysis")
balf_cytokine <- read_excel(file.path(data_dir, "cytokine_multiplex.xlsx"), sheet="BALF_Analysis")
serum_cytokine <- read_excel(file.path(data_dir, "cytokine_multiplex.xlsx"), sheet="Serum_Analysis")
balf_fcm <- read_excel(file.path(data_dir, "FCM_integrated_data_transformed.xlsx"), sheet="BALF_analysis")
pbmc_fcm <- read_excel(file.path(data_dir, "FCM_integrated_data_transformed.xlsx"), sheet="Blood_analysis")
colnames(balf_cytokine)[1] <- colnames(serum_cytokine)[1] <- "Sample_ID"


# ============================================================================
# BayesPrism Deconvolution Data Loading
# ============================================================================
deconv_loaded <- FALSE
bp_theta <- NULL
bp_cell <- NULL

cat("\n=== Loading BayesPrism Deconvolution Results ===\n")

bp_dir <- file.path(data_dir, "output_v3_BayesPrism")

# PRIMARY: Try RDS file for theta (RNA fraction) + CSV for cell fraction
rds_path <- file.path(bp_dir, "full_analysis_results_v3.rds")
if(file.exists(rds_path)) {
  cat("Loading from RDS file...\n")
  tryCatch({
    v3_results <- readRDS(rds_path)
    bp_theta <- v3_results$theta_type
    deconv_loaded <- TRUE
    cat(sprintf("  ✓ bp_theta (RNA fraction): %d samples × %d cell types\n",
                nrow(bp_theta), ncol(bp_theta)))
    cat(sprintf("  Cell types: %s\n", paste(colnames(bp_theta), collapse=", ")))
    if("Plasma" %in% colnames(bp_theta)) {
      cat(sprintf("  ✓ Plasma confirmed (mean: %.4f%%)\n", mean(bp_theta[, "Plasma"]) * 100))
    }

    # Cell fraction (RNA content corrected) from CSV
    csv_cell <- file.path(bp_dir, "celltype_proportions_cellFraction.csv")
    if(file.exists(csv_cell)) {
      bp_cell <- as.matrix(read.csv(csv_cell, row.names=1, check.names=FALSE))
      cat(sprintf("  ✓ bp_cell (cell fraction): %d samples × %d cell types\n",
                  nrow(bp_cell), ncol(bp_cell)))
    } else {
      cat("  ⚠ Cell fraction CSV not found, using RNA fraction as fallback\n")
      bp_cell <- bp_theta
    }
  }, error = function(e) {
    cat(sprintf("  ✗ Error loading RDS: %s\n", e$message))
    deconv_loaded <- FALSE
  })
}

# FALLBACK: Try CSV files
if(!deconv_loaded) {
  cat("Falling back to CSV files...\n")
  for(search_dir in c(bp_dir, "./output_v3_BayesPrism", data_dir)) {
    csv_theta <- file.path(search_dir, "celltype_proportions_BayesPrism.csv")
    if(file.exists(csv_theta)) {
      bp_theta <- as.matrix(read.csv(csv_theta, row.names=1, check.names=FALSE))
      csv_cell <- file.path(search_dir, "celltype_proportions_cellFraction.csv")
      bp_cell <- if(file.exists(csv_cell)) {
        as.matrix(read.csv(csv_cell, row.names=1, check.names=FALSE))
      } else {
        bp_theta
      }
      deconv_loaded <- TRUE
      cat(sprintf("  ✓ CSV loaded: %d samples × %d cell types\n",
                  nrow(bp_theta), ncol(bp_theta)))
      if(!"Plasma" %in% colnames(bp_theta)) {
        cat("  ⚠ Plasma not in CSV data (likely 9 cell type version)\n")
      }
      break
    }
  }
}

if(!deconv_loaded) {
  cat("  ✗ No deconvolution data found\n")
  cat("  Deconvolution analysis will be skipped\n")
}

exclude_samples <- character(0)
rnaseq_raw <- rnaseq_raw[, !colnames(rnaseq_raw) %in% exclude_samples]
microbiome_raw <- microbiome_raw[, !colnames(microbiome_raw) %in% exclude_samples]
clinical_data <- clinical_data[!clinical_data$Sample_ID %in% exclude_samples,]
ct_data <- ct_data[!ct_data$Sample_ID %in% exclude_samples,]
balf_cytokine <- balf_cytokine[!balf_cytokine$Sample_ID %in% exclude_samples,]
serum_cytokine <- serum_cytokine[!serum_cytokine$Sample_ID %in% exclude_samples,]
balf_fcm <- balf_fcm[!balf_fcm$Sample_ID %in% exclude_samples,]
pbmc_fcm <- pbmc_fcm[!pbmc_fcm$Sample_ID %in% exclude_samples,]
if(deconv_loaded) {bp_theta <- bp_theta[!rownames(bp_theta) %in% exclude_samples,]
bp_cell <- bp_cell[!rownames(bp_cell) %in% exclude_samples,]}

all_samples <- unique(c(colnames(rnaseq_raw), grep("^(RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)))
master_data <- data.frame(Sample_ID=all_samples, 
                          Sample_Group=ifelse(grepl("^Sarcoidosis",all_samples),"Control","RA")) %>%
  left_join(clinical_data, by="Sample_ID") %>% left_join(ct_data, by="Sample_ID")

balf_renamed <- balf_cytokine
cytokine_cols <- setdiff(colnames(balf_cytokine), c("Sample_ID","Disease","Infection"))
colnames(balf_renamed)[colnames(balf_renamed) %in% cytokine_cols] <- paste0("BALF_", cytokine_cols)
serum_renamed <- serum_cytokine
serum_cols <- setdiff(colnames(serum_cytokine), c("Sample_ID","Disease","Infection"))
colnames(serum_renamed)[colnames(serum_renamed) %in% serum_cols] <- paste0("Serum_", serum_cols)

master_data <- master_data %>%
  left_join(balf_renamed %>% dplyr::select(-any_of(c("Disease","Infection"))), by="Sample_ID") %>%
  left_join(serum_renamed %>% dplyr::select(-any_of(c("Disease","Infection"))), by="Sample_ID") %>%
  left_join(balf_fcm, by="Sample_ID") %>% left_join(pbmc_fcm, by="Sample_ID")

if("respiratory_infection" %in% colnames(master_data)) {
  master_data$Infection <- master_data$respiratory_infection
  master_data$Infection_Group <- factor(
    ifelse(master_data$Infection==1,"Infection_Positive","Infection_Negative"),
    levels=c("Infection_Negative","Infection_Positive"))
}
if("Total_Score" %in% colnames(master_data)) {
  master_data$ILD_Score <- as.numeric(master_data$Total_Score)
  
  # Redefine the ILD column based on Total_Score (overwrites the old ILD column from ct_data)
  master_data$ILD <- ifelse(master_data$ILD_Score > 0, 1, 0)

  # ILD_Group: Total_Score > 0 is ILD-positive
  master_data$ILD_Group <- factor(
    ifelse(master_data$ILD_Score > 0, "ILD_Positive", "ILD_Negative"),
    levels = c("ILD_Negative", "ILD_Positive"))
  
  # ILD subgroups within RA (for 3-group comparison)
  master_data$Subgroup <- dplyr::case_when(
    master_data$Sample_Group == "Control" ~ "Control",
    master_data$ILD_Score > 0             ~ "RA_ILD",
    TRUE                                  ~ "RA_nonILD"
  )
  master_data$Subgroup <- factor(master_data$Subgroup, 
                                  levels = c("Control", "RA_nonILD", "RA_ILD"))
  
  cat(sprintf("  ILD classification (Total_Score > 0):\n"))
  cat(sprintf("    RA-ILD: %d, RA-nonILD: %d, Control: %d\n",
              sum(master_data$Subgroup == "RA_ILD"),
              sum(master_data$Subgroup == "RA_nonILD"),
              sum(master_data$Subgroup == "Control")))
}

safe_write(master_data, "results/tables/master_data.csv", row.names=FALSE)
cat("✓ Data loaded\n")

# Gene Annotation
cat("\n=== Gene Annotation ===\n")
gene_annotations <- data.frame(Original_ID=rownames(rnaseq_raw), 
                               ENSEMBL=gsub("\\..*","",rownames(rnaseq_raw))) %>%
  left_join(tryCatch(AnnotationDbi::select(org.Hs.eg.db, keys=gsub("\\..*","",rownames(rnaseq_raw)),
                                           columns=c("SYMBOL","GENENAME","ENTREZID"), keytype="ENSEMBL") %>%
                       group_by(ENSEMBL) %>% slice(1) %>% ungroup(), 
                     error=function(e) data.frame()), by="ENSEMBL") %>%
  mutate(Gene_Symbol=ifelse(is.na(SYMBOL), Original_ID, SYMBOL))
dup_names <- gene_annotations$Gene_Symbol[duplicated(gene_annotations$Gene_Symbol)]
for(dn in unique(dup_names)) {
  idx <- which(gene_annotations$Gene_Symbol==dn)
  gene_annotations$Gene_Symbol[idx[-1]] <- paste0(gene_annotations$Gene_Symbol[idx[-1]],"_",
                                                  gene_annotations$ENSEMBL[idx[-1]])
}
rnaseq_annotated <- rnaseq_raw; rownames(rnaseq_annotated) <- gene_annotations$Gene_Symbol
safe_write(gene_annotations, "results/tables/gene_annotations.csv", row.names=FALSE)
cat("✓ Genes annotated\n")

# Expression Matrices
cat("\n=== Expression Matrices ===\n")
common_samples <- intersect(colnames(rnaseq_annotated), master_data$Sample_ID)
count_matrix <- as.matrix(rnaseq_annotated[, common_samples])
count_matrix <- count_matrix[rowSums(count_matrix)>0,]
sample_metadata <- master_data[match(common_samples, master_data$Sample_ID),]
rownames(sample_metadata) <- sample_metadata$Sample_ID

dds <- DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_metadata, design=~Sample_Group)
dds <- dds[rowSums(counts(dds)>=10)>=5,]; dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE); expr_matrix <- assay(vsd)
log_cpm <- cpm(count_matrix[rownames(expr_matrix),], log=TRUE, prior.count=1)

gene_vars <- apply(expr_matrix, 1, var)
top_genes <- names(sort(gene_vars, decreasing=TRUE))[1:min(500,length(gene_vars))]
expr_mat <- t(expr_matrix[top_genes, common_samples])

balf_cytokine_cols <- grep("^BALF_", colnames(master_data), value=TRUE)
balf_cytokine_cols <- balf_cytokine_cols[!grepl("Treg|Th|Macro|Neutro|Lympho|CD", balf_cytokine_cols)]
serum_cytokine_cols <- grep("^Serum_", colnames(master_data), value=TRUE)
all_cytokine_cols <- c(balf_cytokine_cols, serum_cytokine_cols)

cyto_data <- master_data[match(common_samples, master_data$Sample_ID), all_cytokine_cols]
cyto_mat <- log2(as.matrix(cyto_data)+1); rownames(cyto_mat) <- common_samples
colnames(cyto_mat) <- gsub("^BALF_|^Serum_", "", colnames(cyto_mat))

balf_fcm_cols <- grep("^BALF_.*(Treg|Th|Macro|Neutro)", colnames(master_data), value=TRUE)
pbmc_fcm_cols <- grep("^PB_", colnames(master_data), value=TRUE)
all_fcm_cols <- c(balf_fcm_cols, pbmc_fcm_cols)
fcm_data <- master_data[match(common_samples, master_data$Sample_ID), all_fcm_cols]
char_cols <- names(fcm_data)[sapply(fcm_data, is.character)]
if(length(char_cols)>0) for(col in char_cols) fcm_data[[col]] <- as.numeric(fcm_data[[col]])
numeric_cols <- sapply(fcm_data, is.numeric)
if(sum(!numeric_cols)>0) fcm_data <- fcm_data[, numeric_cols, drop=FALSE]
fcm_mat <- log2(as.matrix(fcm_data)+0.01); rownames(fcm_mat) <- common_samples

meta <- sample_metadata[common_samples, c("Sample_ID","Sample_Group","ILD_Group"), drop=FALSE]
if("ILD_Score" %in% colnames(sample_metadata)) {
  meta$ILD_Score <- sample_metadata[common_samples, "ILD_Score"]
}
if("Subgroup" %in% colnames(sample_metadata)) {
  meta$Subgroup <- sample_metadata[common_samples, "Subgroup"]
}
if("Infection_Group" %in% colnames(sample_metadata)) {
  meta$Infection_Group <- sample_metadata[common_samples, "Infection_Group"]
} else if("respiratory_infection" %in% colnames(sample_metadata)) {
  resp_inf <- sample_metadata[common_samples, "respiratory_infection"]
  meta$Infection_Group <- factor(ifelse(resp_inf==1,"Infection_Positive","Infection_Negative"),
                                 levels=c("Infection_Negative","Infection_Positive"))
}
cat("✓ Matrices prepared\n")

# DEG Analysis 
cat("\n=== DEG Analysis ===\n")
res <- results(dds, contrast=c("Sample_Group","RA","Control"))
res_df <- as.data.frame(res) %>% mutate(Gene_Symbol=rownames(res)) %>%
  left_join(gene_annotations[,c("Gene_Symbol","ENTREZID","GENENAME")], by="Gene_Symbol") %>% arrange(padj)
res_df$Regulation <- "NS"
res_df$Regulation[res_df$padj<0.05 & res_df$log2FoldChange>1] <- "Up"
res_df$Regulation[res_df$padj<0.05 & res_df$log2FoldChange < (-1)] <- "Down"
safe_write(res_df, "results/tables/DEG_RA_vs_Control.csv", row.names=FALSE)
cat(sprintf("  DEGs: %d up, %d down\n", sum(res_df$Regulation=="Up",na.rm=TRUE),
            sum(res_df$Regulation=="Down",na.rm=TRUE)))

# GSVA
cat("\n=== GSVA ===\n")
gsva_results <- NULL
tryCatch({
  gene_sets <- list(
    OXPHOS=c("NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFA10",
             "NDUFA11","NDUFA12","NDUFA13","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6",
             "NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3",
             "NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NDUFV3","SDHA","SDHB","SDHC",
             "SDHD","UQCRC1","UQCRC2","UQCRFS1","UQCRQ","UQCRH","UQCRB","CYC1","UQCR10","UQCR11","COX4I1",
             "COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7B","COX7C","COX8A","MT-CO1","MT-CO2",
             "MT-CO3","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB","ATP5PD","ATP5PF","ATP5PO",
             "ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5MG","MT-ATP6","MT-ATP8"),
    Mito_Biogenesis=c("TFAM","TFB1M","TFB2M","POLG","POLG2","POLRMT","MRPL1","MRPL2","MRPL3","MRPL4",
                      "MRPL10","MRPL11","MRPL12","MRPL13","MRPS2","MRPS5","MRPS6","MRPS7","MRPS10","MRPS11",
                      "MRPS12","PGC1A","PPARGC1A","PPARGC1B","PPRC1","NRF1","NRF2","GABPA","ESRRA","ESRRG",
                      "YY1","CREB1","SP1","MEF2A","MEF2D"),
    TCA_Cycle=c("CS","ACO1","ACO2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","OGDH","DLST","DLD","SUCLA2",
                "SUCLG1","SUCLG2","SDHA","SDHB","SDHC","SDHD","FH","MDH1","MDH2","PC","PCK1","PCK2","PDHA1",
                "PDHB","DLAT"),
    Surfactant=c("SFTPA1","SFTPA2","SFTPB","SFTPC","SFTPD","SFTA2","SFTA3","ABCA3","NKX2-1","NAPSA","LPCAT1",
                 "PGC","LAMP3","SLC34A2","MUC1","CLDN18","AQP5","AGER","ETV5","HOPX","SCGB1A1","SCGB3A1",
                 "SCGB3A2","BPIFA1","BPIFB1","PIGR"),
    Innate_Myeloid=c("CD14","CD68","CD163","MSR1","MRC1","MARCO","TLR2","TLR4","TLR7","TLR8","TLR9","MYD88",
                     "IRAK1","IRAK4","TRAF6","NLRP3","PYCARD","CASP1","IL1B","IL6","IL8","CXCL8","TNF","IL12A",
                     "IL12B","IL23A","CCL2","CCL3","CCL4","CCL5","CXCL9","CXCL10","CXCL11","S100A8","S100A9",
                     "S100A12","LYZ"),
    Complement=c("C1QA","C1QB","C1QC","C1R","C1S","C2","C3","C4A","C4B","C5","C6","C7","C8A","C8B","C8G",
                 "C9","CFB","CFD","CFH","CFI","CFHR1","CFHR3","CFHR5","CR1","CR2","CD46","CD55","CD59"),
    Adaptive_Tcell=c("CD3D","CD3E","CD3G","CD4","CD8A","CD8B","IL2","IL4","IL5","IL13","IL17A","IL17F","IL21",
                     "IL22","IFNG","TNF","GZMB","PRF1","TBX21","GATA3","RORC","FOXP3","BCL6","CCR4","CCR5","CCR6",
                     "CCR7","CXCR3","CXCR5","CTLA4","PDCD1","LAG3","TIGIT","HAVCR2","CD28","ICOS","CD40LG","CD27",
                     "CD69","CD25","IL2RA"),
    Cytotoxicity=c("PRF1","GZMA","GZMB","GZMH","GZMK","GZMM","NKG7","GNLY","KLRK1","KLRD1","KLRB1","FASLG",
                   "FAS","TNFSF10","TNFRSF10A","TNFRSF10B","IFNG","TNF","EOMES","TBX21","NCR1","NCR3","CD244","CD160"),
    Bcell_Humoral=c("CD19","MS4A1","CD79A","CD79B","CR2","IGHM","IGHD","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1",
                    "IGHA2","IGKC","IGLC1","IGLC2","IGLC3","PRDM1","XBP1","IRF4","PAX5","BCL6","CXCR5","CXCR4"),
    IFN_Response=c("ISG15","ISG20","MX1","MX2","OAS1","OAS2","OAS3","IFIT1","IFIT2","IFIT3","IFIT5","IFI6",
                   "IFI16","IFI27","IFI35","IFI44","IFI44L","IFITM1","IFITM2","IFITM3","IRF1","IRF3","IRF7",
                   "IRF9","STAT1","STAT2","JAK1","TYK2","RSAD2","HERC5","USP18"),
    TNFa_NFkB=c("TNFAIP3","TNFAIP6","NFKBIA","NFKB1","NFKB2","RELA","RELB","CCL2","CCL3","CCL4","CCL5","CCL20",
                "CXCL1","CXCL2","CXCL3","CXCL8","CXCL10","ICAM1","VCAM1","SELE","SELP","IL1B","IL6","IL8","TNF",
                "PTGS2","BIRC3","BCL2A1","SOD2"),
    Th17_Pathway=c("RORC","IL17A","IL17F","IL17RA","IL17RC","IL21","IL22","IL23R","IL23A","IL12B","CCR6","CCL20",
                   "CXCR6","STAT3","BATF","IRF4","RUNX1","TGFB1","TGFBR1","TGFBR2","IL6","IL6R","IL6ST","AHR","HIF1A"),
    Fibrosis_ECM=c("COL1A1","COL1A2","COL3A1","COL5A1","COL6A1","COL6A2","FN1","VIM","ACTA2","TAGLN","MYH11",
                   "TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","CTGF","PDGFA","PDGFB","PDGFRB","LOX","LOXL1","LOXL2",
                   "LOXL4","MMP1","MMP2","MMP9","MMP14","TIMP1","TIMP2"),
    Autophagy=c("ATG5","ATG7","ATG12","ATG16L1","BECN1","MAP1LC3A","MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2",
                "SQSTM1","NBR1","OPTN","CALCOCO2","ULK1","ULK2","RB1CC1","ATG13","ATG101","PIK3C3","PIK3R4","UVRAG",
                "PINK1","PRKN","MFN1","MFN2","OPA1","DNM1L"),
    Antigen_Presentation=c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","HLA-DRA","HLA-DRB1","HLA-DRB3",
                           "HLA-DRB4","HLA-DRB5","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","B2M","TAP1","TAP2",
                           "TAPBP","PSMB8","PSMB9","PSMB10","CD74","CIITA","NLRC5","RFX5")
  )
  
  gene_sets_filtered <- list()
  for(gs_name in names(gene_sets)) {
    genes_in_set <- intersect(gene_sets[[gs_name]], rownames(log_cpm))
    if(length(genes_in_set)>=5) gene_sets_filtered[[gs_name]] <- genes_in_set
  }
  
  gsva_mat <- log_cpm[,common_samples]; gsva_scores <- NULL; gsva_method <- "unknown"
  if(is.null(gsva_scores)) tryCatch({gsva_param <- gsvaParam(gsva_mat, gene_sets_filtered, kcdf="Gaussian", 
                                                             mx.diff=TRUE); gsva_scores <- gsva(gsva_param)
                                                             gsva_method <- "gsva_new_api"}, error=function(e) NULL)
  if(is.null(gsva_scores)) tryCatch({gsva_scores <- gsva(gsva_mat, gene_sets_filtered, method="gsva", 
                                                         kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE)
  gsva_method <- "gsva_legacy_api"}, error=function(e) NULL)
  if(is.null(gsva_scores)) tryCatch({gsva_scores <- gsva(gsva_mat, gene_sets_filtered, method="ssgsea", 
                                                         verbose=FALSE); gsva_method <- "ssgsea"}, error=function(e) NULL)
  if(is.null(gsva_scores)) {
    gsva_scores <- matrix(NA, nrow=length(gene_sets_filtered), ncol=ncol(gsva_mat))
    rownames(gsva_scores) <- names(gene_sets_filtered); colnames(gsva_scores) <- colnames(gsva_mat)
    for(i in seq_along(gene_sets_filtered)) {
      gs_genes <- intersect(gene_sets_filtered[[i]], rownames(gsva_mat))
      if(length(gs_genes)>0) gsva_scores[i,] <- colMeans(gsva_mat[gs_genes,,drop=FALSE], na.rm=TRUE)
    }
    gsva_scores <- t(scale(t(gsva_scores))); gsva_method <- "zscore_fallback"
  }
  
  ra_samples <- common_samples[meta$Sample_Group=="RA"]; ra_meta <- meta[ra_samples,]
  gsva_infection_comp <- data.frame()
  if("Infection_Group" %in% colnames(ra_meta)) {
    inf_neg <- ra_samples[!is.na(ra_meta$Infection_Group) & ra_meta$Infection_Group=="Infection_Negative"]
    inf_pos <- ra_samples[!is.na(ra_meta$Infection_Group) & ra_meta$Infection_Group=="Infection_Positive"]
    if(length(inf_neg)>=3 && length(inf_pos)>=3) {
      for(gs in rownames(gsva_scores)) {
        wt <- wilcox.test(gsva_scores[gs,inf_neg], gsva_scores[gs,inf_pos], exact=TRUE)
        gsva_infection_comp <- rbind(gsva_infection_comp, data.frame(
          GeneSet=gs, Mean_Inf_Neg=mean(gsva_scores[gs,inf_neg],na.rm=TRUE),
          Mean_Inf_Pos=mean(gsva_scores[gs,inf_pos],na.rm=TRUE),
          Delta=mean(gsva_scores[gs,inf_pos],na.rm=TRUE)-mean(gsva_scores[gs,inf_neg],na.rm=TRUE),
          P_value=wt$p.value, stringsAsFactors=FALSE))
      }
      gsva_infection_comp$P_adjusted <- p.adjust(gsva_infection_comp$P_value, method="BH")
      gsva_infection_comp <- gsva_infection_comp %>% arrange(P_value)
      safe_write(gsva_infection_comp, "results/tables/GSVA_infection_comparison.csv", row.names=FALSE)
    }
  }
  
  ctrl_samples <- common_samples[meta$Sample_Group=="Control"]
  ra_all_samples <- common_samples[meta$Sample_Group=="RA"]; gsva_disease_comp <- data.frame()
  if(length(ctrl_samples)>=3 && length(ra_all_samples)>=3) {
    for(gs in rownames(gsva_scores)) {
      wt <- wilcox.test(gsva_scores[gs,ctrl_samples], gsva_scores[gs,ra_all_samples], exact=TRUE)
      gsva_disease_comp <- rbind(gsva_disease_comp, data.frame(
        GeneSet=gs, Mean_Control=mean(gsva_scores[gs,ctrl_samples],na.rm=TRUE),
        Mean_RA=mean(gsva_scores[gs,ra_all_samples],na.rm=TRUE),
        Delta=mean(gsva_scores[gs,ra_all_samples],na.rm=TRUE)-mean(gsva_scores[gs,ctrl_samples],na.rm=TRUE),
        P_value=wt$p.value, stringsAsFactors=FALSE))
    }
    gsva_disease_comp$P_adjusted <- p.adjust(gsva_disease_comp$P_value, method="BH")
    gsva_disease_comp <- gsva_disease_comp %>% arrange(P_value)
    safe_write(gsva_disease_comp, "results/tables/GSVA_disease_comparison.csv", row.names=FALSE)
  }
  
  external_cor <- data.frame()
  if(ncol(fcm_mat)>0) {
    fcm_samples <- intersect(ra_samples, rownames(fcm_mat))
    if(length(fcm_samples)>=10) {
      for(gs in rownames(gsva_scores)) for(fcm_col in colnames(fcm_mat)) {
        valid_idx <- !is.na(gsva_scores[gs,fcm_samples]) & !is.na(fcm_mat[fcm_samples,fcm_col]) & 
          is.finite(gsva_scores[gs,fcm_samples]) & is.finite(fcm_mat[fcm_samples,fcm_col])
        if(sum(valid_idx)>=10) {
          sp_test <- cor.test(gsva_scores[gs,fcm_samples][valid_idx], fcm_mat[fcm_samples,fcm_col][valid_idx], 
                              method="spearman")
          external_cor <- rbind(external_cor, data.frame(GeneSet=gs, External_Marker=fcm_col, Type="FCM",
                                                         N=sum(valid_idx), Spearman_rho=sp_test$estimate, 
                                                         P_value=sp_test$p.value, stringsAsFactors=FALSE))
        }
      }
    }
  }
  if(ncol(cyto_mat)>0) {
    cyto_samples <- intersect(ra_samples, rownames(cyto_mat))
    if(length(cyto_samples)>=10) {
      for(gs in rownames(gsva_scores)) for(cyto_col in colnames(cyto_mat)) {
        valid_idx <- !is.na(gsva_scores[gs,cyto_samples]) & !is.na(cyto_mat[cyto_samples,cyto_col]) & 
          is.finite(gsva_scores[gs,cyto_samples]) & is.finite(cyto_mat[cyto_samples,cyto_col])
        if(sum(valid_idx)>=10) {
          sp_test <- cor.test(gsva_scores[gs,cyto_samples][valid_idx], cyto_mat[cyto_samples,cyto_col][valid_idx], 
                              method="spearman")
          external_cor <- rbind(external_cor, data.frame(GeneSet=gs, External_Marker=cyto_col, Type="Cytokine",
                                                         N=sum(valid_idx), Spearman_rho=sp_test$estimate, 
                                                         P_value=sp_test$p.value, stringsAsFactors=FALSE))
        }
      }
    }
  }
  if(nrow(external_cor)>0) {
    external_cor$P_adjusted <- p.adjust(external_cor$P_value, method="BH")
    external_cor <- external_cor %>% arrange(P_value)
    safe_write(external_cor, "results/tables/GSVA_external_correlations.csv", row.names=FALSE)
  }
  
  gsva_scores_df <- data.frame(GeneSet=rownames(gsva_scores), gsva_scores, 
                               check.names=FALSE, stringsAsFactors=FALSE)
  safe_write(gsva_scores_df, "results/tables/GSVA_all_scores.csv", row.names=FALSE)
  
  gsva_results <- list(scores=gsva_scores, score_df=gsva_scores_df, gene_sets=gene_sets_filtered, method=gsva_method,
                       infection_comparison=if(nrow(gsva_infection_comp)>0) gsva_infection_comp else NULL,
                       disease_comparison=if(nrow(gsva_disease_comp)>0) gsva_disease_comp else NULL,
                       external_correlations=if(nrow(external_cor)>0) external_cor else NULL)
  cat(sprintf("  %d gene sets, method=%s\n", length(gene_sets_filtered), gsva_method))
}, error=function(e) cat("  ✗ Skipped\n"))

# Microbiome
cat("\n=== Microbiome ===\n")
microbiome_samples <- grep("^(RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)
abundance_matrix <- as.matrix(microbiome_raw[, microbiome_samples])
rownames(abundance_matrix) <- microbiome_raw$Taxonomy
abundance_matrix[is.na(abundance_matrix)] <- 0
abundance_matrix <- abundance_matrix[rowSums(abundance_matrix)>0,]
diversity_df <- data.frame(Sample_ID=colnames(abundance_matrix),
                           Shannon=apply(abundance_matrix, 2, function(x) {x<-x[x>0]; if(length(x)==0) 0 else {p<-x/sum(x); -sum(p*log(p))}}),
                           Simpson=apply(abundance_matrix, 2, function(x) {x<-x[x>0]; if(length(x)==0) 0 else {p<-x/sum(x); 1-sum(p^2)}})) %>%
  left_join(master_data[,c("Sample_ID","Sample_Group")], by="Sample_ID")
safe_write(diversity_df, "results/tables/Alpha_diversity.csv", row.names=FALSE)
bray_dist <- vegdist(t(abundance_matrix)+1, method="bray")
set.seed(42)
permanova_res <- adonis2(bray_dist~Sample_Group, data=diversity_df %>% dplyr::filter(!is.na(Sample_Group)), permutations=999)
microbiome_results <- list(abundance=abundance_matrix, alpha_diversity=diversity_df, 
                           beta_diversity=bray_dist, permanova=permanova_res)
cat("✓ Complete\n")

# Cytokine/FCM
cat("\n=== Cytokine/FCM ===\n")
compare_groups <- function(data, cols, group_col, g1, g2) {
  results <- data.frame(Variable=cols, p_value=NA, log2FC=NA)
  for(i in seq_along(cols)) {
    v1 <- as.numeric(data[[cols[i]]][data[[group_col]]==g1]); v1 <- v1[!is.na(v1)]
    v2 <- as.numeric(data[[cols[i]]][data[[group_col]]==g2]); v2 <- v2[!is.na(v2)]
    if(length(v1)>=3 && length(v2)>=3) {
      results$p_value[i] <- tryCatch(wilcox.test(v1,v2, exact=TRUE)$p.value, error=function(e) NA)
      results$log2FC[i] <- log2((median(v2)+0.01)/(median(v1)+0.01))
    }
  }
  results$p_adjusted <- p.adjust(results$p_value, method="BH"); results
}
balf_cyto_stats <- compare_groups(master_data, balf_cytokine_cols, "Sample_Group", "Control", "RA")
safe_write(balf_cyto_stats, "results/tables/Cytokines_BALF.csv", row.names=FALSE)
serum_cyto_stats <- compare_groups(master_data, serum_cytokine_cols, "Sample_Group", "Control", "RA")
safe_write(serum_cyto_stats, "results/tables/Cytokines_Serum.csv", row.names=FALSE)
balf_fcm_stats <- compare_groups(master_data, balf_fcm_cols, "Sample_Group", "Control", "RA")
safe_write(balf_fcm_stats, "results/tables/FCM_BALF.csv", row.names=FALSE)
cytokine_fcm_results <- list(BALF_cytokines=balf_cyto_stats, Serum_cytokines=serum_cyto_stats, BALF_FCM=balf_fcm_stats)
cat("✓ Complete\n")

# Th17.1 Correlations
cat("\n=== Th17.1 ===\n")
th17_correlations <- NULL
tryCatch({
  th17_col <- grep("Th17\\.1|Th17_1", colnames(master_data), value=TRUE, ignore.case=TRUE)
  if(length(th17_col)>0) {
    ra_samples_th17 <- common_samples[meta$Sample_Group=="RA"]
    th17_vals <- as.numeric(master_data[[th17_col[1]]][match(ra_samples_th17, master_data$Sample_ID)])
    names(th17_vals) <- ra_samples_th17; th17_vals <- th17_vals[!is.na(th17_vals)]
    if(length(th17_vals)>=10) {
      th17_cor_results <- data.frame()
      common_samples_th17 <- intersect(names(th17_vals), colnames(log_cpm))
      if(length(common_samples_th17)>=10) {
        for(gene in rownames(log_cpm)[1:min(1000,nrow(log_cpm))]) {
          gene_vals <- log_cpm[gene, common_samples_th17]; th17_matched <- th17_vals[common_samples_th17]
          if(sum(!is.na(gene_vals) & !is.na(th17_matched))>=10) {
            sp_test <- cor.test(gene_vals, th17_matched, method="spearman", exact=TRUE)
            th17_cor_results <- rbind(th17_cor_results, data.frame(Gene=gene, Spearman_rho=sp_test$estimate,
                                                                   P_value=sp_test$p.value, stringsAsFactors=FALSE))
          }
        }
        if(nrow(th17_cor_results)>0) {
          th17_cor_results$P_adjusted <- p.adjust(th17_cor_results$P_value, method="BH")
          th17_cor_results <- th17_cor_results %>% arrange(P_value)
          safe_write(th17_cor_results, "results/tables/Th17_1_gene_correlations.csv", row.names=FALSE)
        }
      }
      th17_correlations <- th17_cor_results
    }
  }
  cat("✓ Complete\n")
}, error=function(e) cat("✗ Skipped\n"))

# Deconvolution
cat("\n=== Deconvolution ===\n")
deconv_integration_results <- NULL
if(deconv_loaded) {
  bp_samples <- rownames(bp_theta); deconv_common <- intersect(bp_samples, common_samples)
  if(length(deconv_common)>=10) {
    deconv_mat <- bp_theta[deconv_common,]*100; deconv_cell <- bp_cell[deconv_common,]*100
    cell_types_all <- colnames(bp_theta)
    deconv_unique_types <- setdiff(cell_types_all, c("Macrophage","T_cell","Neutrophil","Other","Epithelial"))
    
    # 3-group label (Total_Score based)
    deconv_subgroup <- factor(
      meta$Subgroup[match(deconv_common, meta$Sample_ID)],
      levels = c("Control", "RA_nonILD", "RA_ILD"))
    
    # 2-group label (RA vs Control)
    deconv_group <- factor(
      ifelse(grepl("^Sarcoidosis", deconv_common), "Control", "RA"),
      levels = c("Control", "RA"))
    
    deconv_export <- data.frame(Sample_ID=deconv_common, Group=as.character(deconv_group),
                                Subgroup=as.character(deconv_subgroup), deconv_mat[deconv_common,])
    safe_write(deconv_export, "results/tables/Deconvolution_proportions_all.csv", row.names=FALSE)
    
    # Comparison 1: RA vs Control
    cat("  Comparison 1: RA vs Control\n")
    ra_vs_ctrl_comparison <- data.frame()
    for(ct in deconv_unique_types) {
      if(ct %in% colnames(deconv_mat)) {
        vals_ra <- deconv_mat[deconv_group=="RA",ct]; vals_ctrl <- deconv_mat[deconv_group=="Control",ct]
        if(length(vals_ra)>=3 && length(vals_ctrl)>=3) {
          wt <- wilcox.test(vals_ra, vals_ctrl, exact=TRUE)
          ra_vs_ctrl_comparison <- rbind(ra_vs_ctrl_comparison, data.frame(
            CellType=ct, Median_Control=median(vals_ctrl,na.rm=TRUE), Median_RA=median(vals_ra,na.rm=TRUE),
            P_value=wt$p.value, stringsAsFactors=FALSE))
        }
      }
    }
    if(nrow(ra_vs_ctrl_comparison)>0) {
      ra_vs_ctrl_comparison$P_adjusted <- p.adjust(ra_vs_ctrl_comparison$P_value, method="BH")
      ra_vs_ctrl_comparison <- ra_vs_ctrl_comparison %>% arrange(P_value)
      safe_write(ra_vs_ctrl_comparison, "results/tables/Deconv_RA_vs_Control.csv", row.names=FALSE)
      cat(sprintf("    Sig (padj<0.05): %d cell types\n", sum(ra_vs_ctrl_comparison$P_adjusted<0.05)))
    }
    
    # Comparison 2: RA-ILD vs RA-nonILD
    cat("  Comparison 2: RA-ILD vs RA-nonILD\n")
    ild_vs_nonild_comparison <- data.frame()
    for(ct in deconv_unique_types) {
      if(ct %in% colnames(deconv_mat)) {
        vals_ild <- deconv_mat[deconv_subgroup=="RA_ILD",ct]; vals_nonild <- deconv_mat[deconv_subgroup=="RA_nonILD",ct]
        if(length(vals_ild)>=3 && length(vals_nonild)>=3) {
          wt <- wilcox.test(vals_ild, vals_nonild, exact=TRUE)
          ild_vs_nonild_comparison <- rbind(ild_vs_nonild_comparison, data.frame(
            CellType=ct, Median_RA_nonILD=median(vals_nonild,na.rm=TRUE), Median_RA_ILD=median(vals_ild,na.rm=TRUE),
            P_value=wt$p.value, stringsAsFactors=FALSE))
        }
      }
    }
    if(nrow(ild_vs_nonild_comparison)>0) {
      ild_vs_nonild_comparison$P_adjusted <- p.adjust(ild_vs_nonild_comparison$P_value, method="BH")
      ild_vs_nonild_comparison <- ild_vs_nonild_comparison %>% arrange(P_value)
      safe_write(ild_vs_nonild_comparison, "results/tables/Deconv_RA_ILD_vs_RA_nonILD.csv", row.names=FALSE)
      cat(sprintf("    Sig (padj<0.05): %d cell types\n", sum(ild_vs_nonild_comparison$P_adjusted<0.05)))
    }
    
    # Comparison 3: RA-ILD vs Control
    cat("  Comparison 3: RA-ILD vs Control\n")
    ild_vs_ctrl_comparison <- data.frame()
    for(ct in deconv_unique_types) {
      if(ct %in% colnames(deconv_mat)) {
        vals_ild <- deconv_mat[deconv_subgroup=="RA_ILD",ct]; vals_ctrl <- deconv_mat[deconv_subgroup=="Control",ct]
        if(length(vals_ild)>=3 && length(vals_ctrl)>=3) {
          wt <- wilcox.test(vals_ild, vals_ctrl, exact=TRUE)
          ild_vs_ctrl_comparison <- rbind(ild_vs_ctrl_comparison, data.frame(
            CellType=ct, Median_Control=median(vals_ctrl,na.rm=TRUE), Median_RA_ILD=median(vals_ild,na.rm=TRUE),
            P_value=wt$p.value, stringsAsFactors=FALSE))
        }
      }
    }
    if(nrow(ild_vs_ctrl_comparison)>0) {
      ild_vs_ctrl_comparison$P_adjusted <- p.adjust(ild_vs_ctrl_comparison$P_value, method="BH")
      ild_vs_ctrl_comparison <- ild_vs_ctrl_comparison %>% arrange(P_value)
      safe_write(ild_vs_ctrl_comparison, "results/tables/Deconv_RA_ILD_vs_Control.csv", row.names=FALSE)
      cat(sprintf("    Sig (padj<0.05): %d cell types\n", sum(ild_vs_ctrl_comparison$P_adjusted<0.05)))
    }
    
    deconv_integration_results <- list(
      deconv_theta=deconv_mat, deconv_cell_fraction=deconv_cell,
      deconv_samples=deconv_common, deconv_group=deconv_group, deconv_subgroup=deconv_subgroup,
      cell_types=cell_types_all, unique_cell_types=deconv_unique_types,
      ra_vs_ctrl_comparison=if(nrow(ra_vs_ctrl_comparison)>0) ra_vs_ctrl_comparison else NULL,
      ild_vs_nonild_comparison=if(nrow(ild_vs_nonild_comparison)>0) ild_vs_nonild_comparison else NULL,
      ild_vs_ctrl_comparison=if(nrow(ild_vs_ctrl_comparison)>0) ild_vs_ctrl_comparison else NULL)
    save(master_data, gene_annotations, count_matrix, log_cpm, expr_matrix, expr_mat, cyto_mat, fcm_mat, meta,
         common_samples, deconv_mat, deconv_cell, deconv_common, deconv_group, deconv_subgroup,
         file="results/RA_ILD_Core_Data.RData")
    cat(sprintf("  Total: %d samples (Ctrl:%d, RA-nonILD:%d, RA-ILD:%d), %d types\n",
                length(deconv_common), sum(deconv_subgroup=="Control"),
                sum(deconv_subgroup=="RA_nonILD"), sum(deconv_subgroup=="RA_ILD"), length(cell_types_all)))
  }
} else cat("  Data not loaded\n")

# Surfactant
cat("\n=== Surfactant ===\n")
surfactant_analysis <- NULL
tryCatch({
  if("Infection_Group" %in% colnames(meta)) {
    surfactant_genes <- c("SFTPA1","SFTPA2","SFTPB","SFTPC","SFTPD","SFTA2","SFTA3","ABCA3","NKX2-1",
                          "NAPSA","LPCAT1","PGC","LAMP3","SLC34A2","MUC1","CLDN18","AQP5")
    available_surf <- intersect(surfactant_genes, rownames(log_cpm))
    if(length(available_surf)>=5) {
      ra_samples_surf <- common_samples[meta$Sample_Group=="RA"]; ra_meta_surf <- meta[ra_samples_surf,]
      inf_neg_surf <- ra_samples_surf[!is.na(ra_meta_surf$Infection_Group) & 
                                        ra_meta_surf$Infection_Group=="Infection_Negative"]
      inf_pos_surf <- ra_samples_surf[!is.na(ra_meta_surf$Infection_Group) & 
                                        ra_meta_surf$Infection_Group=="Infection_Positive"]
      if(length(inf_neg_surf)>=3 && length(inf_pos_surf)>=3) {
        surf_comparison <- data.frame()
        for(gene in available_surf) {
          vals_neg <- log_cpm[gene, inf_neg_surf]; vals_pos <- log_cpm[gene, inf_pos_surf]
          wt <- wilcox.test(vals_neg, vals_pos, exact=TRUE)
          surf_comparison <- rbind(surf_comparison, data.frame(Gene=gene, Mean_Inf_Neg=mean(vals_neg,na.rm=TRUE),
                                                               Mean_Inf_Pos=mean(vals_pos,na.rm=TRUE),
                                                               Log2FC=log2((mean(vals_pos,na.rm=TRUE)+0.1)/
                                                                             (mean(vals_neg,na.rm=TRUE)+0.1)),
                                                               P_value=wt$p.value, stringsAsFactors=FALSE))
        }
        surf_comparison$P_adjusted <- p.adjust(surf_comparison$P_value, method="BH")
        surf_comparison <- surf_comparison %>% arrange(P_value)
        safe_write(surf_comparison, "results/tables/Surfactant_genes_infection.csv", row.names=FALSE)
        surfactant_analysis <- list(genes=available_surf, comparison=surf_comparison)
        cat(sprintf("  %d genes\n", nrow(surf_comparison)))
      }
    }
  }
}, error=function(e) cat("✗ Skipped\n"))

cat("\n=== Saving ===\n")
results_list <- list(
  master_data=master_data, gene_annotations=gene_annotations, count_matrix=count_matrix, log_cpm=log_cpm,
  expr_matrix=expr_matrix, expr_mat=expr_mat, cyto_mat=cyto_mat, fcm_mat=fcm_mat, meta=meta,
  common_samples=common_samples, DEG_RA_vs_Control=res_df,
  GSVA=gsva_results, Microbiome=microbiome_results,
  Cytokines_FCM=cytokine_fcm_results, Th17_Correlations=th17_correlations,
  Deconvolution=deconv_integration_results, Surfactant_Analysis=surfactant_analysis, 
  analysis_date=Sys.time(), R_version=R.version.string
)
save(results_list, file="results/RA_ILD_Multiomics_Results.RData")
save.image(file="results/RA_ILD_Workspace.RData")
cat("✓ Saved\n")

# Summary
cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║              ANALYSIS COMPLETE                         ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Samples: %d (RA: %d [ILD: %d, nonILD: %d], Control: %d)\n", nrow(master_data), 
            sum(master_data$Sample_Group=="RA"),
            sum(master_data$Subgroup=="RA_ILD"), sum(master_data$Subgroup=="RA_nonILD"),
            sum(master_data$Sample_Group=="Control")))
cat(sprintf("Genes: %d\n", nrow(count_matrix)))
cat(sprintf("DEGs: %d\n", sum(res_df$Regulation!="NS", na.rm=TRUE)))
if(deconv_loaded && !is.null(deconv_integration_results)) 
  cat(sprintf("Deconvolution: %d types\n", length(deconv_integration_results$cell_types)))
cat("\n═══════════════════════════════════════════════════════════\n")