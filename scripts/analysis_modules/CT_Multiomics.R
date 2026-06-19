# ============================================================================
# CT Quantitative Integration + Multi-omics Analysis
# ============================================================================
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  CT Quantitative + Multi-omics Integration                   ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)


suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(DESeq2); library(edgeR)
  library(pROC); library(randomForest); library(caret)
  library(GSVA); library(pheatmap); library(effsize); library(glmnet)
})
set.seed(42)

out_dir <- "results"

# ============================================================================
# PART 1: CT Data — All Timepoints + Delta
# ============================================================================
cat("=== PART 1: CT Data Integration ===\n")

# Sample-ID to PatientID mapping table (rename the released data file to match)
mapping <- read_excel("sample_id_mapping.xlsx", sheet="Sheet1")
colnames(mapping)[1] <- "row_num"
id_map <- data.frame(
  Sample_ID = as.character(mapping[["No."]]),
  PatientID = as.character(mapping[["ID"]]),
  stringsAsFactors = FALSE
)
id_map <- id_map[grepl("^(KY|RA[0-9]|Sarcoidosis)", id_map$Sample_ID), ]

ct_quant <- read_excel("ct_lung_analysis_merged.xlsx", sheet="Merged_Data")
ct_quant$PatientID <- as.character(ct_quant$PatientID)
ct_quant$StudyDate <- as.character(ct_quant$StudyDate)

ct_all <- ct_quant %>% inner_join(id_map, by="PatientID") %>%
  filter(!Sample_ID %in% character(0))  # QC-excluded samples are absent from the released n=35 dataset

ct_vars <- c("total_lung_cm3","Air_pct","Healthy_Lung_pct","GGO_pct",
             "Consolidation_pct","Dense_Tissue_pct","lung_involvement_pct",
             "lung_mean_HU","lung_std_HU")

# T1 = earliest, T2 = latest per patient
ct_t1 <- ct_all %>% group_by(Sample_ID) %>%
  arrange(StudyDate) %>% dplyr::slice(1) %>% ungroup() %>%
  dplyr::select(Sample_ID, StudyDate, all_of(ct_vars))
colnames(ct_t1) <- c("Sample_ID","CT_T1_Date", paste0("CT_T1_", ct_vars))

ct_t2 <- ct_all %>% group_by(Sample_ID) %>%
  arrange(desc(StudyDate)) %>% dplyr::slice(1) %>% ungroup() %>%
  dplyr::select(Sample_ID, StudyDate, all_of(ct_vars))
colnames(ct_t2) <- c("Sample_ID","CT_T2_Date", paste0("CT_T2_", ct_vars))

# Patients with 2+ timepoints (T1 != T2)
two_tp <- ct_t1 %>% inner_join(ct_t2[, c("Sample_ID","CT_T2_Date")], by="Sample_ID") %>%
  filter(CT_T1_Date != CT_T2_Date)
two_tp_ids <- two_tp$Sample_ID
cat(sprintf("  Two-timepoint data: %d patients\n", length(two_tp_ids)))

# Delta = T2 - T1
delta_df <- data.frame(Sample_ID = two_tp_ids)
for(v in ct_vars) {
  t1_vals <- as.numeric(ct_t1[[paste0("CT_T1_", v)]][match(two_tp_ids, ct_t1$Sample_ID)])
  t2_vals <- as.numeric(ct_t2[[paste0("CT_T2_", v)]][match(two_tp_ids, ct_t2$Sample_ID)])
  delta_df[[paste0("CT_Delta_", v)]] <- t2_vals - t1_vals
}
# Days between
d1 <- as.Date(ct_t1$CT_T1_Date[match(two_tp_ids, ct_t1$Sample_ID)], format="%Y%m%d")
d2 <- as.Date(ct_t2$CT_T2_Date[match(two_tp_ids, ct_t2$Sample_ID)], format="%Y%m%d")
delta_df$CT_Delta_Days <- as.numeric(d2 - d1)

# Annualized rate of change
for(v in ct_vars) {
  delta_df[[paste0("CT_Rate_", v)]] <- delta_df[[paste0("CT_Delta_", v)]] / delta_df$CT_Delta_Days * 365.25
}

cat(sprintf("  Delta computed: %d patients, %d days mean (range %d-%d)\n",
    nrow(delta_df), round(mean(delta_df$CT_Delta_Days)),
    min(delta_df$CT_Delta_Days), max(delta_df$CT_Delta_Days)))

# Merge all CT data (1 row per patient)
ct_integrated <- ct_t1 %>%
  left_join(ct_t2, by="Sample_ID") %>%
  left_join(delta_df, by="Sample_ID")

write.csv(ct_integrated, file.path(out_dir, "tables/CT_quantitative_all_timepoints.csv"), row.names=FALSE)

# ============================================================================
# PART 2: master_data + CT integration
# ============================================================================
cat("\n=== PART 2: Master Data Integration ===\n")

master_ext <- master_data %>% left_join(ct_integrated, by="Sample_ID")
cat(sprintf("  master_data: %d rows x %d cols\n", nrow(master_ext), ncol(master_ext)))
cat(sprintf("  CT T1 present: %d, T2 present: %d, Delta present: %d\n",
    sum(!is.na(master_ext$CT_T1_GGO_pct)),
    sum(!is.na(master_ext$CT_T2_GGO_pct)),
    sum(!is.na(master_ext$CT_Delta_GGO_pct))))

# ============================================================================
# PART 3: CT parameters — Group comparisons with effect sizes
# ============================================================================
cat("\n=== PART 3: CT Group Comparisons ===\n")

ct_comparison <- data.frame()
ct_cols_all <- grep("^CT_T1_|^CT_T2_|^CT_Delta_|^CT_Rate_", colnames(master_ext), value=TRUE)
ct_cols_all <- ct_cols_all[!grepl("Date|Days", ct_cols_all)]

for(col in ct_cols_all) {
  for(comp in list(c("RA","Control","Sample_Group"),
                    c("RA_ILD","RA_nonILD","Subgroup"),
                    c("RA_ILD","Control","Subgroup"))) {
    g1_vals <- as.numeric(master_ext[[col]][master_ext[[comp[3]]] == comp[1]])
    g2_vals <- as.numeric(master_ext[[col]][master_ext[[comp[3]]] == comp[2]])
    g1_vals <- g1_vals[!is.na(g1_vals)]; g2_vals <- g2_vals[!is.na(g2_vals)]
    if(length(g1_vals) >= 3 && length(g2_vals) >= 3) {
      wt <- wilcox.test(g1_vals, g2_vals, exact=TRUE)
      cd <- tryCatch(cliff.delta(g1_vals, g2_vals), error=function(e) list(estimate=NA, magnitude="NA"))
      ct_comparison <- rbind(ct_comparison, data.frame(
        Variable=col, Comparison=paste(comp[1:2], collapse=" vs "),
        Mean_G1=mean(g1_vals), Mean_G2=mean(g2_vals),
        N_G1=length(g1_vals), N_G2=length(g2_vals),
        P_value=wt$p.value, Cliff_delta=cd$estimate, Magnitude=cd$magnitude
      ))
    }
  }
}

# BH correction within comparison
for(comp_name in unique(ct_comparison$Comparison)) {
  idx <- ct_comparison$Comparison == comp_name
  ct_comparison$P_adjusted[idx] <- p.adjust(ct_comparison$P_value[idx], method="BH")
}
ct_comparison <- ct_comparison %>% arrange(P_value)
write.csv(ct_comparison, file.path(out_dir, "tables/CT_group_comparisons.csv"), row.names=FALSE)

cat("  Significant (p<0.05 unadjusted):\n")
sig_ct <- ct_comparison[ct_comparison$P_value < 0.05, ]
if(nrow(sig_ct) > 0) {
  for(i in 1:min(15, nrow(sig_ct)))
    cat(sprintf("    %-35s %-20s p=%.4f delta=%.3f\n",
        sig_ct$Variable[i], sig_ct$Comparison[i], sig_ct$P_value[i], sig_ct$Cliff_delta[i]))
}

# ============================================================================
# PART 4: CT vs ILD_Score correlation (RA patients)
# ============================================================================
cat("\n=== PART 4: CT vs ILD_Score Correlation ===\n")

ra_ext <- master_ext[master_ext$Sample_Group == "RA", ]
ct_ild_cor <- data.frame()

for(col in ct_cols_all) {
  vals <- as.numeric(ra_ext[[col]])
  ild <- as.numeric(ra_ext$ILD_Score)
  valid <- !is.na(vals) & !is.na(ild)
  if(sum(valid) >= 8) {
    sp <- cor.test(vals[valid], ild[valid], method="spearman", exact=TRUE)
    ct_ild_cor <- rbind(ct_ild_cor, data.frame(
      Variable=col, Spearman_rho=sp$estimate, P_value=sp$p.value, N=sum(valid)))
  }
}
ct_ild_cor$P_adjusted <- p.adjust(ct_ild_cor$P_value, method="BH")
ct_ild_cor <- ct_ild_cor %>% arrange(P_value)
write.csv(ct_ild_cor, file.path(out_dir, "tables/CT_vs_ILD_Score_correlation.csv"), row.names=FALSE)

cat("  Top correlations:\n")
for(i in 1:min(10, nrow(ct_ild_cor)))
  cat(sprintf("    %-35s rho=%+.3f p=%.4f\n",
      ct_ild_cor$Variable[i], ct_ild_cor$Spearman_rho[i], ct_ild_cor$P_value[i]))

# ============================================================================
# PART 5: CT vs Deconvolution / GSVA / Cytokines
# ============================================================================
cat("\n=== PART 5: CT vs Multi-omics Correlations (RA only) ===\n")

ct_multiomics_cor <- data.frame()
# Filter to RA only — sarcoidosis excluded from CT progression analysis
ra_ct_ext <- master_ext[master_ext$Sample_Group == "RA", ]
cat(sprintf("  RA patients for CT correlations: %d\n", nrow(ra_ct_ext)))

# Use T1 CT variables + Delta/Rate for available patients
ct_key_vars <- c("CT_T1_GGO_pct","CT_T1_Consolidation_pct","CT_T1_Dense_Tissue_pct",
                  "CT_T1_lung_involvement_pct","CT_T1_Healthy_Lung_pct","CT_T1_lung_mean_HU",
                  "CT_Delta_GGO_pct","CT_Delta_lung_involvement_pct","CT_Delta_Healthy_Lung_pct",
                  "CT_Rate_GGO_pct","CT_Rate_lung_involvement_pct")
ct_key_vars <- intersect(ct_key_vars, colnames(ra_ct_ext))

# vs Deconvolution
if(exists("deconv_mat")) {
  deconv_samples_ct <- intersect(rownames(deconv_mat), ra_ct_ext$Sample_ID[!is.na(ra_ct_ext$CT_T1_GGO_pct)])
  for(ct_v in ct_key_vars) {
    for(cell_t in colnames(deconv_mat)) {
      ct_vals <- as.numeric(ra_ct_ext[[ct_v]][match(deconv_samples_ct, ra_ct_ext$Sample_ID)])
      dc_vals <- deconv_mat[deconv_samples_ct, cell_t]  # already in %
      valid <- !is.na(ct_vals) & !is.na(dc_vals)
      if(sum(valid) >= 8) {
        sp <- cor.test(ct_vals[valid], dc_vals[valid], method="spearman", exact=TRUE)
        ct_multiomics_cor <- rbind(ct_multiomics_cor, data.frame(
          CT_Variable=ct_v, Omics_Variable=cell_t, Category="Deconvolution",
          Spearman_rho=sp$estimate, P_value=sp$p.value, N=sum(valid)))
      }
    }
  }
}

# vs GSVA
if(exists("gsva_scores")) {
  gsva_samples_ct <- intersect(colnames(gsva_scores), ra_ct_ext$Sample_ID[!is.na(ra_ct_ext$CT_T1_GGO_pct)])
  for(ct_v in ct_key_vars) {
    for(gs in rownames(gsva_scores)) {
      ct_vals <- as.numeric(ra_ct_ext[[ct_v]][match(gsva_samples_ct, ra_ct_ext$Sample_ID)])
      gs_vals <- gsva_scores[gs, gsva_samples_ct]
      valid <- !is.na(ct_vals) & !is.na(gs_vals)
      if(sum(valid) >= 8) {
        sp <- cor.test(ct_vals[valid], gs_vals[valid], method="spearman", exact=TRUE)
        ct_multiomics_cor <- rbind(ct_multiomics_cor, data.frame(
          CT_Variable=ct_v, Omics_Variable=gs, Category="GSVA_Pathway",
          Spearman_rho=sp$estimate, P_value=sp$p.value, N=sum(valid)))
      }
    }
  }
}

# vs Cytokines (key serum)
cyto_key <- grep("^BALF_|^Serum_|^PB_", colnames(ra_ct_ext), value=TRUE)
cyto_key <- cyto_key[!grepl("Treg|Th|Macro|Neutro|Lympho|CD|Percent|Activated", cyto_key)]
# Also add FACS variables
facs_key <- grep("^BALF_.*Percent|^BALF_.*Macrophage|^BALF_.*Activated|^PB_.*Th|^PB_.*Treg|^PB_.*Neutro|^BALF_.*Th", colnames(ra_ct_ext), value=TRUE)
all_biomarkers <- unique(c(cyto_key, facs_key))
for(ct_v in ct_key_vars) {
  for(cy in all_biomarkers) {
    ct_vals <- as.numeric(ra_ct_ext[[ct_v]])
    cy_vals <- as.numeric(ra_ct_ext[[cy]])
    valid <- !is.na(ct_vals) & !is.na(cy_vals)
    if(sum(valid) >= 8) {
      sp <- cor.test(ct_vals[valid], cy_vals[valid], method="spearman", exact=TRUE)
      ct_multiomics_cor <- rbind(ct_multiomics_cor, data.frame(
        CT_Variable=ct_v, Omics_Variable=cy, Category="FACS/Cytokine",
        Spearman_rho=sp$estimate, P_value=sp$p.value, N=sum(valid)))
    }
  }
}

ct_multiomics_cor$P_adjusted <- p.adjust(ct_multiomics_cor$P_value, method="BH")
ct_multiomics_cor <- ct_multiomics_cor %>% arrange(P_value)
write.csv(ct_multiomics_cor, file.path(out_dir, "tables/CT_vs_Multiomics_correlations.csv"), row.names=FALSE)

cat(sprintf("  Total correlations tested: %d\n", nrow(ct_multiomics_cor)))
cat(sprintf("  Significant (padj<0.05): %d\n", sum(ct_multiomics_cor$P_adjusted < 0.05)))
cat("\n  Top 15 correlations:\n")
for(i in 1:min(15, nrow(ct_multiomics_cor)))
  cat(sprintf("    %-30s × %-25s rho=%+.3f p=%.4f padj=%.4f\n",
      ct_multiomics_cor$CT_Variable[i], ct_multiomics_cor$Omics_Variable[i],
      ct_multiomics_cor$Spearman_rho[i], ct_multiomics_cor$P_value[i], ct_multiomics_cor$P_adjusted[i]))

# ============================================================================
# PART 6: Elastic Net — ILD_Score prediction with CT + omics
# ============================================================================
cat("\n=== PART 6: Elastic Net with CT + Omics ===\n")

ra_samples_en <- intersect(common_samples[meta$Sample_Group=="RA"],
                            master_ext$Sample_ID[!is.na(master_ext$CT_T1_GGO_pct)])
cat(sprintf("  RA samples with CT + omics: %d\n", length(ra_samples_en)))

if(length(ra_samples_en) >= 12) {
  # Feature matrix: CT T1 + Delta + top expression + cytokines + FACS + deconvolution + GSVA
  feature_list <- list()

  # CT T1
  ct_t1_cols <- grep("^CT_T1_", colnames(master_ext), value=TRUE)
  ct_t1_cols <- ct_t1_cols[!grepl("Date", ct_t1_cols)]
  ct_t1_data <- master_ext[match(ra_samples_en, master_ext$Sample_ID), ct_t1_cols]
  for(col in colnames(ct_t1_data)) ct_t1_data[[col]] <- as.numeric(ct_t1_data[[col]])
  feature_list$CT_T1 <- ct_t1_data

  # CT Delta (if available)
  ct_delta_cols <- grep("^CT_Delta_", colnames(master_ext), value=TRUE)
  ct_delta_cols <- ct_delta_cols[!grepl("Days", ct_delta_cols)]
  if(length(ct_delta_cols) > 0) {
    ct_d_data <- master_ext[match(ra_samples_en, master_ext$Sample_ID), ct_delta_cols]
    for(col in colnames(ct_d_data)) ct_d_data[[col]] <- as.numeric(ct_d_data[[col]])
    feature_list$CT_Delta <- ct_d_data
  }

  # Expression (top 50)
  expr_en <- expr_mat[ra_samples_en, 1:min(50, ncol(expr_mat))]
  feature_list$Expression <- as.data.frame(expr_en)

  # Cytokines
  cyto_en <- cyto_mat[ra_samples_en, ]
  feature_list$Cytokines <- as.data.frame(cyto_en)

  # FACS
  facs_en <- facs_mat[ra_samples_en, ]
  feature_list$FACS <- as.data.frame(facs_en)

  # Deconvolution
  if(exists("deconv_mat") && all(ra_samples_en %in% rownames(deconv_mat))) {
    deconv_en <- as.data.frame(deconv_mat[ra_samples_en, ])  # already in %
    colnames(deconv_en) <- paste0("Deconv_", colnames(deconv_en))
    feature_list$Deconvolution <- deconv_en
  }

  # GSVA
  if(exists("gsva_scores") && all(ra_samples_en %in% colnames(gsva_scores))) {
    gsva_en <- as.data.frame(t(gsva_scores[, ra_samples_en]))
    colnames(gsva_en) <- paste0("GSVA_", colnames(gsva_en))
    feature_list$GSVA <- gsva_en
  }

  # Combine
  X_all <- do.call(cbind, feature_list)
  X_all[is.na(X_all)] <- 0
  X_all[!is.finite(as.matrix(X_all))] <- 0

  # Remove zero-variance
  fv <- apply(X_all, 2, var, na.rm=TRUE)
  X_all <- X_all[, !is.na(fv) & fv > 0]

  y <- as.numeric(master_ext$ILD_Score[match(ra_samples_en, master_ext$Sample_ID)])

  cat(sprintf("  Features: %d (CT:%d, Expr:%d, Cyto:%d, FACS:%d, Deconv:%d, GSVA:%d)\n",
      ncol(X_all),
      sum(grepl("^CT_", colnames(X_all))),
      sum(grepl("^[A-Z].*\\.", colnames(X_all)) & !grepl("^CT_|^BALF_|^Serum_|^Deconv_|^GSVA_", colnames(X_all))),
      sum(grepl("^BALF_|^Serum_", colnames(X_all)) & !grepl("Treg|Th|Macro|Neutro", colnames(X_all))),
      sum(grepl("^BALF_.*Treg|^BALF_.*Th|^PB_", colnames(X_all))),
      sum(grepl("^Deconv_", colnames(X_all))),
      sum(grepl("^GSVA_", colnames(X_all)))))

  # Elastic Net (alpha=0.5)
  X_mat <- as.matrix(X_all)
  cv_fit <- cv.glmnet(X_mat, y, alpha=0.5, nfolds=min(10, length(y)))
  coefs <- coef(cv_fit, s="lambda.min")
  nonzero <- which(coefs[-1, 1] != 0)
  en_features <- data.frame(
    Feature = rownames(coefs)[-1][nonzero],
    Coefficient = coefs[-1, 1][nonzero]
  ) %>% arrange(desc(abs(Coefficient)))

  # Feature type
  en_features$Type <- case_when(
    grepl("^CT_", en_features$Feature) ~ "CT",
    grepl("^Deconv_", en_features$Feature) ~ "Deconvolution",
    grepl("^GSVA_", en_features$Feature) ~ "GSVA",
    grepl("^BALF_|^Serum_", en_features$Feature) ~ "Cytokine/FACS",
    TRUE ~ "Expression"
  )

  write.csv(en_features, file.path(out_dir, "tables/ElasticNet_CT_Multiomics_features.csv"), row.names=FALSE)

  cat(sprintf("  Selected features: %d\n", nrow(en_features)))
  cat("  Feature type breakdown:\n")
  print(table(en_features$Type))

  cat("\n  Top 15 features:\n")
  for(i in 1:min(15, nrow(en_features)))
    cat(sprintf("    %-40s %+8.4f [%s]\n", en_features$Feature[i], en_features$Coefficient[i], en_features$Type[i]))

  # LOOCV
  loocv_pred <- numeric(length(y))
  for(i in 1:length(y)) {
    cv_i <- cv.glmnet(X_mat[-i, ], y[-i], alpha=0.5, nfolds=min(10, length(y)-1))
    loocv_pred[i] <- predict(cv_i, newx=X_mat[i, , drop=FALSE], s="lambda.min")
  }
  loocv_cor <- cor.test(y, loocv_pred, method="spearman", exact=TRUE)
  loocv_r2 <- 1 - sum((y - loocv_pred)^2) / sum((y - mean(y))^2)

  cat(sprintf("\n  LOOCV: rho=%.3f (p=%.4f), R2=%.3f\n",
      loocv_cor$estimate, loocv_cor$p.value, loocv_r2))

  report_en <- list(n_features=nrow(en_features), features=en_features,
                     loocv_rho=loocv_cor$estimate, loocv_p=loocv_cor$p.value, loocv_r2=loocv_r2)
}

# ============================================================================
# PART 7: CT Progression Heatmap
# ============================================================================
cat("\n=== PART 7: CT Progression Heatmap ===\n")

tryCatch({
  # Delta data for visualization
  delta_viz <- delta_df %>%
    left_join(master_ext[, c("Sample_ID","Sample_Group","Subgroup")], by="Sample_ID")

  delta_mat <- delta_viz %>%
    dplyr::select(Sample_ID, CT_Delta_GGO_pct, CT_Delta_Consolidation_pct, CT_Delta_Dense_Tissue_pct,
           CT_Delta_lung_involvement_pct, CT_Delta_Healthy_Lung_pct) %>%
    column_to_rownames("Sample_ID") %>% as.matrix()
  colnames(delta_mat) <- gsub("CT_Delta_", "", colnames(delta_mat))

  # Annotation
  ann <- data.frame(
    Group = delta_viz$Subgroup[match(rownames(delta_mat), delta_viz$Sample_ID)],
    row.names = rownames(delta_mat)
  )
  ann_colors <- list(Group = c("Control"="#4DAF4A", "RA_nonILD"="#377EB8", "RA_ILD"="#E41A1C"))

  pdf(file.path(out_dir, "figures/CT_Delta_heatmap.pdf"), width=10, height=12)
  pheatmap(delta_mat, cluster_rows=TRUE, cluster_cols=TRUE,
           annotation_row=ann, annotation_colors=ann_colors,
           main="CT Change (T2 - T1)", fontsize_row=8,
           color=colorRampPalette(c("blue","white","red"))(100),
           breaks=seq(-20, 20, length.out=101))
  dev.off()
  cat("  ✓ Heatmap saved\n")
}, error=function(e) { tryCatch(dev.off(), error=function(x) NULL); cat(sprintf("  ✗ %s\n", e$message)) })

# ============================================================================
# PART 8: Summary
# ============================================================================
cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  Complete                                                      ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

save(master_ext, ct_integrated, delta_df, ct_all, ct_comparison,
     ct_ild_cor, ct_multiomics_cor,
     file=file.path(out_dir, "CT_Multiomics_Results.RData"))

cat("Output files:\n")
for(f in list.files(file.path(out_dir, "tables"), pattern="CT_"))
  cat(sprintf("  ✓ tables/%s\n", f))
for(f in list.files(file.path(out_dir, "figures"), pattern="CT_"))
  cat(sprintf("  ✓ figures/%s\n", f))
