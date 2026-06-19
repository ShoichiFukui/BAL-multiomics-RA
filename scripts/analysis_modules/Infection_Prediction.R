# ============================================================================
# Outcome 2: Future Respiratory Infection Prediction
# Multi-omics integrated analysis
# ============================================================================
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)


suppressPackageStartupMessages({
  library(tidyverse); library(DESeq2); library(edgeR)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(pROC); library(randomForest); library(caret)
  library(effsize); library(glmnet); library(GSVA)
})
set.seed(42)

out_dir <- "results"
for(d in c(file.path(out_dir,"tables"), file.path(out_dir,"figures")))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  Future Infection Prediction — Multi-omics Analysis          ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Setup: RA patients with infection outcome
# ============================================================================
ra_all <- master_data[master_data$Sample_Group == "RA", ]
ra_inf <- ra_all[!is.na(ra_all$respiratory_infection), ]
ra_inf$Infection <- factor(ifelse(ra_inf$respiratory_infection==1, "Inf_pos", "Inf_neg"),
                            levels=c("Inf_neg","Inf_pos"))
ra_common <- intersect(ra_inf$Sample_ID, common_samples)
ra_inf_c <- ra_inf[match(ra_common, ra_inf$Sample_ID), ]

cat(sprintf("RA with infection outcome + multi-omics: %d\n", length(ra_common)))
cat(sprintf("  Future infection(+): %d, infection(-): %d\n",
    sum(ra_inf_c$Infection=="Inf_pos"), sum(ra_inf_c$Infection=="Inf_neg")))

# ============================================================================
# SECTION 1: Comprehensive group comparison with effect sizes
# ============================================================================
cat("\n=== SECTION 1: Infection(+) vs (-) — All Omics Layers ===\n")

inf_comparison <- data.frame()

compare_inf <- function(vals_pos, vals_neg, variable, category) {
  if(length(vals_pos) < 3 || length(vals_neg) < 3) return(NULL)
  wt <- tryCatch(wilcox.test(vals_pos, vals_neg, exact=TRUE), error=function(e) NULL)
  cd <- tryCatch(cliff.delta(vals_pos, vals_neg), error=function(e) list(estimate=NA, magnitude="NA"))
  if(is.null(wt)) return(NULL)
  data.frame(Variable=variable, Category=category,
             Median_InfPos=median(vals_pos, na.rm=TRUE),
             Median_InfNeg=median(vals_neg, na.rm=TRUE),
             N_pos=length(vals_pos), N_neg=length(vals_neg),
             P_value=wt$p.value, Cliff_delta=cd$estimate, Magnitude=cd$magnitude)
}

# 1a. Cytokines (BALF + Serum) — prefix retained to distinguish compartments
cyto_cols <- grep("^BALF_|^Serum_", colnames(ra_inf_c), value=TRUE)
cyto_cols <- cyto_cols[!grepl("Treg|Th1|Th2|Th17|Macro|Neutro|Lympho|Eosin|Baso|Plasma_Cell|CD[0-9]|Percent|Activated|Ratio|CD14|CD86|CD66", cyto_cols)]
cat(sprintf("  Cytokine columns: BALF=%d, Serum=%d\n",
    sum(grepl("^BALF_",cyto_cols)), sum(grepl("^Serum_",cyto_cols))))
for(col in cyto_cols) {
  vals <- as.numeric(ra_inf_c[[col]])
  pos <- vals[ra_inf_c$Infection=="Inf_pos"]; neg <- vals[ra_inf_c$Infection=="Inf_neg"]
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  category <- ifelse(grepl("^Serum_",col), "Serum_Cytokine", "BALF_Cytokine")
  res <- compare_inf(pos, neg, col, category)
  if(!is.null(res)) inf_comparison <- rbind(inf_comparison, res)
}

# 1b. FACS
facs_cols_inf <- grep("^BALF_.*Treg|^BALF_.*Th|^BALF_.*Macro|^BALF_.*Neutro|^BALF_.*Activated|^PB_|CD[0-9]|Ratio|Percent",
                       colnames(ra_inf_c), value=TRUE)
facs_cols_inf <- facs_cols_inf[!grepl("^CT_|^BALF_EGF|^BALF_Eotaxin", facs_cols_inf)]
for(col in facs_cols_inf) {
  vals <- as.numeric(ra_inf_c[[col]])
  pos <- vals[ra_inf_c$Infection=="Inf_pos"]; neg <- vals[ra_inf_c$Infection=="Inf_neg"]
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  res <- compare_inf(pos, neg, col, "FCM")
  if(!is.null(res)) inf_comparison <- rbind(inf_comparison, res)
}

# 1c. Deconvolution
if(exists("deconv_mat")) {
  deconv_ids <- intersect(ra_common, rownames(deconv_mat))
  inf_deconv <- ra_inf_c$Infection[match(deconv_ids, ra_inf_c$Sample_ID)]
  for(ct in colnames(deconv_mat)) {
    # deconv_mat is already *100 (percentage)
    pos <- deconv_mat[deconv_ids[inf_deconv=="Inf_pos"], ct]
    neg <- deconv_mat[deconv_ids[inf_deconv=="Inf_neg"], ct]
    res <- compare_inf(pos, neg, ct, "Deconvolution")
    if(!is.null(res)) inf_comparison <- rbind(inf_comparison, res)
  }
}

# 1d. GSVA
if(exists("gsva_scores")) {
  gsva_ids <- intersect(ra_common, colnames(gsva_scores))
  inf_gsva <- ra_inf_c$Infection[match(gsva_ids, ra_inf_c$Sample_ID)]
  for(gs in rownames(gsva_scores)) {
    pos <- gsva_scores[gs, gsva_ids[inf_gsva=="Inf_pos"]]
    neg <- gsva_scores[gs, gsva_ids[inf_gsva=="Inf_neg"]]
    res <- compare_inf(pos, neg, gs, "GSVA_Pathway")
    if(!is.null(res)) inf_comparison <- rbind(inf_comparison, res)
  }
}

# 1e. Microbiome diversity
if(exists("microbiome_results") && !is.null(microbiome_results$alpha_diversity)) {
  div <- microbiome_results$alpha_diversity
  div_ids <- intersect(ra_common, div$Sample_ID)
  inf_div <- ra_inf_c$Infection[match(div_ids, ra_inf_c$Sample_ID)]
  for(dv in c("Shannon","Simpson")) {
    vals <- div[[dv]][match(div_ids, div$Sample_ID)]
    pos <- vals[inf_div=="Inf_pos"]; neg <- vals[inf_div=="Inf_neg"]
    pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
    res <- compare_inf(pos, neg, dv, "Microbiome")
    if(!is.null(res)) inf_comparison <- rbind(inf_comparison, res)
  }
}

inf_comparison$P_adjusted <- p.adjust(inf_comparison$P_value, method="BH")
inf_comparison <- inf_comparison %>% arrange(P_value)
write.csv(inf_comparison, file.path(out_dir, "tables/Infection_AllOmics_Comparison.csv"), row.names=FALSE)

cat(sprintf("  Total variables tested: %d\n", nrow(inf_comparison)))
cat(sprintf("  Nominal p<0.05: %d\n", sum(inf_comparison$P_value < 0.05)))
cat(sprintf("  BH padj<0.1: %d\n", sum(inf_comparison$P_adjusted < 0.1)))
cat("\n  Top 20:\n")
for(i in 1:min(20, nrow(inf_comparison))) {
  sig <- ifelse(inf_comparison$P_adjusted[i]<0.05,"**FDR",
         ifelse(inf_comparison$P_value[i]<0.01,"**",ifelse(inf_comparison$P_value[i]<0.05,"*","")))
  cat(sprintf("    %-40s [%-15s] delta=%+.3f p=%.4f %s\n",
      inf_comparison$Variable[i], inf_comparison$Category[i],
      inf_comparison$Cliff_delta[i], inf_comparison$P_value[i], sig))
}

# ============================================================================
# SECTION 2: DEG — Future infection(+) vs (-)
# ============================================================================
cat("\n=== SECTION 2: DEG — Infection(+) vs (-) ===\n")

count_inf <- count_matrix[, ra_common]
count_inf <- count_inf[rowSums(count_inf >= 10) >= 3, ]
coldata_inf <- data.frame(row.names=ra_common,
                           Infection=ra_inf_c$Infection[match(ra_common, ra_inf_c$Sample_ID)])

dds_inf <- DESeqDataSetFromMatrix(countData=count_inf, colData=coldata_inf, design=~Infection)
dds_inf <- DESeq(dds_inf)
res_inf <- results(dds_inf, contrast=c("Infection","Inf_pos","Inf_neg"))

# Shrunken LFC
res_inf_shrunk <- lfcShrink(dds_inf, coef="Infection_Inf_pos_vs_Inf_neg", type="normal")
res_inf_df <- as.data.frame(res_inf_shrunk) %>%
  mutate(Gene_Symbol=rownames(res_inf_shrunk)) %>%
  left_join(gene_annotations[,c("Gene_Symbol","ENTREZID","GENENAME")], by="Gene_Symbol") %>%
  arrange(pvalue)

res_inf_df$Regulation <- "NS"
res_inf_df$Regulation[res_inf_df$padj < 0.05 & res_inf_df$log2FoldChange > 0.585] <- "Up_in_InfPos"
res_inf_df$Regulation[res_inf_df$padj < 0.05 & res_inf_df$log2FoldChange < -0.585] <- "Down_in_InfPos"

write.csv(res_inf_df, file.path(out_dir, "tables/DEG_Infection_pos_vs_neg.csv"), row.names=FALSE)
n_sig <- sum(res_inf_df$padj < 0.05, na.rm=TRUE)
cat(sprintf("  DEGs (padj<0.05): %d\n", n_sig))
cat(sprintf("  Up in Inf(+): %d, Down: %d\n",
    sum(res_inf_df$Regulation=="Up_in_InfPos",na.rm=T),
    sum(res_inf_df$Regulation=="Down_in_InfPos",na.rm=T)))

# Top genes
cat("  Top 15 genes by p-value:\n")
top_inf_genes <- head(res_inf_df[!is.na(res_inf_df$pvalue),], 15)
for(i in 1:nrow(top_inf_genes))
  cat(sprintf("    %-15s LFC=%+.3f padj=%s\n",
      top_inf_genes$Gene_Symbol[i], top_inf_genes$log2FoldChange[i],
      ifelse(is.na(top_inf_genes$padj[i]),"NA",sprintf("%.4f",top_inf_genes$padj[i]))))

# ============================================================================
# SECTION 3: GSEA — Infection
# ============================================================================
cat("\n=== SECTION 3: GSEA — Infection ===\n")

ranked_inf <- res_inf_df$log2FoldChange
names(ranked_inf) <- res_inf_df$ENTREZID
ranked_inf <- ranked_inf[!is.na(names(ranked_inf)) & names(ranked_inf)!=""]
ranked_inf <- ranked_inf[!duplicated(names(ranked_inf))]
ranked_inf <- sort(ranked_inf, decreasing=TRUE)

gsea_inf_go <- tryCatch(
  gseGO(geneList=ranked_inf, OrgDb=org.Hs.eg.db, ont="BP",
        pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=15, maxGSSize=500, verbose=FALSE),
  error=function(e) { cat(sprintf("  GO GSEA: %s\n", e$message)); NULL })

if(!is.null(gsea_inf_go) && nrow(gsea_inf_go@result) > 0) {
  write.csv(gsea_inf_go@result, file.path(out_dir, "tables/GSEA_Infection_GO_BP.csv"), row.names=FALSE)
  cat(sprintf("  GO BP: %d enriched (padj<0.05)\n", sum(gsea_inf_go@result$p.adjust<0.05)))
  cat("  Top 10:\n")
  top <- head(gsea_inf_go@result[order(gsea_inf_go@result$p.adjust),], 10)
  for(i in 1:nrow(top))
    cat(sprintf("    NES=%+.2f padj=%.4f %s\n", top$NES[i], top$p.adjust[i], substr(top$Description[i],1,50)))
}

gsea_inf_kegg <- tryCatch(
  gseKEGG(geneList=ranked_inf, organism="hsa", pvalueCutoff=0.1,
          pAdjustMethod="BH", minGSSize=15, maxGSSize=500, verbose=FALSE),
  error=function(e) { cat(sprintf("  KEGG GSEA: %s\n", e$message)); NULL })

if(!is.null(gsea_inf_kegg) && nrow(gsea_inf_kegg@result) > 0) {
  write.csv(gsea_inf_kegg@result, file.path(out_dir, "tables/GSEA_Infection_KEGG.csv"), row.names=FALSE)
  cat(sprintf("  KEGG: %d enriched (padj<0.05)\n", sum(gsea_inf_kegg@result$p.adjust<0.05)))
}

# ============================================================================
# SECTION 4: Microbiome — Infection
# ============================================================================
cat("\n=== SECTION 4: Microbiome — Infection ===\n")

if(exists("abundance_matrix")) {
  micro_ids <- intersect(ra_common, colnames(abundance_matrix))
  inf_micro <- ra_inf_c$Infection[match(micro_ids, ra_inf_c$Sample_ID)]
  cat(sprintf("  Microbiome samples: %d (Inf+:%d, Inf-:%d)\n",
      length(micro_ids), sum(inf_micro=="Inf_pos"), sum(inf_micro=="Inf_neg")))

  if(sum(inf_micro=="Inf_pos") >= 3 && sum(inf_micro=="Inf_neg") >= 3) {
    # Relative abundance
    abund <- abundance_matrix[, micro_ids]
    abund_rel <- sweep(abund, 2, colSums(abund), "/") * 100

    # Filter: present in at least 30% of samples
    prevalence <- rowSums(abund_rel > 0) / ncol(abund_rel)
    abund_filt <- abund_rel[prevalence >= 0.3, ]
    cat(sprintf("  Taxa after prevalence filter (>30%%): %d / %d\n", nrow(abund_filt), nrow(abund_rel)))

    # Differential abundance (Wilcoxon)
    taxa_diff <- data.frame()
    for(taxon in rownames(abund_filt)) {
      pos <- abund_filt[taxon, inf_micro=="Inf_pos"]
      neg <- abund_filt[taxon, inf_micro=="Inf_neg"]
      wt <- tryCatch(wilcox.test(pos, neg, exact=TRUE), error=function(e) NULL)
      if(!is.null(wt)) {
        taxa_diff <- rbind(taxa_diff, data.frame(
          Taxon=taxon,
          Mean_InfPos=mean(pos, na.rm=TRUE), Mean_InfNeg=mean(neg, na.rm=TRUE),
          Median_InfPos=median(pos, na.rm=TRUE), Median_InfNeg=median(neg, na.rm=TRUE),
          Log2FC=log2((mean(pos,na.rm=TRUE)+0.001)/(mean(neg,na.rm=TRUE)+0.001)),
          P_value=wt$p.value))
      }
    }
    taxa_diff$P_adjusted <- p.adjust(taxa_diff$P_value, method="BH")
    taxa_diff <- taxa_diff %>% arrange(P_value)
    write.csv(taxa_diff, file.path(out_dir, "tables/Microbiome_Infection_DiffAbundance.csv"), row.names=FALSE)

    cat(sprintf("  Differentially abundant taxa (p<0.05): %d\n", sum(taxa_diff$P_value<0.05)))
    cat("  Top 10:\n")
    for(i in 1:min(10, nrow(taxa_diff))) {
      sig <- ifelse(taxa_diff$P_value[i]<0.05,"*","")
      # Extract genus/species from taxonomy
      tax_short <- tail(strsplit(taxa_diff$Taxon[i], ";")[[1]], 2)
      cat(sprintf("    LFC=%+.2f p=%.4f %s %s\n",
          taxa_diff$Log2FC[i], taxa_diff$P_value[i], sig, paste(tax_short, collapse=";")))
    }
  }
}

# ============================================================================
# SECTION 5: Multi-omics Infection Prediction Model (LOOCV)
# ============================================================================
cat("\n=== SECTION 5: Multi-omics Infection Prediction ===\n")

# Assemble feature matrix
feature_blocks <- list()

# Expression top 100
expr_inf <- expr_mat[ra_common, ]
expr_vars_inf <- apply(expr_inf, 2, var, na.rm=TRUE)
top_expr <- names(sort(expr_vars_inf, decreasing=TRUE))[1:min(100, sum(!is.na(expr_vars_inf)))]
feature_blocks$Expression <- as.data.frame(expr_inf[, top_expr])

# Cytokines — taken directly from master_data (BALF/Serum prefix retained)
cyto_cols_feat <- grep("^BALF_|^Serum_", colnames(ra_inf_c), value=TRUE)
cyto_cols_feat <- cyto_cols_feat[!grepl("Treg|Th1|Th2|Th17|Macro|Neutro|Lympho|Eosin|Baso|Plasma_Cell|CD[0-9]|Percent|Activated|Ratio|CD14|CD86|CD66", cyto_cols_feat)]
cyto_feat <- ra_inf_c[match(ra_common, ra_inf_c$Sample_ID), cyto_cols_feat]
for(col in colnames(cyto_feat)) cyto_feat[[col]] <- as.numeric(cyto_feat[[col]])
cyto_feat <- log2(cyto_feat + 1)
rownames(cyto_feat) <- ra_common
feature_blocks$Cytokines <- cyto_feat

# FACS — taken directly from master_data (BALF/PB prefix retained)
facs_cols_feat <- grep("^BALF_.*Treg|^BALF_.*Th|^BALF_.*Macro|^BALF_.*Activated|^BALF_.*Neutro.*CD|^PB_", colnames(ra_inf_c), value=TRUE)
facs_feat <- ra_inf_c[match(ra_common, ra_inf_c$Sample_ID), facs_cols_feat]
for(col in colnames(facs_feat)) facs_feat[[col]] <- as.numeric(facs_feat[[col]])
facs_feat <- log2(facs_feat + 0.01)
rownames(facs_feat) <- ra_common
feature_blocks$FACS <- facs_feat

# Deconvolution
if(exists("deconv_mat")) {
  deconv_inf_ids <- intersect(ra_common, rownames(deconv_mat))
  if(length(deconv_inf_ids) == length(ra_common)) {
    deconv_inf <- as.data.frame(deconv_mat[ra_common, ])
    colnames(deconv_inf) <- paste0("Deconv_", colnames(deconv_inf))
    feature_blocks$Deconvolution <- deconv_inf
  }
}

# GSVA
if(exists("gsva_scores")) {
  gsva_inf_ids <- intersect(ra_common, colnames(gsva_scores))
  if(length(gsva_inf_ids) == length(ra_common)) {
    gsva_inf <- as.data.frame(t(gsva_scores[, ra_common]))
    colnames(gsva_inf) <- paste0("GSVA_", colnames(gsva_inf))
    feature_blocks$GSVA <- gsva_inf
  }
}

# Combine
X_inf <- do.call(cbind, feature_blocks)
X_inf[is.na(X_inf)] <- 0
X_inf[!is.finite(as.matrix(X_inf))] <- 0
fv_inf <- apply(X_inf, 2, var, na.rm=TRUE)
X_inf <- X_inf[, !is.na(fv_inf) & fv_inf > 0]
y_inf <- ra_inf_c$Infection[match(ra_common, ra_inf_c$Sample_ID)]

cat(sprintf("  Features: %d\n", ncol(X_inf)))
cat(sprintf("  Samples: %d (Inf+: %d, Inf-: %d)\n",
    length(y_inf), sum(y_inf=="Inf_pos"), sum(y_inf=="Inf_neg")))

# Nested LOOCV with class-weight balancing
set.seed(42)
n_inf_lo <- length(y_inf)
pred_prob <- numeric(n_inf_lo)
pred_class <- character(n_inf_lo)

cat("  Running nested LOOCV...\n")
for(i in 1:n_inf_lo) {
  train_X <- X_inf[-i, ]
  train_y <- y_inf[-i]
  test_X <- X_inf[i, , drop=FALSE]

  # Feature selection within fold (top 30 by Wilcoxon)
  pvals_inner <- apply(train_X, 2, function(x) {
    tryCatch(wilcox.test(x ~ train_y, exact=TRUE)$p.value, error=function(e) 1)
  })
  pvals_inner[is.na(pvals_inner)] <- 1
  n_valid <- sum(pvals_inner < 1, na.rm=TRUE)
  top_feats <- names(sort(pvals_inner))[1:min(30, max(5, n_valid))]

  # Class-weighted RF
  class_weights <- c("Inf_neg"=1, "Inf_pos"=sum(train_y=="Inf_neg")/sum(train_y=="Inf_pos"))
  rf_inner <- randomForest(x=train_X[, top_feats, drop=FALSE], y=train_y,
                            ntree=1000, classwt=class_weights)
  pred_prob[i] <- stats::predict(rf_inner, test_X[, top_feats, drop=FALSE], type="prob")[,"Inf_pos"]
  pred_class[i] <- as.character(predict(rf_inner, test_X[, top_feats, drop=FALSE]))
}

roc_inf_multi <- roc(y_inf, pred_prob, quiet=TRUE)
auc_inf <- as.numeric(auc(roc_inf_multi))
ci_inf <- ci.auc(roc_inf_multi, method="bootstrap", boot.n=2000, quiet=TRUE)

cm_inf <- table(Predicted=pred_class, Actual=y_inf)
accuracy_inf <- sum(diag(cm_inf)) / sum(cm_inf)
sens_inf <- cm_inf["Inf_pos","Inf_pos"] / sum(cm_inf[,"Inf_pos"])
spec_inf <- cm_inf["Inf_neg","Inf_neg"] / sum(cm_inf[,"Inf_neg"])

cat(sprintf("\n  Multi-omics Infection Prediction:\n"))
cat(sprintf("    AUC: %.3f (95%% CI: %.3f-%.3f)\n", auc_inf, ci_inf[1], ci_inf[3]))
cat(sprintf("    Accuracy: %.3f\n", accuracy_inf))
cat(sprintf("    Sensitivity: %.3f\n", sens_inf))
cat(sprintf("    Specificity: %.3f\n", spec_inf))
print(cm_inf)

# Feature importance (full model)
pvals_full <- apply(X_inf, 2, function(x) {
  tryCatch(wilcox.test(x ~ y_inf, exact=TRUE)$p.value, error=function(e) 1)
})
top_full <- names(sort(pvals_full))[1:min(30, length(pvals_full))]
rf_full <- randomForest(x=X_inf[, top_full], y=y_inf, ntree=1000,
                         classwt=c("Inf_neg"=1, "Inf_pos"=sum(y_inf=="Inf_neg")/sum(y_inf=="Inf_pos")),
                         importance=TRUE)
imp_inf <- data.frame(Feature=rownames(importance(rf_full)),
                       MeanDecreaseAccuracy=importance(rf_full)[,"MeanDecreaseAccuracy"],
                       MeanDecreaseGini=importance(rf_full)[,"MeanDecreaseGini"]) %>%
  arrange(desc(MeanDecreaseGini))
imp_inf$Type <- case_when(
  grepl("^Deconv_", imp_inf$Feature) ~ "Deconvolution",
  grepl("^GSVA_", imp_inf$Feature) ~ "GSVA",
  grepl("^BALF_|^Serum_|^PB_", imp_inf$Feature) ~ "Cytokine/FACS",
  TRUE ~ "Expression"
)
write.csv(imp_inf, file.path(out_dir, "tables/Infection_Multiomics_FeatureImportance.csv"), row.names=FALSE)

cat("\n  Top 15 features:\n")
for(i in 1:min(15, nrow(imp_inf)))
  cat(sprintf("    %-40s [%-15s] Gini=%.3f\n", imp_inf$Feature[i], imp_inf$Type[i], imp_inf$MeanDecreaseGini[i]))

# ============================================================================
# SECTION 6: Permutation test
# ============================================================================
cat("\n=== SECTION 6: Permutation Test ===\n")
set.seed(42)
n_perm <- 500
perm_aucs <- numeric(n_perm)
for(p in 1:n_perm) {
  perm_y <- sample(y_inf)
  perm_prob <- numeric(n_inf_lo)
  for(i in 1:n_inf_lo) {
    train_X <- X_inf[-i,]; train_y_p <- perm_y[-i]
    pvals_p <- apply(train_X, 2, function(x) tryCatch(wilcox.test(x~train_y_p,exact=TRUE)$p.value,error=function(e)1))
    pvals_p[is.na(pvals_p)] <- 1
    tp <- names(sort(pvals_p))[1:min(30, max(5, sum(pvals_p<1,na.rm=TRUE)))]
    cw <- c("Inf_neg"=1,"Inf_pos"=max(1,sum(train_y_p=="Inf_neg")/max(1,sum(train_y_p=="Inf_pos"))))
    rf_p <- randomForest(x=train_X[,tp,drop=FALSE],y=train_y_p,ntree=500,classwt=cw)
    perm_prob[i] <- stats::predict(rf_p, X_inf[i,tp,drop=FALSE], type="prob")[,"Inf_pos"]
  }
  perm_roc <- tryCatch(roc(y_inf, perm_prob, quiet=TRUE), error=function(e) NULL)
  perm_aucs[p] <- if(!is.null(perm_roc)) as.numeric(auc(perm_roc)) else 0.5
  if(p %% 100 == 0) cat(sprintf("    %d/%d\n", p, n_perm))
}
perm_p <- (sum(perm_aucs >= auc_inf) + 1) / (n_perm + 1)
cat(sprintf("  Permutation p-value: %.4f\n", perm_p))

# ============================================================================
# SECTION 7: Figure — Infection Prediction
# ============================================================================
cat("\n=== SECTION 7: Generating Figures ===\n")

suppressPackageStartupMessages(library(ggplot2))
theme_natmed <- function(base_size=10) {
  theme_classic(base_size=base_size) %+replace%
    theme(text=element_text(family="Helvetica",color="black"),
          axis.text=element_text(size=base_size-1,color="black"),
          axis.title=element_text(size=base_size,face="bold"),
          plot.title=element_text(size=base_size+1,face="bold",hjust=0),
          panel.border=element_rect(fill=NA,color="black",linewidth=0.5),
          axis.line=element_blank(),
          plot.margin=unit(c(5,5,5,5),"mm"))
}

# 7a ROC
roc_df <- data.frame(Sens=roc_inf_multi$sensitivities, Spec=1-roc_inf_multi$specificities)
p7a <- ggplot(roc_df, aes(Spec, Sens)) +
  geom_line(color="#C0392B", linewidth=1) +
  geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60") +
  annotate("text",x=0.55,y=0.15,size=3,hjust=0,
           label=sprintf("AUC = %.3f\n95%% CI: %.3f\u2013%.3f\nPerm p = %.4f\nn = %d (Inf+:%d)",
                          auc_inf,ci_inf[1],ci_inf[3],perm_p,n_inf_lo,sum(y_inf=="Inf_pos"))) +
  labs(x="1 - Specificity",y="Sensitivity",title="a  Infection prediction (multi-omics LOOCV)") +
  coord_equal() + theme_natmed()

# 7b Top discriminators
top_disc <- inf_comparison %>% filter(P_value < 0.1) %>%
  mutate(Label=gsub("^BALF_|^Serum_|^GSVA_|^Deconv_","",Variable)) %>%
  arrange(Cliff_delta)
if(nrow(top_disc) > 0) {
  top_disc$Label <- fct_inorder(top_disc$Label)
  p7b <- ggplot(top_disc, aes(Cliff_delta, Label, fill=Category)) +
    geom_col(width=0.6) +
    scale_fill_brewer(palette="Set2") +
    geom_vline(xintercept=0, linewidth=0.4) +
    labs(x="Cliff's delta (Inf+ vs Inf-)", y="", title="b  Infection-associated features") +
    theme_natmed(base_size=8) + theme(legend.position="bottom", legend.key.size=unit(3,"mm"))
}

# 7c GSVA boxplots (top pathways)
gsva_inf_comp <- inf_comparison %>% filter(Category=="GSVA_Pathway") %>% arrange(P_value)
if(nrow(gsva_inf_comp) >= 3) {
  top_gs <- head(gsva_inf_comp$Variable, 6)
  gsva_box_data <- data.frame()
  for(gs in top_gs) {
    vals <- gsva_scores[gs, ra_common]
    gsva_box_data <- rbind(gsva_box_data, data.frame(
      Pathway=gs, Score=vals, Infection=y_inf))
  }
  p7c <- ggplot(gsva_box_data, aes(Pathway, Score, fill=Infection)) +
    geom_boxplot(outlier.size=0.8, width=0.6, linewidth=0.3) +
    scale_fill_manual(values=c("Inf_neg"="#3498DB","Inf_pos"="#E74C3C"),
                      labels=c("Inf_neg"="No infection","Inf_pos"="Future infection")) +
    labs(x="", y="GSVA score", title="c  Pathway activity by infection outcome") +
    theme_natmed(base_size=8) + theme(axis.text.x=element_text(angle=30,hjust=1),
                                       legend.position=c(0.8,0.9))
}

# 7d Permutation histogram
perm_df <- data.frame(AUC=perm_aucs)
p7d <- ggplot(perm_df, aes(AUC)) +
  geom_histogram(bins=30, fill="grey70", color="grey50", linewidth=0.2) +
  geom_vline(xintercept=auc_inf, color="#C0392B", linewidth=1, linetype="dashed") +
  annotate("text", x=auc_inf, y=Inf, vjust=2, hjust=-0.1, size=3, color="#C0392B",
           label=sprintf("Observed\nAUC=%.3f", auc_inf)) +
  labs(x="Permuted AUC", y="Count", title="d  Permutation distribution") +
  theme_natmed()

ggsave(file.path(out_dir, "figures/Fig_Infection_Prediction.pdf"),
       gridExtra::arrangeGrob(p7a, p7b, p7c, p7d, ncol=2),
       width=180, height=170, units="mm", dpi=300)
cat("  ✓ Fig_Infection_Prediction.pdf\n")

# ============================================================================
# SECTION 8: Summary
# ============================================================================
cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  Infection Analysis Complete                                  ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Samples: %d RA (Inf+:%d, Inf-:%d)\n",
    n_inf_lo, sum(y_inf=="Inf_pos"), sum(y_inf=="Inf_neg")))
cat(sprintf("  All-omics comparison: %d tested, %d nominal p<0.05\n",
    nrow(inf_comparison), sum(inf_comparison$P_value<0.05)))
cat(sprintf("  DEGs (Inf+ vs Inf-): %d (padj<0.05)\n", n_sig))
cat(sprintf("  GSEA GO: %d terms\n",
    if(!is.null(gsea_inf_go)) sum(gsea_inf_go@result$p.adjust<0.05) else 0))
cat(sprintf("  Multi-omics LOOCV AUC: %.3f (CI: %.3f-%.3f, perm p=%.4f)\n",
    auc_inf, ci_inf[1], ci_inf[3], perm_p))

save(inf_comparison, res_inf_df, imp_inf, auc_inf, ci_inf, perm_p,
     file=file.path(out_dir, "Infection_Prediction_Results.RData"))

# ============================================================================
# SECTION 8: Single Biomarker ROC (for Figure 5a)
# ============================================================================
cat("\n=== SECTION 8: Single Biomarker ROC ===\n")

ra_c_roc <- master_data[master_data$Sample_Group=="RA" & !is.na(master_data$respiratory_infection) &
                          master_data$Sample_ID %in% common_samples, ]
y_roc <- factor(ifelse(ra_c_roc$respiratory_infection==1,"Inf_pos","Inf_neg"), levels=c("Inf_neg","Inf_pos"))

biomarkers_roc <- list()
for(col in c("BALF_Th17.1_CCR6pos_CXCR3pos","BALF_Th17_CCR6pos_CXCR3neg",
             "BALF_Treg_CD127neg_CD25pos","BALF_Th1_CCR6neg_CXCR3pos",
             "BALF_Activated_Macrophage_CD86pos_CD14pos","BALF_MCP3","BALF_GMCSF","BALF_MDC",
             "Serum_IL6"))
  if(col %in% colnames(ra_c_roc)) biomarkers_roc[[col]] <- as.numeric(ra_c_roc[[col]])

# BAL/PB ratios
if("BALF_Th17.1_CCR6pos_CXCR3pos" %in% colnames(ra_c_roc) & "PB_Th17.1_CCR6pos_CXCR3pos" %in% colnames(ra_c_roc))
  biomarkers_roc[["Ratio_Th17.1_BAL_PB"]] <- as.numeric(ra_c_roc[["BALF_Th17.1_CCR6pos_CXCR3pos"]]) /
    (as.numeric(ra_c_roc[["PB_Th17.1_CCR6pos_CXCR3pos"]]) + 0.01)
if("BALF_Treg_CD127neg_CD25pos" %in% colnames(ra_c_roc) & "PB_Treg_CD127neg_CD25pos" %in% colnames(ra_c_roc))
  biomarkers_roc[["Ratio_Treg_BAL_PB"]] <- as.numeric(ra_c_roc[["BALF_Treg_CD127neg_CD25pos"]]) /
    (as.numeric(ra_c_roc[["PB_Treg_CD127neg_CD25pos"]]) + 0.01)

# GSVA
if(exists("gsva_scores")) for(gs in c("OXPHOS","TCA_Cycle","Autophagy"))
  if(gs %in% rownames(gsva_scores) & all(ra_c_roc$Sample_ID %in% colnames(gsva_scores)))
    biomarkers_roc[[paste0("GSVA_",gs)]] <- gsva_scores[gs, ra_c_roc$Sample_ID]

# Deconv
if(exists("deconv_mat") & all(ra_c_roc$Sample_ID %in% rownames(deconv_mat)))
  for(ct in c("NK","DC","Mast","Neutrophil","Macrophage"))
    biomarkers_roc[[paste0("Deconv_",ct)]] <- deconv_mat[ra_c_roc$Sample_ID, ct]

# Gene
if(all(ra_c_roc$Sample_ID %in% colnames(log_cpm)))
  for(gene in c("SFTPB","ETV5","NDUFS5","CA2"))
    if(gene %in% rownames(log_cpm)) biomarkers_roc[[paste0("Gene_",gene)]] <- log_cpm[gene, ra_c_roc$Sample_ID]

roc_results <- data.frame()
for(nm in names(biomarkers_roc)) {
  vv <- biomarkers_roc[[nm]]; valid <- !is.na(vv) & !is.na(y_roc)
  if(sum(valid) >= 10 & length(unique(y_roc[valid]))==2) {
    ro <- tryCatch(roc(y_roc[valid], vv[valid], quiet=TRUE, direction="auto"), error=function(e) NULL)
    if(!is.null(ro)) {
      auc_v <- as.numeric(auc(ro))
      ci_v <- tryCatch(ci.auc(ro, method="bootstrap", boot.n=2000, quiet=TRUE), error=function(e) c(NA,NA,NA))
      mn_pos <- mean(vv[valid & y_roc=="Inf_pos"], na.rm=TRUE)
      mn_neg <- mean(vv[valid & y_roc=="Inf_neg"], na.rm=TRUE)
      roc_results <- rbind(roc_results, data.frame(
        Biomarker=nm, AUC=auc_v, CI_lo=ci_v[1], CI_hi=ci_v[3],
        Direction=ifelse(mn_pos>mn_neg,"Up in Inf+","Down in Inf+"), N=sum(valid)))
    }
  }
}
roc_results <- roc_results %>% arrange(desc(AUC))
write.csv(roc_results, file.path(out_dir, "tables/Infection_SingleBiomarker_ROC.csv"), row.names=FALSE)
cat(sprintf("  Single biomarker ROC: %d biomarkers computed\n", nrow(roc_results)))
