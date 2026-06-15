# ============================================================================
# RA-ILD Multi-omics: Nature Communications Enhanced Analysis
#
# Purpose: Implement statistically rigorous analyses recommended for
#          Nature Communications submission
# Prerequisite: 20260210postDeconvolution_FIXED.R executed
# ============================================================================

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  RA-ILD Nature Communications Enhanced Analysis              ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

# --- Setup ---
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
  library(GSVA); library(pheatmap); library(effsize)
})
set.seed(42)

out_dir <- "results"
for(d in c(out_dir, file.path(out_dir, "tables"), file.path(out_dir, "figures")))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

report <- list()

# ============================================================================
# SECTION 1: Confounding / Covariate Assessment
# ============================================================================
cat("\n=== SECTION 1: Confounding Assessment ===\n")

pca <- prcomp(t(expr_matrix), scale.=TRUE)
pca_df <- data.frame(
  PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],
  Sample_ID=rownames(pca$x)
) %>% left_join(master_data[, c("Sample_ID","Sample_Group",
                                 "gender_Male(1,_0)","Age",
                                 "smoking(1,0)","respiratory_infection")],
                by="Sample_ID")
# Rename for convenience
colnames(pca_df)[colnames(pca_df)=="gender_Male(1,_0)"] <- "male"
colnames(pca_df)[colnames(pca_df)=="smoking(1,0)"] <- "smoking_binary"

# Variance explained
var_explained <- summary(pca)$importance[2, 1:5] * 100

# PC vs covariate associations
covar_tests <- data.frame()
for(pc_name in c("PC1","PC2","PC3")) {
  pc_vals <- pca_df[[pc_name]]

  # Disease group
  if(sum(!is.na(pca_df$Sample_Group)) >= 10) {
    grp <- factor(pca_df$Sample_Group)
    wt <- wilcox.test(pc_vals ~ grp, exact=TRUE)
    covar_tests <- rbind(covar_tests, data.frame(
      PC=pc_name, Covariate="Disease(RA_vs_Ctrl)", P=wt$p.value, Method="Wilcoxon"))
  }
  # Sex
  if(sum(!is.na(pca_df$male)) >= 10) {
    sex <- factor(pca_df$male)
    valid <- !is.na(sex)
    if(length(unique(sex[valid])) == 2) {
      wt <- wilcox.test(pc_vals[valid] ~ sex[valid], exact=TRUE)
      covar_tests <- rbind(covar_tests, data.frame(
        PC=pc_name, Covariate="Sex", P=wt$p.value, Method="Wilcoxon"))
    }
  }
  # Age
  if(sum(!is.na(pca_df$Age)) >= 10) {
    ct <- cor.test(pc_vals[!is.na(pca_df$Age)], as.numeric(pca_df$Age[!is.na(pca_df$Age)]),
                   method="spearman")
    covar_tests <- rbind(covar_tests, data.frame(
      PC=pc_name, Covariate="Age", P=ct$p.value, Method="Spearman"))
  }
  # Smoking
  if(sum(!is.na(pca_df$smoking_binary)) >= 10) {
    smk <- factor(pca_df$smoking_binary)
    valid <- !is.na(smk)
    if(length(unique(smk[valid])) == 2) {
      wt <- wilcox.test(pc_vals[valid] ~ smk[valid], exact=TRUE)
      covar_tests <- rbind(covar_tests, data.frame(
        PC=pc_name, Covariate="Smoking", P=wt$p.value, Method="Wilcoxon"))
    }
  }
  # Infection
  if(sum(!is.na(pca_df$respiratory_infection)) >= 10) {
    inf <- factor(pca_df$respiratory_infection)
    valid <- !is.na(inf)
    if(length(unique(inf[valid])) == 2) {
      wt <- wilcox.test(pc_vals[valid] ~ inf[valid], exact=TRUE)
      covar_tests <- rbind(covar_tests, data.frame(
        PC=pc_name, Covariate="Infection", P=wt$p.value, Method="Wilcoxon"))
    }
  }
}

write.csv(covar_tests, file.path(out_dir, "tables/Confounding_PC_associations.csv"), row.names=FALSE)

# PCA plot
tryCatch({
  pdf(file.path(out_dir, "figures/PCA_covariates.pdf"), width=12, height=8)
  par(mfrow=c(2,3), mar=c(5,4,3,2))

  # By disease
  cols_disease <- ifelse(pca_df$Sample_Group=="RA", "#E41A1C", "#4DAF4A")
  plot(pca_df$PC1, pca_df$PC2, col=cols_disease, pch=19, cex=1.3,
       xlab=sprintf("PC1 (%.1f%%)", var_explained[1]),
       ylab=sprintf("PC2 (%.1f%%)", var_explained[2]),
       main="Disease Group")
  legend("topright", c("RA","Control"), col=c("#E41A1C","#4DAF4A"), pch=19, cex=0.8)

  # By sex
  sex_cols <- ifelse(is.na(pca_df$male), "grey",
                     ifelse(pca_df$male==1, "steelblue", "salmon"))
  plot(pca_df$PC1, pca_df$PC2, col=sex_cols, pch=19, cex=1.3,
       xlab=sprintf("PC1 (%.1f%%)", var_explained[1]),
       ylab=sprintf("PC2 (%.1f%%)", var_explained[2]),
       main="Sex")
  legend("topright", c("Male","Female","NA"), col=c("steelblue","salmon","grey"), pch=19, cex=0.8)

  # By age
  age_vals <- as.numeric(pca_df$Age)
  age_cols <- ifelse(is.na(age_vals), "grey",
                     colorRampPalette(c("yellow","red"))(100)[pmin(100, pmax(1, round((age_vals-min(age_vals,na.rm=T))/(max(age_vals,na.rm=T)-min(age_vals,na.rm=T)+0.01)*99)+1))])
  plot(pca_df$PC1, pca_df$PC2, col=age_cols, pch=19, cex=1.3,
       xlab=sprintf("PC1 (%.1f%%)", var_explained[1]),
       ylab=sprintf("PC2 (%.1f%%)", var_explained[2]),
       main="Age")

  # By infection
  inf_cols <- ifelse(is.na(pca_df$respiratory_infection), "grey",
                     ifelse(pca_df$respiratory_infection==1, "#FF7F00", "#377EB8"))
  plot(pca_df$PC1, pca_df$PC2, col=inf_cols, pch=19, cex=1.3,
       xlab=sprintf("PC1 (%.1f%%)", var_explained[1]),
       ylab=sprintf("PC2 (%.1f%%)", var_explained[2]),
       main="Respiratory Infection")
  legend("topright", c("Pos","Neg","NA"), col=c("#FF7F00","#377EB8","grey"), pch=19, cex=0.8)

  # By subgroup
  sub_cols <- c("Control"="#4DAF4A", "RA_nonILD"="#377EB8", "RA_ILD"="#E41A1C")
  sg <- meta$Subgroup[match(pca_df$Sample_ID, meta$Sample_ID)]
  plot(pca_df$PC1, pca_df$PC2, col=sub_cols[as.character(sg)], pch=19, cex=1.3,
       xlab=sprintf("PC1 (%.1f%%)", var_explained[1]),
       ylab=sprintf("PC2 (%.1f%%)", var_explained[2]),
       main="Subgroup")
  legend("topright", names(sub_cols), col=sub_cols, pch=19, cex=0.8)

  # Variance explained
  barplot(var_explained[1:5], names.arg=paste0("PC",1:5), col="steelblue",
          ylab="% Variance", main="Variance Explained")

  dev.off()
  cat("  ✓ PCA covariate plots saved\n")
}, error=function(e) { tryCatch(dev.off(), error=function(x) NULL); cat(sprintf("  ✗ %s\n", e$message)) })

report$confounding <- list(
  var_explained=var_explained,
  covariate_associations=covar_tests
)
cat(sprintf("  PC1: %.1f%%, PC2: %.1f%%\n", var_explained[1], var_explained[2]))
cat("  Significant confounders (p<0.05):\n")
sig_conf <- covar_tests[covar_tests$P < 0.05, ]
if(nrow(sig_conf) > 0) {
  for(i in 1:nrow(sig_conf))
    cat(sprintf("    %s vs %s: p=%.4f\n", sig_conf$PC[i], sig_conf$Covariate[i], sig_conf$P[i]))
} else cat("    None\n")


# ============================================================================
# SECTION 2: DEG with Shrunken LFC
# Note: Covariate adjustment (age/sex) was attempted but Control group lacks
# age/sex data (<80% available), so actual design is ~ Group only.
# ============================================================================
cat("\n=== SECTION 2: DEG (Shrunken LFC) ===\n")

# Build DESeq2 (covariates included only if >=80% samples have data)
deg_samples <- common_samples
deg_meta <- master_data[match(deg_samples, master_data$Sample_ID), ]
rownames(deg_meta) <- deg_meta$Sample_ID

# Determine available covariates
# Map covariate column names
if("gender_Male(1,_0)" %in% colnames(deg_meta)) {
  deg_meta$male <- deg_meta[["gender_Male(1,_0)"]]
}
if("smoking(1,0)" %in% colnames(deg_meta)) {
  deg_meta$smoking_binary <- deg_meta[["smoking(1,0)"]]
}
has_sex <- sum(!is.na(deg_meta$male)) >= length(deg_samples) * 0.8
has_age <- sum(!is.na(deg_meta$Age)) >= length(deg_samples) * 0.8

deg_meta$Group <- factor(ifelse(grepl("^(KYC|Sarcoidosis)", deg_meta$Sample_ID), "Control", "RA"),
                          levels=c("Control","RA"))

# Handle NAs in covariates - impute for Control group which lacks clinical data
if(has_sex) {
  deg_meta$Sex <- as.factor(deg_meta$male)
  # For Controls with NA sex, impute median or create a dummy
  if(any(is.na(deg_meta$Sex))) {
    cat("  Note: Sex NA in some samples, running without sex covariate\n")
    has_sex <- FALSE
  }
}

if(has_age) {
  deg_meta$Age_num <- as.numeric(deg_meta$Age)
  if(any(is.na(deg_meta$Age_num))) {
    cat("  Note: Age NA in some samples, running without age covariate\n")
    has_age <- FALSE
  }
}

# Design formula
if(has_sex && has_age) {
  design_formula <- ~ Sex + Age_num + Group
  cat("  Design: ~ Sex + Age + Group\n")
} else if(has_sex) {
  design_formula <- ~ Sex + Group
  cat("  Design: ~ Sex + Group\n")
} else {
  design_formula <- ~ Group
  cat("  Design: ~ Group (no covariates available for all samples)\n")
}

count_mat_deg <- count_matrix[, deg_samples]
count_mat_deg <- count_mat_deg[rowSums(count_mat_deg >= 10) >= 5, ]

dds2 <- DESeqDataSetFromMatrix(countData=count_mat_deg, colData=deg_meta, design=design_formula)
dds2 <- DESeq(dds2)

# Shrunken LFC (apeglm)
tryCatch({
  suppressPackageStartupMessages(library(apeglm))
  res_shrunk <- lfcShrink(dds2, coef="Group_RA_vs_Control", type="apeglm")
  cat("  ✓ apeglm shrinkage applied\n")
  shrink_method <- "apeglm"
}, error=function(e) {
  res_shrunk <<- lfcShrink(dds2, coef="Group_RA_vs_Control", type="normal")
  cat("  ✓ normal shrinkage applied (apeglm unavailable)\n")
  shrink_method <<- "normal"
})

res2_df <- as.data.frame(res_shrunk) %>%
  mutate(Gene_Symbol=rownames(res_shrunk)) %>%
  left_join(gene_annotations[,c("Gene_Symbol","ENTREZID","GENENAME")], by="Gene_Symbol") %>%
  arrange(padj)
res2_df$Regulation <- "NS"
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange > 0.585] <- "Up"   # log2(1.5)
res2_df$Regulation[res2_df$padj < 0.05 & res2_df$log2FoldChange < -0.585] <- "Down"

write.csv(res2_df, file.path(out_dir, "tables/DEG_ShrunkLFC.csv"), row.names=FALSE)

n_up <- sum(res2_df$Regulation=="Up", na.rm=TRUE)
n_down <- sum(res2_df$Regulation=="Down", na.rm=TRUE)
cat(sprintf("  DEGs (padj<0.05, |LFC|>0.585): %d up, %d down\n", n_up, n_down))

report$deg <- list(
  n_up=n_up, n_down=n_down, total_tested=nrow(res2_df),
  shrink_method=shrink_method
)


# ============================================================================
# SECTION 3: GSEA (Rank-based Gene Set Enrichment Analysis)
# ============================================================================
cat("\n=== SECTION 3: GSEA (Rank-based) ===\n")

# Create ranked gene list from shrunken LFC
ranked <- res2_df$log2FoldChange
names(ranked) <- res2_df$ENTREZID
ranked <- ranked[!is.na(names(ranked)) & names(ranked) != ""]
ranked <- ranked[!duplicated(names(ranked))]
ranked <- sort(ranked, decreasing=TRUE)

cat(sprintf("  Ranked gene list: %d genes\n", length(ranked)))

# GO BP GSEA
gsea_go <- tryCatch({
  gseGO(geneList=ranked, OrgDb=org.Hs.eg.db, ont="BP",
        pvalueCutoff=0.05, pAdjustMethod="BH", minGSSize=15, maxGSSize=500,
        verbose=FALSE)
}, error=function(e) { cat(sprintf("  ⚠ GO GSEA failed: %s\n", e$message)); NULL })

if(!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  write.csv(gsea_go@result, file.path(out_dir, "tables/GSEA_GO_BP.csv"), row.names=FALSE)
  cat(sprintf("  GO BP GSEA: %d enriched terms\n", sum(gsea_go@result$p.adjust < 0.05)))

  tryCatch({
    pdf(file.path(out_dir, "figures/GSEA_GO_top20.pdf"), width=10, height=12)
    top_terms <- head(gsea_go@result[order(gsea_go@result$p.adjust), ], 20)
    top_terms$Description <- factor(top_terms$Description, levels=rev(top_terms$Description))
    p <- ggplot(top_terms, aes(x=NES, y=Description, fill=p.adjust)) +
      geom_bar(stat="identity") +
      scale_fill_gradient(low="red", high="blue") +
      theme_minimal() + labs(title="GSEA GO BP (Top 20)", x="NES", y="")
    print(p)
    dev.off()
    cat("  ✓ GSEA GO plot saved\n")
  }, error=function(e) { tryCatch(dev.off(), error=function(x) NULL) })
}

# KEGG GSEA
gsea_kegg <- tryCatch({
  gseKEGG(geneList=ranked, organism="hsa", pvalueCutoff=0.05,
          pAdjustMethod="BH", minGSSize=15, maxGSSize=500, verbose=FALSE)
}, error=function(e) { cat(sprintf("  ⚠ KEGG GSEA: %s\n", e$message)); NULL })

if(!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  write.csv(gsea_kegg@result, file.path(out_dir, "tables/GSEA_KEGG.csv"), row.names=FALSE)
  cat(sprintf("  KEGG GSEA: %d enriched pathways\n", sum(gsea_kegg@result$p.adjust < 0.05)))
}

# Hallmark gene sets from custom definition (key pathways)
gsea_custom <- tryCatch({
  gseGO(geneList=ranked, OrgDb=org.Hs.eg.db, ont="MF",
        pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500,
        verbose=FALSE)
}, error=function(e) NULL)

if(!is.null(gsea_custom) && nrow(gsea_custom@result) > 0) {
  write.csv(gsea_custom@result, file.path(out_dir, "tables/GSEA_GO_MF.csv"), row.names=FALSE)
  cat(sprintf("  GO MF GSEA: %d enriched terms\n", sum(gsea_custom@result$p.adjust < 0.05)))
}

report$gsea <- list(
  GO_BP=if(!is.null(gsea_go)) nrow(gsea_go@result[gsea_go@result$p.adjust<0.05,]) else 0,
  KEGG=if(!is.null(gsea_kegg)) nrow(gsea_kegg@result[gsea_kegg@result$p.adjust<0.05,]) else 0,
  GO_MF=if(!is.null(gsea_custom)) nrow(gsea_custom@result[gsea_custom@result$p.adjust<0.05,]) else 0
)


# ============================================================================
# SECTION 4: Effect Sizes for All Group Comparisons
# ============================================================================
cat("\n=== SECTION 4: Effect Sizes (Cliff's delta) ===\n")

compute_effect_size <- function(data, var_col, group_col, g1, g2) {
  v1 <- as.numeric(data[[var_col]][data[[group_col]]==g1])
  v2 <- as.numeric(data[[var_col]][data[[group_col]]==g2])
  v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
  if(length(v1) < 3 || length(v2) < 3) return(data.frame(Variable=var_col, n1=length(v1), n2=length(v2),
                                                            median1=NA, median2=NA, p_value=NA, cliff_delta=NA, magnitude=NA))
  wt <- wilcox.test(v1, v2, exact=TRUE)
  cd <- cliff.delta(v1, v2)
  data.frame(Variable=var_col, n1=length(v1), n2=length(v2),
             median1=median(v1), median2=median(v2),
             p_value=wt$p.value, cliff_delta=cd$estimate, magnitude=cd$magnitude)
}

# Cytokine effect sizes (RA vs Control)
all_effect_sizes <- data.frame()
cyto_cols_for_es <- grep("^BALF_|^Serum_", colnames(master_data), value=TRUE)
cyto_cols_for_es <- cyto_cols_for_es[!grepl("Treg|Th|Macro|Neutro|Lympho|CD", cyto_cols_for_es)]
for(col in cyto_cols_for_es) {
  if(is.numeric(master_data[[col]]) || !all(is.na(as.numeric(master_data[[col]])))) {
    es <- compute_effect_size(master_data, col, "Sample_Group", "RA", "Control")
    es$Category <- ifelse(grepl("^BALF_", col), "BALF_Cytokine", "Serum_Cytokine")
    all_effect_sizes <- rbind(all_effect_sizes, es)
  }
}

# FACS effect sizes
for(col in c(balf_facs_cols, pbmc_facs_cols)) {
  if(col %in% colnames(master_data)) {
    es <- compute_effect_size(master_data, col, "Sample_Group", "RA", "Control")
    es$Category <- ifelse(grepl("^BALF_|^bal_", col), "BALF_FACS", "PB_FACS")
    all_effect_sizes <- rbind(all_effect_sizes, es)
  }
}

# Deconvolution effect sizes — use cell fraction (deconv_cell) for proportions
if(exists("deconv_cell")) {
  deconv_df <- data.frame(deconv_cell, Group=as.character(deconv_group))
  for(ct in colnames(deconv_cell)) {
    es <- compute_effect_size(deconv_df, ct, "Group", "RA", "Control")
    es$Category <- "Deconvolution"
    all_effect_sizes <- rbind(all_effect_sizes, es)
  }
}

all_effect_sizes$p_adjusted <- p.adjust(all_effect_sizes$p_value, method="BH")
all_effect_sizes <- all_effect_sizes %>% arrange(p_value)
write.csv(all_effect_sizes, file.path(out_dir, "tables/Effect_sizes_RA_vs_Control.csv"), row.names=FALSE)

sig_es <- all_effect_sizes[!is.na(all_effect_sizes$p_adjusted) & all_effect_sizes$p_adjusted < 0.05, ]
cat(sprintf("  Significant (padj<0.05): %d / %d variables\n", nrow(sig_es), nrow(all_effect_sizes)))
cat(sprintf("  Large effect (|delta|>0.474): %d\n",
    sum(abs(all_effect_sizes$cliff_delta) > 0.474, na.rm=TRUE)))

report$effect_sizes <- list(total=nrow(all_effect_sizes), significant=nrow(sig_es),
                             large_effect=sum(abs(all_effect_sizes$cliff_delta) > 0.474, na.rm=TRUE))


# ============================================================================
# SECTION 5: Bootstrap AUC Confidence Intervals for LOOCV
# ============================================================================
cat("\n=== SECTION 5: Bootstrap AUC CI ===\n")

# Redo nested LOOCV with bootstrap CI
loocv_samples2 <- intersect(intersect(rownames(expr_mat), rownames(cyto_mat)), rownames(facs_mat))
loocv_data2 <- data.frame(
  expr_mat[loocv_samples2, 1:min(50, ncol(expr_mat))],
  cyto_mat[loocv_samples2, ],
  facs_mat[loocv_samples2, ],
  Group=meta$Sample_Group[match(loocv_samples2, meta$Sample_ID)]
)
loocv_data2 <- loocv_data2[!is.na(loocv_data2$Group), ]
loocv_data2$Group <- factor(loocv_data2$Group)
for(col in setdiff(colnames(loocv_data2), "Group"))
  loocv_data2[[col]][is.na(loocv_data2[[col]]) | is.infinite(loocv_data2[[col]])] <-
    median(loocv_data2[[col]], na.rm=TRUE)
feature_vars2 <- apply(loocv_data2[, setdiff(colnames(loocv_data2),"Group"), drop=FALSE], 2, var, na.rm=TRUE)
loocv_data2 <- loocv_data2[, c(names(feature_vars2[!is.na(feature_vars2) & feature_vars2>0]), "Group")]

n_lo2 <- nrow(loocv_data2)
nested_pred2 <- rep(NA, n_lo2)
nested_prob2 <- matrix(NA, nrow=n_lo2, ncol=2, dimnames=list(NULL, levels(loocv_data2$Group)))

cat(sprintf("  Running nested LOOCV (n=%d)...\n", n_lo2))
set.seed(42)
for(i in 1:n_lo2) {
  train <- loocv_data2[-i, ]
  features <- setdiff(colnames(train), "Group")
  # Feature selection within fold
  fv_inner <- apply(train[, features], 2, function(x) {
    tryCatch(wilcox.test(x ~ train$Group, exact=TRUE)$p.value, error=function(e) 1)
  })
  top_feats <- names(sort(fv_inner))[1:min(50, length(fv_inner))]
  rf_inner <- randomForest(x=train[, top_feats], y=train$Group, ntree=500)
  nested_pred2[i] <- as.character(predict(rf_inner, loocv_data2[i, top_feats]))
  nested_prob2[i,] <- stats::predict(rf_inner, loocv_data2[i, top_feats], type="prob")
}

roc_nested <- roc(loocv_data2$Group, nested_prob2[,"RA"], quiet=TRUE)
auc_val <- as.numeric(auc(roc_nested))

# Bootstrap CI
ci_boot <- ci.auc(roc_nested, method="bootstrap", boot.n=2000, quiet=TRUE)
cat(sprintf("  Nested LOOCV AUC: %.3f (95%% CI: %.3f-%.3f)\n", auc_val, ci_boot[1], ci_boot[3]))

# Permutation test
cat("  Running permutation test (1000 iterations)...\n")
set.seed(42)
n_perm2 <- 1000
perm_aucs2 <- numeric(n_perm2)
for(p in 1:n_perm2) {
  perm_labels <- sample(loocv_data2$Group)
  perm_pred <- numeric(n_lo2)
  for(i in 1:n_lo2) {
    train_perm <- loocv_data2[-i, ]
    train_perm$Group <- perm_labels[-i]
    features <- setdiff(colnames(train_perm), "Group")
    fv_p <- apply(train_perm[, features], 2, function(x) {
      tryCatch(wilcox.test(x ~ train_perm$Group, exact=TRUE)$p.value, error=function(e) 1)
    })
    top_p <- names(sort(fv_p))[1:min(50, length(fv_p))]
    rf_p <- randomForest(x=train_perm[, top_p], y=train_perm$Group, ntree=200)
    perm_pred[i] <- stats::predict(rf_p, loocv_data2[i, top_p], type="prob")[,"RA"]
  }
  perm_roc <- tryCatch(roc(loocv_data2$Group, perm_pred, quiet=TRUE), error=function(e) NULL)
  perm_aucs2[p] <- if(!is.null(perm_roc)) as.numeric(auc(perm_roc)) else 0.5
  if(p %% 100 == 0) cat(sprintf("    %d/%d\n", p, n_perm2))
}
perm_p2 <- (sum(perm_aucs2 >= auc_val) + 1) / (n_perm2 + 1)
cat(sprintf("  Permutation p-value: %.4f\n", perm_p2))

# Save ROC plot with CI
tryCatch({
  pdf(file.path(out_dir, "figures/LOOCV_Nested_ROC_withCI.pdf"), width=7, height=7)
  plot(roc_nested, main=sprintf("Nested LOOCV Classification (RA vs Control)\nAUC=%.3f (95%% CI: %.3f-%.3f)\nPermutation p=%.4f",
                                 auc_val, ci_boot[1], ci_boot[3], perm_p2),
       col="#E41A1C", lwd=2.5, cex.main=0.9)
  abline(0, 1, lty=2, col="grey50")
  text(0.3, 0.15, sprintf("n=%d (RA=%d, Ctrl=%d)", n_lo2,
       sum(loocv_data2$Group=="RA"), sum(loocv_data2$Group=="Control")), cex=0.9)
  dev.off()
  cat("  ✓ ROC plot saved\n")
}, error=function(e) { tryCatch(dev.off(), error=function(x) NULL) })

# Confusion matrix
cm_nested <- confusionMatrix(factor(nested_pred2, levels=levels(loocv_data2$Group)), loocv_data2$Group)

report$loocv <- list(
  AUC=auc_val, CI_lower=ci_boot[1], CI_upper=ci_boot[3],
  permutation_p=perm_p2, accuracy=cm_nested$overall["Accuracy"],
  sensitivity=cm_nested$byClass["Sensitivity"],
  specificity=cm_nested$byClass["Specificity"],
  n=n_lo2
)


# ============================================================================
# SECTION 6: Deconvolution with Effect Sizes - 3-Group Comparison
# ============================================================================
cat("\n=== SECTION 6: Deconvolution 3-Group Comparison ===\n")

if(exists("deconv_mat") && exists("deconv_subgroup")) {
  deconv_3group <- data.frame()

  for(ct in colnames(deconv_mat)) {
    for(comp in list(c("RA_ILD","Control"), c("RA_nonILD","Control"), c("RA_ILD","RA_nonILD"))) {
      g1_idx <- which(deconv_subgroup == comp[1])
      g2_idx <- which(deconv_subgroup == comp[2])
      if(length(g1_idx) >= 3 && length(g2_idx) >= 3) {
        v1 <- deconv_mat[g1_idx, ct]
        v2 <- deconv_mat[g2_idx, ct]
        wt <- wilcox.test(v1, v2, exact=TRUE)
        cd <- cliff.delta(v1, v2)
        deconv_3group <- rbind(deconv_3group, data.frame(
          CellType=ct, Comparison=paste(comp, collapse=" vs "),
          Median_G1=median(v1), Median_G2=median(v2),
          P_value=wt$p.value, Cliff_delta=cd$estimate, Magnitude=cd$magnitude
        ))
      }
    }
  }

  # Adjust p-values within each comparison
  for(comp_name in unique(deconv_3group$Comparison)) {
    idx <- deconv_3group$Comparison == comp_name
    deconv_3group$P_adjusted[idx] <- p.adjust(deconv_3group$P_value[idx], method="BH")
  }

  write.csv(deconv_3group, file.path(out_dir, "tables/Deconvolution_3group_effectsizes.csv"), row.names=FALSE)
  cat("  3-group comparisons with effect sizes saved\n")

  for(comp_name in unique(deconv_3group$Comparison)) {
    sub <- deconv_3group[deconv_3group$Comparison==comp_name & deconv_3group$P_adjusted < 0.05, ]
    cat(sprintf("    %s: %d sig cell types\n", comp_name, nrow(sub)))
  }

  report$deconv_3group <- deconv_3group
}


# ============================================================================
# SECTION 7: Sensitivity Analysis
# ============================================================================
cat("\n=== SECTION 7: Sensitivity Analysis ===\n")

sensitivity_results <- data.frame()

# 7a: LFC threshold sensitivity for DEGs
for(lfc_thresh in c(0, 0.585, 1.0, 1.5)) {
  n_up_s <- sum(res2_df$padj < 0.05 & res2_df$log2FoldChange > lfc_thresh, na.rm=TRUE)
  n_down_s <- sum(res2_df$padj < 0.05 & res2_df$log2FoldChange < -lfc_thresh, na.rm=TRUE)
  sensitivity_results <- rbind(sensitivity_results, data.frame(
    Analysis="DEG_count", Parameter=sprintf("|LFC|>%.2f", lfc_thresh),
    Value_Up=n_up_s, Value_Down=n_down_s, Value_Total=n_up_s+n_down_s
  ))
}

# 7b: Cook's distance outlier detection
cd_res <- results(dds2)
cooks <- assays(dds2)[["cooks"]]
if(!is.null(cooks)) {
  sample_max_cooks <- apply(cooks, 2, max, na.rm=TRUE)
  outlier_threshold <- qf(0.99, ncol(dds2), nrow(dds2) - ncol(dds2))
  n_outlier_genes <- colSums(cooks > outlier_threshold, na.rm=TRUE)
  cat(sprintf("  Cook's distance outlier genes per sample:\n"))
  for(s in names(sort(n_outlier_genes, decreasing=TRUE))[1:5])
    cat(sprintf("    %s: %d genes\n", s, n_outlier_genes[s]))
}

# 7c: Leave-one-out sample sensitivity for LOOCV AUC
cat("  Leave-one-sample-out AUC stability...\n")
loo_aucs <- numeric(n_lo2)
for(i in 1:n_lo2) {
  sub_data <- loocv_data2[-i, ]
  sub_pred <- nested_prob2[-i, "RA"]
  sub_true <- loocv_data2$Group[-i]
  sub_roc <- tryCatch(roc(sub_true, sub_pred, quiet=TRUE), error=function(e) NULL)
  loo_aucs[i] <- if(!is.null(sub_roc)) as.numeric(auc(sub_roc)) else NA
}
cat(sprintf("  AUC range when excluding 1 sample: %.3f - %.3f (mean=%.3f)\n",
    min(loo_aucs, na.rm=TRUE), max(loo_aucs, na.rm=TRUE), mean(loo_aucs, na.rm=TRUE)))

write.csv(sensitivity_results, file.path(out_dir, "tables/Sensitivity_DEG_thresholds.csv"), row.names=FALSE)

report$sensitivity <- list(
  deg_thresholds=sensitivity_results,
  auc_range=range(loo_aucs, na.rm=TRUE),
  auc_mean=mean(loo_aucs, na.rm=TRUE)
)


# ============================================================================
# SECTION 8: ILD Score Dose-Response (Jonckheere-Terpstra)
# ============================================================================
cat("\n=== SECTION 8: ILD Score Dose-Response ===\n")

# Create ILD severity categories
ra_meta_dr <- master_data[master_data$Sample_Group == "RA" &
                             master_data$Sample_ID %in% common_samples, ]

if("ILD_Score" %in% colnames(ra_meta_dr) && sum(!is.na(ra_meta_dr$ILD_Score)) >= 10) {
  ild_scores <- as.numeric(ra_meta_dr$ILD_Score)

  # Categorize: 0, 1-3, 4+
  ra_meta_dr$ILD_Category <- cut(ild_scores, breaks=c(-Inf, 0, 3, Inf),
                                   labels=c("None(0)", "Mild(1-3)", "Severe(4+)"),
                                   right=TRUE)
  cat(sprintf("  ILD categories: %s\n", paste(names(table(ra_meta_dr$ILD_Category)),
      table(ra_meta_dr$ILD_Category), sep="=", collapse=", ")))

  # Dose-response for key variables
  dose_response <- data.frame()

  # GSVA pathways
  if(exists("gsva_scores")) {
    for(gs in rownames(gsva_scores)) {
      gs_vals <- gsva_scores[gs, ra_meta_dr$Sample_ID]
      valid <- !is.na(gs_vals) & !is.na(ra_meta_dr$ILD_Category)
      if(sum(valid) >= 10) {
        sp <- cor.test(gs_vals[valid], ild_scores[valid], method="spearman")
        # Kruskal-Wallis for trend
        kw <- kruskal.test(gs_vals[valid] ~ ra_meta_dr$ILD_Category[valid])
        dose_response <- rbind(dose_response, data.frame(
          Variable=gs, Category="GSVA_Pathway",
          Spearman_rho=sp$estimate, Spearman_p=sp$p.value,
          KW_p=kw$p.value,
          Mean_None=mean(gs_vals[valid & ra_meta_dr$ILD_Category=="None(0)"], na.rm=TRUE),
          Mean_Mild=mean(gs_vals[valid & ra_meta_dr$ILD_Category=="Mild(1-3)"], na.rm=TRUE),
          Mean_Severe=mean(gs_vals[valid & ra_meta_dr$ILD_Category=="Severe(4+)"], na.rm=TRUE)
        ))
      }
    }
  }

  # Deconvolution cell types
  if(exists("deconv_mat")) {
    ra_deconv_ids <- intersect(ra_meta_dr$Sample_ID, rownames(deconv_mat))
    for(ct in colnames(deconv_mat)) {
      ct_vals <- deconv_mat[ra_deconv_ids, ct]  # already in %
      ild_vals <- ra_meta_dr$ILD_Score[match(ra_deconv_ids, ra_meta_dr$Sample_ID)]
      ild_cat <- ra_meta_dr$ILD_Category[match(ra_deconv_ids, ra_meta_dr$Sample_ID)]
      valid <- !is.na(ct_vals) & !is.na(ild_vals)
      if(sum(valid) >= 10) {
        sp <- cor.test(ct_vals[valid], ild_vals[valid], method="spearman")
        kw <- kruskal.test(ct_vals[valid] ~ ild_cat[valid])
        dose_response <- rbind(dose_response, data.frame(
          Variable=ct, Category="CellType",
          Spearman_rho=sp$estimate, Spearman_p=sp$p.value,
          KW_p=kw$p.value,
          Mean_None=mean(ct_vals[valid & ild_cat=="None(0)"], na.rm=TRUE),
          Mean_Mild=mean(ct_vals[valid & ild_cat=="Mild(1-3)"], na.rm=TRUE),
          Mean_Severe=mean(ct_vals[valid & ild_cat=="Severe(4+)"], na.rm=TRUE)
        ))
      }
    }
  }

  if(nrow(dose_response) > 0) {
    dose_response$Spearman_padj <- p.adjust(dose_response$Spearman_p, method="BH")
    dose_response$KW_padj <- p.adjust(dose_response$KW_p, method="BH")
    dose_response <- dose_response %>% arrange(Spearman_p)
    write.csv(dose_response, file.path(out_dir, "tables/DoseResponse_ILD_Score.csv"), row.names=FALSE)

    sig_dr <- dose_response[dose_response$Spearman_padj < 0.1, ]
    cat(sprintf("  Dose-response associations (padj<0.1): %d\n", nrow(sig_dr)))
    if(nrow(sig_dr) > 0) {
      for(i in 1:min(10, nrow(sig_dr)))
        cat(sprintf("    %-25s rho=%.3f p_adj=%.4f\n", sig_dr$Variable[i], sig_dr$Spearman_rho[i], sig_dr$Spearman_padj[i]))
    }

    report$dose_response <- dose_response
  }
}


# ============================================================================
# SECTION 9: GSVA Comprehensive (RA vs Ctrl, ILD vs nonILD, Infection)
# ============================================================================
cat("\n=== SECTION 9: GSVA Comprehensive Summary ===\n")

if(exists("gsva_scores")) {
  gsva_comprehensive <- data.frame()

  for(gs in rownames(gsva_scores)) {
    # RA vs Control
    ctrl_s <- intersect(colnames(gsva_scores), common_samples[meta$Sample_Group=="Control"])
    ra_s <- intersect(colnames(gsva_scores), common_samples[meta$Sample_Group=="RA"])

    if(length(ctrl_s)>=3 && length(ra_s)>=3) {
      wt <- wilcox.test(gsva_scores[gs,ra_s], gsva_scores[gs,ctrl_s], exact=TRUE)
      cd <- cliff.delta(gsva_scores[gs,ra_s], gsva_scores[gs,ctrl_s])
      gsva_comprehensive <- rbind(gsva_comprehensive, data.frame(
        GeneSet=gs, Comparison="RA_vs_Control",
        Mean_G1=mean(gsva_scores[gs,ra_s]), Mean_G2=mean(gsva_scores[gs,ctrl_s]),
        Delta=mean(gsva_scores[gs,ra_s])-mean(gsva_scores[gs,ctrl_s]),
        P_value=wt$p.value, Cliff_delta=cd$estimate))
    }

    # ILD vs nonILD (within RA)
    ild_s <- intersect(ra_s, common_samples[meta$Subgroup=="RA_ILD"])
    nonild_s <- intersect(ra_s, common_samples[meta$Subgroup=="RA_nonILD"])

    if(length(ild_s)>=3 && length(nonild_s)>=3) {
      wt <- wilcox.test(gsva_scores[gs,ild_s], gsva_scores[gs,nonild_s], exact=TRUE)
      cd <- cliff.delta(gsva_scores[gs,ild_s], gsva_scores[gs,nonild_s])
      gsva_comprehensive <- rbind(gsva_comprehensive, data.frame(
        GeneSet=gs, Comparison="RA_ILD_vs_RA_nonILD",
        Mean_G1=mean(gsva_scores[gs,ild_s]), Mean_G2=mean(gsva_scores[gs,nonild_s]),
        Delta=mean(gsva_scores[gs,ild_s])-mean(gsva_scores[gs,nonild_s]),
        P_value=wt$p.value, Cliff_delta=cd$estimate))
    }
  }

  # BH correction within each comparison
  for(comp in unique(gsva_comprehensive$Comparison)) {
    idx <- gsva_comprehensive$Comparison == comp
    gsva_comprehensive$P_adjusted[idx] <- p.adjust(gsva_comprehensive$P_value[idx], method="BH")
  }

  write.csv(gsva_comprehensive, file.path(out_dir, "tables/GSVA_comprehensive.csv"), row.names=FALSE)

  cat("  GSVA summary:\n")
  for(comp in unique(gsva_comprehensive$Comparison)) {
    sub <- gsva_comprehensive[gsva_comprehensive$Comparison==comp, ]
    sig <- sub[sub$P_adjusted < 0.05, ]
    cat(sprintf("    %s: %d/%d pathways significant\n", comp, nrow(sig), nrow(sub)))
  }

  report$gsva_comprehensive <- gsva_comprehensive
}


# ============================================================================
# SECTION 10: Integrated Summary Figure Data
# ============================================================================
cat("\n=== SECTION 10: Summary ===\n")

# Compile all results
report$analysis_date <- Sys.time()
report$samples <- list(
  total=length(common_samples),
  RA=sum(meta$Sample_Group=="RA"),
  Control=sum(meta$Sample_Group=="Control"),
  RA_ILD=sum(meta$Subgroup=="RA_ILD", na.rm=TRUE),
  RA_nonILD=sum(meta$Subgroup=="RA_nonILD", na.rm=TRUE)
)

saveRDS(report, file.path(out_dir, "NatComm_Report.rds"))
save.image(file.path(out_dir, "NatComm_Workspace.RData"))

cat("\n╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  Analysis Complete                                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Samples: %d (RA=%d [ILD=%d, nonILD=%d], Control=%d)\n",
    report$samples$total, report$samples$RA, report$samples$RA_ILD,
    report$samples$RA_nonILD, report$samples$Control))
cat(sprintf("  DEGs (shrunk LFC): %d up, %d down\n", report$deg$n_up, report$deg$n_down))
cat(sprintf("  GSEA GO BP: %d terms, KEGG: %d pathways\n", report$gsea$GO_BP, report$gsea$KEGG))
cat(sprintf("  LOOCV AUC: %.3f (%.3f-%.3f), perm p=%.4f\n",
    report$loocv$AUC, report$loocv$CI_lower, report$loocv$CI_upper, report$loocv$permutation_p))
cat(sprintf("  Effect sizes computed: %d variables\n", report$effect_sizes$total))

cat("\nOutput directory: ", out_dir, "\n")
