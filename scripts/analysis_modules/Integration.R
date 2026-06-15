# ============================================================================
# Multi-omics Integration Analysis (MOFA2代替)
# 1. Cross-layer correlation network
# 2. Multi-omics integrated heatmap
# 3. Multi-omics factor analysis (PCA-based surrogate)
# 4. Layer contribution to classification
# ============================================================================
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)

load("results/CT_Multiomics_Results.RData")

suppressPackageStartupMessages({
  library(tidyverse); library(pheatmap); library(ggplot2); library(gridExtra)
  library(pROC); library(randomForest); library(corrplot)
})
set.seed(42)
out_dir <- "results"
fig_dir <- "results/figures_final"

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  Multi-omics Integration Analysis                        ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

tn <- function(bs=7.5) {
  theme_classic(base_size=bs) %+replace%
    theme(text=element_text(color="black"),
          axis.text=element_text(size=rel(.9),color="black"),
          axis.title=element_text(size=rel(.95),face="bold"),
          plot.title=element_text(size=rel(1.05),face="bold",hjust=0),
          legend.text=element_text(size=rel(.78)),
          legend.key.size=unit(2.5,"mm"),
          legend.background=element_rect(fill=alpha("white",.92),color=NA),
          panel.border=element_rect(fill=NA,color="black",linewidth=.4),
          axis.line=element_blank(),
          plot.margin=unit(c(3,4,3,3),"mm"))
}

samples <- common_samples
group <- ifelse(grepl("^(KYC|Sarcoidosis)", samples), "Sarcoidosis", "RA")

# ============================================================================
# Prepare representative features per layer (prefix preserved)
# ============================================================================
cat("=== Data Preparation ===\n")

# Expression: top 20 variable genes
ev <- apply(expr_mat[samples,], 2, var, na.rm=TRUE)
X_expr <- expr_mat[samples, names(sort(ev, decreasing=TRUE))[1:20]]

# BALF Cytokines (from master_data)
balf_cols <- grep("^BALF_", colnames(master_data), value=TRUE)
balf_cols <- balf_cols[!grepl("Treg|Th1|Th2|Th17|Macro|Neutro|Lympho|Eosin|Baso|Plasma_Cell|CD[0-9]|Percent|Activated|Ratio|CD14|CD86|CD66", balf_cols)]
X_balf <- master_data[match(samples, master_data$Sample_ID), balf_cols]
for(col in colnames(X_balf)) X_balf[[col]] <- as.numeric(X_balf[[col]])
X_balf <- log2(as.matrix(X_balf)+1); rownames(X_balf) <- samples
# Top 15 by variance
bv <- apply(X_balf, 2, var, na.rm=TRUE)
X_balf <- X_balf[, names(sort(bv, decreasing=TRUE))[1:min(15, sum(!is.na(bv)&bv>0))]]

# Serum Cytokines
serum_cols <- grep("^Serum_", colnames(master_data), value=TRUE)
X_serum <- master_data[match(samples, master_data$Sample_ID), serum_cols]
for(col in colnames(X_serum)) X_serum[[col]] <- as.numeric(X_serum[[col]])
X_serum <- log2(as.matrix(X_serum)+1); rownames(X_serum) <- samples
sv <- apply(X_serum, 2, var, na.rm=TRUE)
X_serum <- X_serum[, names(sort(sv, decreasing=TRUE))[1:min(15, sum(!is.na(sv)&sv>0))]]

# FACS (all)
X_facs <- facs_mat[samples, ]

# Deconvolution
X_deconv <- deconv_mat[samples, ]

# GSVA
X_gsva <- t(gsva_scores[, samples])

cat(sprintf("  Expression: %d, BALF Cyto: %d, Serum Cyto: %d, FACS: %d, Deconv: %d, GSVA: %d\n",
    ncol(X_expr), ncol(X_balf), ncol(X_serum), ncol(X_facs), ncol(X_deconv), ncol(X_gsva)))

# ============================================================================
# 1. Cross-layer correlation network
# ============================================================================
cat("\n=== 1. Cross-layer Correlation Network ===\n")

# Combine all layers
all_mat <- cbind(X_expr, X_balf, X_serum, X_facs, X_deconv, X_gsva)
all_mat[is.na(all_mat) | !is.finite(all_mat)] <- 0
layer_labels <- c(rep("Expression", ncol(X_expr)),
                   rep("BALF_Cytokine", ncol(X_balf)),
                   rep("Serum_Cytokine", ncol(X_serum)),
                   rep("FCM", ncol(X_facs)),
                   rep("Deconvolution", ncol(X_deconv)),
                   rep("GSVA", ncol(X_gsva)))

# Spearman correlation matrix
cor_mat <- cor(all_mat, method="spearman", use="pairwise.complete.obs")

# Cross-layer edges (|rho| > 0.5)
edges <- data.frame()
for(i in 1:(ncol(all_mat)-1)) {
  for(j in (i+1):ncol(all_mat)) {
    if(layer_labels[i] != layer_labels[j] && abs(cor_mat[i,j]) > 0.5 && !is.na(cor_mat[i,j])) {
      edges <- rbind(edges, data.frame(
        Node1=colnames(all_mat)[i], Layer1=layer_labels[i],
        Node2=colnames(all_mat)[j], Layer2=layer_labels[j],
        rho=cor_mat[i,j]))
    }
  }
}
edges <- edges %>% arrange(desc(abs(rho)))
write.csv(edges, file.path(out_dir, "tables/CrossLayer_Network_Edges.csv"), row.names=FALSE)

cat(sprintf("  Cross-layer edges (|rho|>0.5): %d\n", nrow(edges)))
cat("  By layer pair:\n")
edge_summary <- edges %>% mutate(Pair=paste(pmin(Layer1,Layer2), pmax(Layer1,Layer2), sep=" - ")) %>%
  group_by(Pair) %>% summarise(N=n(), Mean_abs_rho=mean(abs(rho)), .groups="drop") %>% arrange(desc(N))
for(i in 1:nrow(edge_summary))
  cat(sprintf("    %-40s %3d edges (mean|rho|=%.3f)\n",
      edge_summary$Pair[i], edge_summary$N[i], edge_summary$Mean_abs_rho[i]))

# ============================================================================
# 2. Multi-omics integrated heatmap
# ============================================================================
cat("\n=== 2. Integrated Heatmap ===\n")

# Scale each layer independently
all_scaled <- all_mat
for(layer in unique(layer_labels)) {
  idx <- which(layer_labels == layer)
  all_scaled[, idx] <- scale(all_mat[, idx])
}
all_scaled[!is.finite(all_scaled)] <- 0

# Annotation
ann_col <- data.frame(
  Group = ifelse(grepl("^(KYC|Sarcoidosis)", samples), "Sarcoidosis", "RA"),
  row.names = samples)
ann_row <- data.frame(Layer = layer_labels, row.names = colnames(all_mat))
ann_colors <- list(
  Group = c("RA"="#C0392B", "Sarcoidosis"="#2980B9"),
  Layer = c(Expression="#E74C3C", BALF_Cytokine="#3498DB", Serum_Cytokine="#27AE60",
            FACS="#9B59B6", Deconvolution="#E67E22", GSVA="#F39C12"))

# Order samples by group
sample_order <- samples[order(group)]

tryCatch({
  pdf(file.path(fig_dir, "FigIntegration_Heatmap.pdf"), width=14, height=10)
  pheatmap(t(all_scaled[sample_order, ]),
           cluster_rows=TRUE, cluster_cols=FALSE,
           annotation_col=ann_col[sample_order, , drop=FALSE],
           annotation_row=ann_row,
           annotation_colors=ann_colors,
           show_colnames=TRUE, show_rownames=FALSE,
           color=colorRampPalette(c("#2980B9","white","#C0392B"))(100),
           breaks=seq(-3, 3, length.out=101),
           fontsize_col=6, fontsize_row=4,
           main="Multi-omics Integrated Heatmap (z-scored within each layer)")
  dev.off()
  cat("  Heatmap done\n")
}, error=function(e) { tryCatch(dev.off(), error=function(x) NULL); cat(sprintf("  ✗ %s\n",e$message)) })

# ============================================================================
# 3. Joint PCA (pseudo-MOFA: PCA on concatenated scaled data)
# ============================================================================
cat("\n=== 3. Joint PCA (Multi-omics Factor Analysis) ===\n")

jpca <- prcomp(all_scaled, scale.=FALSE)
jve <- summary(jpca)$importance[2,1:10]*100
jfactors <- jpca$x[, 1:5]

cat("  Variance explained: ")
cat(paste(sprintf("JPC%d=%.1f%%", 1:5, jve[1:5]), collapse=", "), "\n")

# Factor vs disease
cat("\n  Joint factor vs Disease:\n")
for(f in 1:5) {
  wt <- wilcox.test(jfactors[group=="RA", f], jfactors[group=="Sarcoidosis", f], exact=TRUE)
  cat(sprintf("    JPC%d: RA=%+.2f Sarc=%+.2f p=%.4f %s\n",
      f, mean(jfactors[group=="RA",f]), mean(jfactors[group=="Sarcoidosis",f]),
      wt$p.value, ifelse(wt$p.value<0.05,"*","")))
}

# Factor vs infection
inf <- master_data$respiratory_infection[match(samples, master_data$Sample_ID)]
cat("\n  Joint factor vs Infection (RA only):\n")
for(f in 1:5) {
  ra_mask <- group=="RA" & !is.na(inf)
  pos <- jfactors[ra_mask & inf==1, f]; neg <- jfactors[ra_mask & inf==0, f]
  if(length(pos)>=3 && length(neg)>=3) {
    wt <- wilcox.test(pos, neg, exact=TRUE)
    cat(sprintf("    JPC%d: Inf+=%+.2f Inf-=%+.2f p=%.4f %s\n",
        f, mean(pos), mean(neg), wt$p.value, ifelse(wt$p.value<0.05,"*","")))
  }
}

# Loading: which layers contribute to each JPC
cat("\n  Layer contributions to joint factors:\n")
loadings <- jpca$rotation[, 1:5]
layer_contrib <- data.frame()
for(f in 1:5) {
  for(layer in unique(layer_labels)) {
    idx <- which(layer_labels == layer)
    contrib <- sum(loadings[idx, f]^2) / sum(loadings[, f]^2) * 100
    layer_contrib <- rbind(layer_contrib, data.frame(Factor=f, Layer=layer, Contribution=contrib))
  }
}
layer_wide <- layer_contrib %>% pivot_wider(names_from=Factor, values_from=Contribution, names_prefix="JPC")
cat("  "); print(as.data.frame(layer_wide), row.names=FALSE)
write.csv(layer_contrib, file.path(out_dir, "tables/JointPCA_LayerContributions.csv"), row.names=FALSE)

# ============================================================================
# 4. Layer contribution to classification (ablation study)
# ============================================================================
cat("\n=== 4. Layer Ablation Study ===\n")

Y <- factor(group, levels=c("Sarcoidosis","RA"))

# Full model
run_loocv <- function(X, Y, label) {
  set.seed(42)
  n <- length(Y); pred <- numeric(n)
  X[is.na(X)|!is.finite(X)] <- 0
  fv <- apply(X, 2, var, na.rm=TRUE); X <- X[, !is.na(fv) & fv > 0, drop=FALSE]
  for(i in 1:n) {
    tr_X <- X[-i,,drop=F]; tr_Y <- Y[-i]
    pv <- apply(tr_X, 2, function(x) tryCatch(wilcox.test(x~tr_Y, exact=TRUE)$p.value,error=function(e)1))
    tp <- names(sort(pv))[1:min(30,length(pv))]
    rf <- randomForest(x=tr_X[,tp,drop=F], y=tr_Y, ntree=500)
    pred[i] <- stats::predict(rf, X[i,tp,drop=F], type="prob")[,"RA"]
  }
  roc_obj <- tryCatch(roc(Y, pred, quiet=TRUE), error=function(e) NULL)
  auc_val <- if(!is.null(roc_obj)) as.numeric(auc(roc_obj)) else NA
  cat(sprintf("    %-30s AUC=%.3f (p=%d)\n", label, auc_val, ncol(X)))
  return(auc_val)
}

cat("  Single-layer LOOCV:\n")
auc_results <- data.frame()
for(layer_info in list(
  list("Expression", X_expr), list("BALF_Cytokine", X_balf), list("Serum_Cytokine", X_serum),
  list("FCM", X_facs), list("Deconvolution", X_deconv), list("GSVA", X_gsva))) {
  auc_val <- run_loocv(as.matrix(layer_info[[2]]), Y, layer_info[[1]])
  auc_results <- rbind(auc_results, data.frame(Model=layer_info[[1]], AUC=auc_val, Type="Single"))
}

# All combined
cat("\n  Combined models:\n")
auc_all <- run_loocv(all_mat, Y, "All layers")
auc_results <- rbind(auc_results, data.frame(Model="All_combined", AUC=auc_all, Type="Combined"))

# Leave-one-layer-out
cat("\n  Leave-one-layer-out:\n")
for(layer_name in unique(layer_labels)) {
  idx_keep <- which(layer_labels != layer_name)
  auc_loo <- run_loocv(all_mat[, idx_keep], Y, sprintf("Without %s", layer_name))
  auc_results <- rbind(auc_results, data.frame(
    Model=sprintf("Without_%s", layer_name), AUC=auc_loo, Type="Ablation"))
}

write.csv(auc_results, file.path(out_dir, "tables/LayerAblation_AUC.csv"), row.names=FALSE)

# ============================================================================
# 5. Figures
# ============================================================================
cat("\n=== 5. Figures ===\n")

CG <- c("RA"="#C0392B","Sarcoidosis"="#2980B9")

# a: Joint PCA
df_jpca <- data.frame(JPC1=jfactors[,1], JPC2=jfactors[,2], Group=group)
pa <- ggplot(df_jpca, aes(JPC1, JPC2, color=Group)) +
  geom_point(size=2.5, alpha=.8) + stat_ellipse(level=.95, linewidth=.5) +
  scale_color_manual(values=CG) +
  labs(x=sprintf("JPC1 (%.1f%%)", jve[1]), y=sprintf("JPC2 (%.1f%%)", jve[2]),
       title="a  Joint multi-omics PCA", color="") +
  tn() + theme(legend.position="bottom", legend.direction="horizontal")

# b: Layer contributions
lc_df <- layer_contrib %>% filter(Factor <= 3) %>%
  mutate(Factor=sprintf("JPC%d", Factor))
pb <- ggplot(lc_df, aes(Factor, Contribution, fill=Layer)) +
  geom_col(width=.6) +
  scale_fill_manual(values=c(Expression="#E74C3C", BALF_Cytokine="#3498DB",
                              Serum_Cytokine="#27AE60", FACS="#9B59B6",
                              Deconvolution="#E67E22", GSVA="#F39C12")) +
  labs(x="", y="Contribution (%)", title="b  Layer contributions", fill="") +
  tn(7) + theme(legend.position="bottom", legend.text=element_text(size=5.5))

# c: AUC comparison (single + combined + ablation)
auc_results$Model <- fct_reorder(auc_results$Model, auc_results$AUC)
pc <- ggplot(auc_results, aes(AUC, Model, fill=Type)) +
  geom_col(width=.55) +
  scale_fill_manual(values=c(Single="#3498DB", Combined="#C0392B", Ablation="#E67E22")) +
  geom_vline(xintercept=.5, linetype="dashed", color="grey60", linewidth=.2) +
  labs(x="LOOCV AUC", y="", title="c  Layer ablation study", fill="") +
  xlim(0, 1) +
  tn(7) + theme(legend.position="bottom", axis.text.y=element_text(size=6))

# d: Cross-layer network summary
edge_sum2 <- edges %>%
  mutate(Pair=paste(pmin(Layer1,Layer2), "\n", pmax(Layer1,Layer2))) %>%
  group_by(Pair) %>% summarise(N=n(), .groups="drop") %>% arrange(desc(N))
edge_sum2$Pair <- fct_reorder(edge_sum2$Pair, edge_sum2$N)
pd <- ggplot(edge_sum2, aes(N, Pair)) +
  geom_col(fill="#34495E", width=.55) +
  labs(x="Cross-layer edges (|rho|>0.5)", y="", title="d  Layer connectivity") +
  tn(7) + theme(axis.text.y=element_text(size=6))

ggsave(file.path(fig_dir, "FigIntegration_Summary.pdf"),
       arrangeGrob(pa, pb, pc, pd, ncol=2),
       width=200, height=190, units="mm", dpi=300)
cat("  FigIntegration_Summary.pdf done\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n╔═══════════════════════════════════════════════════════════╗\n")
cat("║  Integration Summary                                      ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Cross-layer edges: %d\n", nrow(edges)))
cat(sprintf("  Joint PCA: JPC1=%.1f%%, JPC2=%.1f%%\n", jve[1], jve[2]))
cat(sprintf("  All-layer AUC: %.3f\n", auc_all))
cat(sprintf("  Best single layer: %s (AUC=%.3f)\n",
    auc_results$Model[auc_results$Type=="Single"][which.max(auc_results$AUC[auc_results$Type=="Single"])],
    max(auc_results$AUC[auc_results$Type=="Single"])))

save(edges, jfactors, jve, layer_contrib, auc_results,
     file=file.path(out_dir, "Integration_Results.RData"))
cat("  saved\n")
