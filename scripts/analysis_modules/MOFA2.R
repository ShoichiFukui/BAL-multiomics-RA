# ============================================================================
# MOFA2 Multi-omics Factor Analysis
# Unsupervised integration — identifies shared variation across omics layers
# ============================================================================
library(reticulate)
use_python("/tmp/mofa_venv/bin/python3", required = TRUE)
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)


suppressPackageStartupMessages({
  library(tidyverse); library(MOFA2); library(ggplot2); library(gridExtra)
  library(pROC); library(effsize)
})
set.seed(42)
out_dir <- "results"
fig_dir <- "results/figures_final"

cat("╔═══════════════════════════════════════════════════════════╗\n")
cat("║  MOFA2 Multi-omics Factor Analysis                       ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# 1. Data preparation (same blocks as DIABLO, from master_data)
# ============================================================================
cat("=== 1. Data Preparation ===\n")
samples <- common_samples

# Block 1: Expression (top 100 variable genes)
expr_vars <- apply(expr_mat[samples,], 2, var, na.rm=TRUE)
top100_genes <- names(sort(expr_vars, decreasing=TRUE))[1:100]
X_expr <- expr_mat[samples, top100_genes]
X_expr[is.na(X_expr)] <- 0

# Block 2: BALF Cytokines
balf_cyto_cols <- grep("^BALF_", colnames(master_data), value=TRUE)
balf_cyto_cols <- balf_cyto_cols[!grepl("Treg|Th1|Th2|Th17|Macro|Neutro|Lympho|Eosin|Baso|Plasma_Cell|CD[0-9]|Percent|Activated|Ratio|CD14|CD86|CD66", balf_cyto_cols)]
X_balf <- as.data.frame(master_data[match(samples, master_data$Sample_ID), balf_cyto_cols])
for(col in colnames(X_balf)) X_balf[[col]] <- as.numeric(X_balf[[col]])
X_balf <- log2(as.matrix(X_balf) + 1); rownames(X_balf) <- samples
X_balf[!is.finite(X_balf)] <- NA

# Block 3: Serum Cytokines
serum_cyto_cols <- grep("^Serum_", colnames(master_data), value=TRUE)
X_serum <- as.data.frame(master_data[match(samples, master_data$Sample_ID), serum_cyto_cols])
for(col in colnames(X_serum)) X_serum[[col]] <- as.numeric(X_serum[[col]])
X_serum <- log2(as.matrix(X_serum) + 1); rownames(X_serum) <- samples
X_serum[!is.finite(X_serum)] <- NA

# Block 4: FACS
X_facs <- facs_mat[samples, ]

# Block 5: Cellular (Deconv + GSVA)
X_deconv <- deconv_mat[samples, ]
X_gsva <- t(gsva_scores[, samples])
X_cellular <- cbind(X_deconv, X_gsva)
colnames(X_cellular)[1:ncol(X_deconv)] <- paste0("Deconv_", colnames(X_deconv))
colnames(X_cellular)[(ncol(X_deconv)+1):ncol(X_cellular)] <- paste0("GSVA_", colnames(X_gsva))

# Remove zero-variance
for(nm in c("X_expr","X_balf","X_serum","X_facs","X_cellular")) {
  X <- get(nm); fv <- apply(X, 2, var, na.rm=TRUE)
  assign(nm, X[, !is.na(fv) & fv > 1e-10, drop=FALSE])
}

cat(sprintf("  Expression: %d x %d\n", nrow(X_expr), ncol(X_expr)))
cat(sprintf("  BALF Cytokines: %d x %d\n", nrow(X_balf), ncol(X_balf)))
cat(sprintf("  Serum Cytokines: %d x %d\n", nrow(X_serum), ncol(X_serum)))
cat(sprintf("  FACS: %d x %d\n", nrow(X_facs), ncol(X_facs)))
cat(sprintf("  Cellular: %d x %d\n", nrow(X_cellular), ncol(X_cellular)))

# ============================================================================
# 2. Create MOFA object
# ============================================================================
cat("\n=== 2. Creating MOFA Object ===\n")

data_list <- list(
  Expression = t(X_expr),
  BALF_Cytokine = t(X_balf),
  Serum_Cytokine = t(X_serum),
  FACS = t(X_facs),
  Cellular = t(X_cellular)
)

mofa <- create_mofa(data_list)

# Model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 5  # n=35に対して保守的に

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "medium"
train_opts$maxiter <- 500
train_opts$seed <- 42
train_opts$verbose <- FALSE

data_opts <- get_default_data_options(mofa)

mofa <- prepare_mofa(mofa, model_options=model_opts,
                      training_options=train_opts, data_options=data_opts)

cat("  Running MOFA2...\n")
mofa_trained <- run_mofa(mofa, use_basilisk=FALSE)
cat("  done\n")

# ============================================================================
# 3. Results
# ============================================================================
cat("\n=== 3. Results ===\n")

# Variance explained
r2 <- get_variance_explained(mofa_trained)
cat("  Variance explained per factor per view:\n")
print(round(r2$r2_per_factor[[1]], 2))
cat("\n  Total variance explained per view:\n")
print(round(r2$r2_total[[1]], 2))

# Factor values
factors <- get_factors(mofa_trained)[[1]]
cat(sprintf("\n  Factors extracted: %d\n", ncol(factors)))

# Associate factors with disease group
group <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors)), "Sarcoidosis", "RA")
cat("\n  Factor vs Disease group (Wilcoxon):\n")
factor_disease <- data.frame()
for(f in 1:ncol(factors)) {
  wt <- wilcox.test(factors[group=="RA", f], factors[group=="Sarcoidosis", f], exact=TRUE)
  cat(sprintf("    Factor %d: RA mean=%+.3f, Sarc mean=%+.3f, p=%.4f %s\n",
      f, mean(factors[group=="RA",f]), mean(factors[group=="Sarcoidosis",f]),
      wt$p.value, ifelse(wt$p.value<0.05,"*","")))
  factor_disease <- rbind(factor_disease, data.frame(
    Factor=f, p_disease=wt$p.value,
    RA_mean=mean(factors[group=="RA",f]),
    Sarc_mean=mean(factors[group=="Sarcoidosis",f])))
}

# Associate factors with infection
inf_status <- master_data$respiratory_infection[match(rownames(factors), master_data$Sample_ID)]
cat("\n  Factor vs Infection (Wilcoxon, RA only):\n")
factor_infection <- data.frame()
for(f in 1:ncol(factors)) {
  ra_mask <- group=="RA" & !is.na(inf_status)
  if(sum(ra_mask) >= 10) {
    pos <- factors[ra_mask & inf_status==1, f]
    neg <- factors[ra_mask & inf_status==0, f]
    if(length(pos)>=3 && length(neg)>=3) {
      wt <- wilcox.test(pos, neg, exact=TRUE)
      cat(sprintf("    Factor %d: Inf+ mean=%+.3f, Inf- mean=%+.3f, p=%.4f %s\n",
          f, mean(pos), mean(neg), wt$p.value, ifelse(wt$p.value<0.05,"*","")))
      factor_infection <- rbind(factor_infection, data.frame(
        Factor=f, p_infection=wt$p.value, InfPos_mean=mean(pos), InfNeg_mean=mean(neg)))
    }
  }
}

# ============================================================================
# 4. Figures
# ============================================================================
cat("\n=== 4. Figures ===\n")

tn <- function(bs=8) {
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
CG <- c("RA"="#C0392B","Sarcoidosis"="#2980B9")
CI <- c("No infection"="#3498DB","Future infection"="#E74C3C")

# a: Variance explained heatmap
r2_mat <- r2$r2_per_factor[[1]]
r2_long <- as.data.frame(r2_mat) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to="View", values_to="R2")
r2_long$Factor <- factor(r2_long$Factor, levels=paste0("Factor",1:nrow(r2_mat)))
pa <- ggplot(r2_long, aes(View, Factor, fill=R2)) + geom_tile(color="white") +
  scale_fill_gradient(low="white", high="#C0392B", name="R2 (%)") +
  geom_text(aes(label=sprintf("%.1f",R2)), size=2) +
  labs(x="", y="", title="a  Variance explained") +
  tn(7) + theme(axis.text.x=element_text(angle=35, hjust=1, size=6))

# b: Factor 1 vs Factor 2 (disease)
df_fac <- data.frame(F1=factors[,1], F2=factors[,2], Group=group)
pb <- ggplot(df_fac, aes(F1, F2, color=Group)) + geom_point(size=2.5, alpha=.8) +
  scale_color_manual(values=CG) +
  labs(x="Factor 1", y="Factor 2", title="b  Disease", color="") +
  tn() + theme(legend.position="bottom", legend.direction="horizontal")

# c: Factor vs disease (significant factors)
sig_factors <- factor_disease$Factor[factor_disease$p_disease < 0.1]
if(length(sig_factors) == 0) sig_factors <- 1:min(2, ncol(factors))
fac_box <- data.frame()
for(f in sig_factors) {
  fac_box <- rbind(fac_box, data.frame(Factor=paste0("Factor ",f), Value=factors[,f], Group=group))
}
pc <- ggplot(fac_box, aes(Group, Value, fill=Group)) +
  geom_boxplot(outlier.shape=NA, width=.4, linewidth=.2) +
  geom_jitter(width=.07, size=1, alpha=.5, show.legend=FALSE) +
  scale_fill_manual(values=CG) + facet_wrap(~Factor, scales="free_y") +
  labs(x="", y="Factor value", title="c  Factors by disease") +
  tn(7) + theme(legend.position="none")

# d: Factor vs infection
if(nrow(factor_infection) > 0) {
  sig_inf <- factor_infection$Factor[factor_infection$p_infection < 0.1]
  if(length(sig_inf) == 0) sig_inf <- 1:min(2, ncol(factors))
  fac_inf <- data.frame()
  for(f in sig_inf) {
    ra_mask2 <- group=="RA" & !is.na(inf_status)
    inf_label <- ifelse(inf_status[ra_mask2]==1, "Future infection", "No infection")
    fac_inf <- rbind(fac_inf, data.frame(Factor=paste0("Factor ",f),
                                          Value=factors[ra_mask2,f], InfGroup=inf_label))
  }
  pd <- ggplot(fac_inf, aes(InfGroup, Value, fill=InfGroup)) +
    geom_boxplot(outlier.shape=NA, width=.4, linewidth=.2) +
    geom_jitter(width=.07, size=1, alpha=.5, show.legend=FALSE) +
    scale_fill_manual(values=CI) + facet_wrap(~Factor, scales="free_y") +
    labs(x="", y="Factor value", title="d  Factors by infection") +
    tn(7) + theme(legend.position="none", axis.text.x=element_text(size=6))
} else {
  pd <- ggplot() + geom_blank() + labs(title="d  No significant factors") + tn()
}

ggsave(file.path(fig_dir, "FigMOFA2_Results.pdf"),
       arrangeGrob(pa, pb, pc, pd, ncol=2),
       width=200, height=180, units="mm", dpi=300)
cat("  FigMOFA2_Results.pdf done\n")

# ============================================================================
# 5. Top weights per factor
# ============================================================================
cat("\n=== 5. Top Weights ===\n")
weights <- get_weights(mofa_trained)
all_weights <- data.frame()
for(view in names(weights)) {
  w <- weights[[view]]
  for(f in 1:ncol(w)) {
    top_idx <- order(abs(w[,f]), decreasing=TRUE)[1:min(10, nrow(w))]
    all_weights <- rbind(all_weights, data.frame(
      View=view, Factor=f, Feature=rownames(w)[top_idx], Weight=w[top_idx,f]))
  }
}
write.csv(all_weights, file.path(out_dir, "tables/MOFA2_top_weights.csv"), row.names=FALSE)

# Show top weights for significant factor
best_factor <- if(nrow(factor_disease)>0) factor_disease$Factor[which.min(factor_disease$p_disease)] else 1
cat(sprintf("\n  Top weights for Factor %d (most disease-associated):\n", best_factor))
for(view in names(weights)) {
  w <- weights[[view]][, best_factor]
  top5 <- names(sort(abs(w), decreasing=TRUE))[1:5]
  cat(sprintf("  [%s]\n", view))
  for(f in top5) cat(sprintf("    %-35s %+.4f\n", f, w[f]))
}

# ============================================================================
# 6. Summary
# ============================================================================
cat("\n╔═══════════════════════════════════════════════════════════╗\n")
cat("║  MOFA2 Summary                                            ║\n")
cat("╚═══════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Views: %d\n", length(data_list)))
cat(sprintf("  Samples: %d\n", nrow(factors)))
cat(sprintf("  Factors: %d\n", ncol(factors)))
cat(sprintf("  Total R2 per view:\n"))
for(v in names(r2$r2_total[[1]])) cat(sprintf("    %-20s: %.1f%%\n", v, r2$r2_total[[1]][v]))
cat(sprintf("  Disease-associated factors (p<0.05): %d\n", sum(factor_disease$p_disease<0.05)))
if(nrow(factor_infection)>0)
  cat(sprintf("  Infection-associated factors (p<0.05): %d\n", sum(factor_infection$p_infection<0.05)))

save(mofa_trained, factors, r2, factor_disease, factor_infection, all_weights,
     file=file.path(out_dir, "MOFA2_Results.RData"))
cat("  Results saved\n")
