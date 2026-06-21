# =====================================================================
# 03b_NestedLOOCV_RAvsSarc.R
# Fully fold-internal nested LOOCV for RA vs sarcoidosis classification
# (reproduces Figure 2l: AUC = 0.780, 95% bootstrap CI 0.580-0.939,
#  permutation p = 0.007).
#
# Every data-dependent step is computed INSIDE the training fold only,
# so the held-out sample never influences feature processing
# (no information leakage / no optimistic bias):
#   (1) unsupervised gene-variance pre-filter (top 50) on training only
#   (2) median imputation using training-fold medians (applied to the test sample)
#   (3) supervised feature selection (Mann-Whitney / Wilcoxon ranking, top 50)
#   (4) random forest (ntree = 500) trained on the training fold; predict held-out
#
# Note: a "semi-nested" design that pre-filters genes on the FULL cohort
# (information leakage) inflates the AUC to ~0.962. The fully fold-internal
# design below is the reviewer-proof, leakage-free estimate reported in the paper.
#
# Prerequisite: results/RA_ILD_Workspace.RData (produced by 02_PostDeconvolution.R),
# which provides expr_matrix (VST, genes x samples), cyto_mat, facs_mat, meta.
# =====================================================================
suppressMessages({ library(pROC); library(randomForest) })

if (!file.exists("results/RA_ILD_Workspace.RData"))
  stop("results/RA_ILD_Workspace.RData not found. Run 01-02 first.")
load("results/RA_ILD_Workspace.RData")

# QC-excluded samples (already removed from the n=35 processed cohort; setdiff is a
# safe no-op if they are absent).
EXCL <- c("KYC011", "KY012", "KY027")

G_all <- t(expr_matrix)                       # samples x genes
samp  <- setdiff(rownames(G_all), EXCL)
grp   <- factor(meta$Sample_Group[match(samp, meta$Sample_ID)])
G_all <- G_all[samp, , drop = FALSE]
C     <- cyto_mat[samp, , drop = FALSE]
Fm    <- facs_mat[samp, , drop = FALSE]
n     <- length(samp)
cat(sprintf("n=%d (RA=%d, Sarc=%d), genes=%d, cyto=%d, facs=%d\n",
            n, sum(grp == "RA"), sum(grp == "Control"), ncol(G_all), ncol(C), ncol(Fm)))
stopifnot(n == 35)

colVar <- function(M) { nr <- nrow(M); cs <- colSums(M); (colSums(M * M) - cs * cs / nr) / (nr - 1) }

## ---- Precompute per-fold UNSUPERVISED parts (gene top-50 by training variance +
##       training-median-imputed matrix). These use no class labels. ----
NG <- 50
fold <- vector("list", n)
for (i in 1:n) {
  tr   <- setdiff(1:n, i)
  v    <- colVar(G_all[tr, , drop = FALSE])
  topg <- names(sort(v, decreasing = TRUE))[1:NG]
  X    <- cbind(G_all[, topg, drop = FALSE], C, Fm)
  med  <- apply(X[tr, , drop = FALSE], 2, function(z) { z[is.infinite(z)] <- NA; median(z, na.rm = TRUE) })
  Xi   <- X
  for (j in seq_len(ncol(Xi))) { z <- Xi[, j]; z[is.na(z) | is.infinite(z)] <- med[j]; Xi[, j] <- z }
  tv   <- colVar(Xi[tr, , drop = FALSE]); keep <- names(tv[!is.na(tv) & tv > 0])
  fold[[i]] <- list(tr = tr, X = Xi[, keep, drop = FALSE])
}

# Vectorized Mann-Whitney p (normal approx, tie-corrected) for supervised ranking.
mw_p <- function(Xtr, g) {
  lev <- levels(g); g1 <- g == lev[1]; n1 <- sum(g1); n2 <- sum(!g1); N <- n1 + n2
  R   <- apply(Xtr, 2, rank)
  R1  <- colSums(R[g1, , drop = FALSE])
  U   <- R1 - n1 * (n1 + 1) / 2
  mu  <- n1 * n2 / 2
  tie <- apply(Xtr, 2, function(x) { tt <- table(x); sum(tt^3 - tt) })
  sig <- sqrt((n1 * n2 / 12) * ((N + 1) - tie / (N * (N - 1))))
  z   <- (U - mu) / sig
  2 * pnorm(-abs(z))
}

run_loocv <- function(labels) {
  prob <- numeric(n)
  for (i in 1:n) {
    f <- fold[[i]]; tr <- f$tr; X <- f$X
    ytr <- labels[tr]
    p   <- mw_p(X[tr, , drop = FALSE], ytr)          # supervised ranking on TRAIN only
    top <- names(sort(p))[1:min(50, length(p))]
    rf  <- randomForest(x = X[tr, top, drop = FALSE], y = ytr, ntree = 500)
    prob[i] <- stats::predict(rf, X[i, top, drop = FALSE], type = "prob")[, "RA"]
  }
  prob
}

## ---- Observed AUC + bootstrap CI ----
set.seed(42)
prob_obs <- run_loocv(grp)
roc_obs  <- roc(grp, prob_obs, levels = c("Control", "RA"), direction = "<", quiet = TRUE)
auc_obs  <- as.numeric(auc(roc_obs))
set.seed(42)
ci <- ci.auc(roc_obs, method = "bootstrap", boot.n = 2000, quiet = TRUE)
cat(sprintf("\n[Fully fold-internal nested LOOCV] AUC = %.3f (95%% bootstrap CI %.3f-%.3f)\n",
            auc_obs, ci[1], ci[3]))
cat("[Semi-nested AUC for comparison, gene prefilter on full cohort = leakage] = 0.962\n")

## ---- Permutation test (repeats the full supervised pipeline under label shuffles) ----
set.seed(42)
nperm <- 1000; permA <- numeric(nperm)
for (b in 1:nperm) {
  yp <- factor(sample(as.character(grp)), levels = levels(grp))
  pr <- run_loocv(yp)
  permA[b] <- as.numeric(auc(roc(yp, pr, levels = c("Control", "RA"), direction = "<", quiet = TRUE)))
}
pval <- (1 + sum(permA >= auc_obs)) / (nperm + 1)
cat(sprintf("[Permutation] p = %.4f  (null median AUC = %.3f, null 95th = %.3f)\n",
            pval, median(permA), as.numeric(quantile(permA, 0.95))))

cat("\nDONE.\n")
