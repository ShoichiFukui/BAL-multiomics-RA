# 04_disease_MOFA2.R -- disease MOFA2, six views (Fig. 5a-c; Suppl. Data 8).
# Produces: diseaseMOFA2_*.csv and the model.
# Inputs: results/RA_ILD_Workspace.RData.
PROJECT_ROOT <- Sys.getenv("PROJECT_ROOT", getwd()); setwd(PROJECT_ROOT)
MOFA_PYTHON <- Sys.getenv("MOFA_PYTHON", "/usr/bin/python3"); Sys.setenv(RETICULATE_PYTHON = MOFA_PYTHON)
suppressPackageStartupMessages({
  library(reticulate); use_python(MOFA_PYTHON, required = TRUE)
  library(MOFA2); library(pROC); library(DESeq2)
})
source("scripts/00_cohort_config.R")
cfg <- set_cohort()
tag <- ""
out <- function(f) file.path(cfg$tables, paste0("diseaseMOFA2", tag, "_", f))

load("results/RA_ILD_Workspace.RData")
DROP <- c("BALF_Neutrophil_CD66bpos_CD16pos", "PB_Neutrophil_CD66bpos_CD16pos"); stopifnot(all(DROP %in% colnames(facs_mat))); facs_mat <- facs_mat[, !colnames(facs_mat) %in% DROP]

S <- setdiff(common_samples, cfg$EXCL)
check_n(S, cfg, "MOFA samples")

# ---------------------------------------------------------------- layers ----
top_var <- function(X, k) { v <- apply(X, 2, var, na.rm = TRUE)
  X[, names(sort(v, decreasing = TRUE))[1:min(k, sum(!is.na(v) & v > 0))], drop = FALSE] }

cm <- as.matrix(rnaseq_annotated[, S]); cm <- cm[rowSums(cm) > 0, ]
smd <- master_data[match(S, master_data$Sample_ID), ]; rownames(smd) <- S
dds <- DESeqDataSetFromMatrix(countData = cm, colData = smd, design = ~Sample_Group)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- suppressMessages(DESeq(dds, quiet = TRUE))
vs <- assay(vst(dds, blind = FALSE))
X_expr <- t(vs)[S, ]
X_expr <- top_var(X_expr, 20)

cyto <- function(prefix, exclude = NULL) {
  cc <- grep(prefix, colnames(master_data), value = TRUE)
  if (!is.null(exclude)) cc <- cc[!grepl(exclude, cc)]
  X <- master_data[match(S, master_data$Sample_ID), cc]
  for (c in colnames(X)) X[[c]] <- as.numeric(X[[c]])
  X <- log2(as.matrix(X) + 1); rownames(X) <- S; X
}
X_balf <- cyto("^BALF_", "Treg|Th1|Th2|Th17|Macro|Neutro|Lympho|Eosin|Baso|Plasma_Cell|CD[0-9]|Percent|Activated|Ratio|CD14|CD86|CD66")
X_serum <- cyto("^Serum_")
X_balf <- top_var(X_balf, 15); X_serum <- top_var(X_serum, 15)

X_facs <- facs_mat[S, ]

parts <- strsplit(microbiome_raw$Taxonomy, "; ")
key <- sapply(parts, function(p) if (length(p) >= 2) p[length(p) - 1] else p[1])
A <- as.matrix(microbiome_raw[, S]); A[is.na(A)] <- 0
agg <- rowsum(A, key); rel <- sweep(agg, 2, colSums(agg), "/") * 100
mv <- apply(rel, 1, var)
taxa <- names(sort(mv, decreasing = TRUE))[1:20]
X_micro <- cbind(t(log2(rel[taxa, S, drop = FALSE] + 0.001)),
                 Shannon = diversity_df$Shannon[match(S, diversity_df$Sample_ID)],
                 Simpson = diversity_df$Simpson[match(S, diversity_df$Sample_ID)])

clean <- function(X) { X <- as.matrix(X); X[is.na(X) | is.infinite(X) | is.nan(X)] <- 0
  X[, apply(X, 2, var) > 1e-10, drop = FALSE] }
layers <- list(Expression = X_expr, BALF_Cytokine = X_balf, Serum_Cytokine = X_serum,
               FACS = X_facs, Microbiome = X_micro)
layers <- lapply(layers, clean)

facs <- layers$FACS
d6 <- list(Expression = t(layers$Expression), BALF_Cytokine = t(layers$BALF_Cytokine),
           Serum_Cytokine = t(layers$Serum_Cytokine),
           BALF_FACS = t(facs[, grep("^BALF_", colnames(facs))]),
           PB_FACS = t(facs[, grep("^PB_", colnames(facs))]),
           Microbiome = t(layers$Microbiome))
feat_tab <- do.call(rbind, lapply(names(d6), function(v) data.frame(View = v, Feature = rownames(d6[[v]]))))
write.csv(feat_tab, out("input_features.csv"), row.names = FALSE)

# ------------------------------------------------------------------ MOFA ----
td <- list(V1 = matrix(rnorm(50), 5, 10, dimnames = list(paste0("a", 1:5), paste0("s", 1:10))),
           V2 = matrix(rnorm(30), 3, 10, dimnames = list(paste0("b", 1:3), paste0("s", 1:10))))
mt <- create_mofa(td)
mo_t <- get_default_model_options(mt); mo_t$num_factors <- 2
to_t <- get_default_training_options(mt); to_t$maxiter <- 10; to_t$seed <- 42; to_t$verbose <- FALSE
mt <- suppressWarnings(prepare_mofa(mt, model_options = mo_t, training_options = to_t))
mt <- run_mofa(mt, outfile = tempfile(fileext = ".hdf5"), use_basilisk = FALSE)

mofa <- create_mofa(d6)
mo <- get_default_model_options(mofa); mo$num_factors <- 5
to <- get_default_training_options(mofa)
to$maxiter <- 500; to$seed <- 42; to$verbose <- FALSE; to$convergence_mode <- "medium"
mofa <- suppressWarnings(prepare_mofa(mofa, model_options = mo, training_options = to))
hdf <- file.path(cfg$tables, paste0("diseaseMOFA2", tag, ".hdf5"))
if (file.exists(hdf)) file.remove(hdf)
m6 <- run_mofa(mofa, outfile = hdf, use_basilisk = FALSE)

r2 <- get_variance_explained(m6)
R2f <- r2$r2_per_factor[[1]]; R2t <- r2$r2_total[[1]]
Fm <- get_factors(m6)[[1]]; W <- get_weights(m6)
check_n(rownames(Fm), cfg, "MOFA factors")

# disease factor: the factor with the largest share of serum cytokine variance (Factor 1 with seed 42), sign as returned
kD <- which.max(R2f[, "Serum_Cytokine"]); sgn <- 1
fD <- sgn * Fm[, kD]

# ------------------------------------------------------------ statistics ----
y <- ifelse(is_sarc(rownames(Fm)), 0, 1)
wt <- wilcox.test(fD[y == 1], fD[y == 0], exact = TRUE)
pf <- sapply(1:ncol(Fm), function(k) wilcox.test(Fm[y == 1, k], Fm[y == 0, k], exact = TRUE)$p.value)
rocF <- roc(y, fD, quiet = TRUE, direction = "auto")
set.seed(42); ciF <- ci.auc(rocF, method = "bootstrap", boot.n = 2000)
serum <- master_data[match(rownames(Fm), master_data$Sample_ID), ]
sa <- sapply(grep("^Serum_", colnames(serum), value = TRUE), function(cc) {
  v <- as.numeric(serum[[cc]]); if (sum(!is.na(v)) < 0.85 * length(v)) return(NA)
  tryCatch(as.numeric(auc(roc(y, v, quiet = TRUE, direction = "auto"))), error = function(e) NA) })
sa <- sort(sa, decreasing = TRUE)

met <- data.frame(
  Metric = c("n RA", "n sarcoidosis", "Disease factor (retrained index)", "Disease factor sign applied", "Disease factor R2 Serum_Cytokine (%)",
             "Disease factor R2 BALF_Cytokine (%)", "Disease factor R2 Expression (%)",
             "Disease factor R2 BALF_FACS (%)", "Disease factor R2 PB_FACS (%)", "Disease factor R2 Microbiome (%)",
             "Disease factor Wilcoxon p (exact)", "Disease factor AUC", "Disease factor AUC 95% CI low",
             "Disease factor AUC 95% CI high", "Best single serum cytokine", "Best single serum cytokine AUC"),
  Value = c(sum(y == 1), sum(y == 0), kD, sgn,
            sprintf("%.1f", R2f[kD, c("Serum_Cytokine", "BALF_Cytokine", "Expression", "BALF_FACS", "PB_FACS", "Microbiome")]),
            sprintf("%.3g", wt$p.value), sprintf("%.3f", auc(rocF)), sprintf("%.3f", ciF[1]), sprintf("%.3f", ciF[3]),
            names(sa)[1], sprintf("%.3f", sa[1])))
write.csv(met, out("metrics.csv"), row.names = FALSE)
print(met, row.names = FALSE)

r2df <- data.frame(Factor = rownames(R2f), round(R2f, 2), check.names = FALSE)
r2df <- rbind(r2df, data.frame(Factor = "Total", t(round(R2t, 2)), check.names = FALSE))
write.csv(r2df, out("R2_per_factor_view.csv"), row.names = FALSE)
print(r2df, row.names = FALSE)
write.csv(data.frame(Factor = 1:ncol(Fm), p_disease_exact = pf), out("factor_disease_p.csv"), row.names = FALSE)

write.csv(data.frame(Sample_ID = rownames(Fm), Group = ifelse(y == 1, "RA", "Sarcoidosis"),
                     Fm, Disease_factor_aligned = fD), out("factor_scores.csv"), row.names = FALSE)
aw <- do.call(rbind, lapply(names(W), function(v) { w <- W[[v]]
  do.call(rbind, lapply(1:ncol(w), function(f) { ti <- order(abs(w[, f]), decreasing = TRUE)[1:min(10, nrow(w))]
    data.frame(View = v, Factor = f, Feature = rownames(w)[ti], Weight = w[ti, f]) })) }))
write.csv(aw, out("weights_top10_SD8format.csv"), row.names = FALSE)
allw <- do.call(rbind, lapply(names(W), function(v) data.frame(View = v, Feature = rownames(W[[v]]), W[[v]], check.names = FALSE)))
write.csv(allw, out("weights_all.csv"), row.names = FALSE)
roc_c <- data.frame(Model = "DiseaseFactor", Spec = rocF$specificities, Sens = rocF$sensitivities)
write.csv(roc_c, out("ROC_coords.csv"), row.names = FALSE)
write.csv(data.frame(Serum = names(sa), AUC = sa), out("single_serum_AUC.csv"), row.names = FALSE)
save(m6, Fm, W, r2, layers, d6, kD, sgn, met, file = out("model.RData"))
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))

