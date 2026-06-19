# ============================================================================
# RA-ILD Multi-omics Integrated Analysis Pipeline
# Primary Files Only - For Nature Communications Submission
# ============================================================================

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║    RA-ILD Multi-omics Analysis - Complete Pipeline        ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Environment Setup
suppressPackageStartupMessages({
  pkgs <- c("tidyverse", "readxl", "DESeq2", "edgeR", "WGCNA", "clusterProfiler", 
            "enrichplot", "org.Hs.eg.db", "ReactomePA", "GSVA", "vegan", "mixOmics",
            "pheatmap", "randomForest", "caret", "pROC", "ConsensusClusterPlus")
  for(p in pkgs) {
    if(!requireNamespace(p, quietly=TRUE)) {
      if(p %in% c("DESeq2","edgeR","WGCNA","clusterProfiler","enrichplot","org.Hs.eg.db",
                  "ReactomePA","GSVA","ConsensusClusterPlus")) {
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

rnaseq_raw <- read.csv(file.path(data_dir, "Hosogaya RNAseq/RawCount_heatmap_table.csv"), 
                       row.names=1, check.names=FALSE)
if("Cluster_ID" %in% colnames(rnaseq_raw)) rnaseq_raw <- rnaseq_raw[, !colnames(rnaseq_raw) %in% "Cluster_ID"]

microbiome_raw <- read.csv(file.path(data_dir, "OTU abundance table-2 (merged abundance table).csv"),
                           check.names=FALSE)
sample_cols <- grep("^(KY|RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)
new_names <- gsub("^(KYC?\\d{3,4}).*", "\\1", sample_cols)
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
clinical_data <- read_excel(file.path(data_dir, "RA_BALF_microbiome_clinical_20240430.xlsx"), sheet="Analysis")
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
balf_cytokine <- read_excel(file.path(data_dir, "Hosogaya_multiplex_20250422.xlsx"), sheet="BALF_Analysis")
serum_cytokine <- read_excel(file.path(data_dir, "Hosogaya_multiplex_20250422.xlsx"), sheet="Serum_Analysis")
balf_facs <- read_excel(file.path(data_dir, "FCM_integrated_data_transformed.xlsx"), sheet="BALF_analysis")
pbmc_facs <- read_excel(file.path(data_dir, "FCM_integrated_data_transformed.xlsx"), sheet="Blood_analysis")
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
balf_facs <- balf_facs[!balf_facs$Sample_ID %in% exclude_samples,]
pbmc_facs <- pbmc_facs[!pbmc_facs$Sample_ID %in% exclude_samples,]
if(deconv_loaded) {bp_theta <- bp_theta[!rownames(bp_theta) %in% exclude_samples,]
bp_cell <- bp_cell[!rownames(bp_cell) %in% exclude_samples,]}

all_samples <- unique(c(colnames(rnaseq_raw), grep("^(KY|RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)))
master_data <- data.frame(Sample_ID=all_samples, 
                          Sample_Group=ifelse(grepl("^(KYC|Sarcoidosis)",all_samples),"Control","RA")) %>%
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
  left_join(balf_facs, by="Sample_ID") %>% left_join(pbmc_facs, by="Sample_ID")

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

balf_facs_cols <- grep("^BALF_.*(Treg|Th|Macro|Neutro)", colnames(master_data), value=TRUE)
pbmc_facs_cols <- grep("^PB_", colnames(master_data), value=TRUE)
all_facs_cols <- c(balf_facs_cols, pbmc_facs_cols)
facs_data <- master_data[match(common_samples, master_data$Sample_ID), all_facs_cols]
char_cols <- names(facs_data)[sapply(facs_data, is.character)]
if(length(char_cols)>0) for(col in char_cols) facs_data[[col]] <- as.numeric(facs_data[[col]])
numeric_cols <- sapply(facs_data, is.numeric)
if(sum(!numeric_cols)>0) facs_data <- facs_data[, numeric_cols, drop=FALSE]
facs_mat <- log2(as.matrix(facs_data)+0.01); rownames(facs_mat) <- common_samples

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

# DIABLO
cat("\n=== DIABLO ===\n")
diablo_results <- NULL
tryCatch({
  Y <- as.factor(meta$Sample_Group); names(Y) <- rownames(meta)
  expr_clean <- expr_mat; expr_clean[is.na(expr_clean)|is.infinite(expr_clean)] <- 0
  cyto_clean <- cyto_mat; cyto_clean[is.na(cyto_clean)|is.infinite(cyto_clean)] <- 0
  cyto_clean <- cyto_clean[, colSums(!is.na(cyto_clean))>nrow(cyto_clean)/2]
  facs_clean <- facs_mat; facs_clean[is.na(facs_clean)|is.infinite(facs_clean)] <- 0
  facs_clean <- facs_clean[, colSums(!is.na(facs_clean))>nrow(facs_clean)/2]
  data_list <- list(Expression=expr_clean, Cytokines=cyto_clean, FACS=facs_clean)
  design <- matrix(0.1, ncol=3, nrow=3, dimnames=list(names(data_list),names(data_list)))
  diag(design) <- 0
  diablo_res <- block.splsda(X=data_list, Y=Y, ncomp=2, design=design)
  diablo_results <- list(model=diablo_res, features=selectVar(diablo_res,comp=1), design=design)
  safe_pdf("results/figures/DIABLO_samples.pdf", 
           {plotIndiv(diablo_res, comp=c(1,2), group=Y, legend=TRUE, title="DIABLO")}, width=8, height=7)
  cat("  ✓ Complete\n")
}, error=function(e) cat("  ✗ Skipped\n"))

# WGCNA
cat("\n=== WGCNA ===\n")
wgcna_results <- MEs <- NULL
tryCatch({
  wgcna_expr <- t(expr_matrix)
  gsg <- goodSamplesGenes(wgcna_expr, verbose=3)
  if(!gsg$allOK) wgcna_expr <- wgcna_expr[gsg$goodSamples, gsg$goodGenes]
  sft <- pickSoftThreshold(wgcna_expr, powerVector=c(seq(1,10,1),seq(12,20,2)), verbose=0)
  softPower <- ifelse(is.na(sft$powerEstimate), 6, sft$powerEstimate)
  net <- blockwiseModules(wgcna_expr, power=softPower, TOMType="unsigned", minModuleSize=30,
                          reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE, 
                          pamRespectsDendro=FALSE, verbose=0)
  moduleColors <- labels2colors(net$colors); MEs <- net$MEs
  traits <- data.frame(Disease=as.numeric(meta$Sample_Group=="RA"),
                       Infection=ifelse(is.na(meta$Infection_Group),0,
                                        as.numeric(meta$Infection_Group=="Infection_Positive")))
  moduleTraitCor <- cor(MEs, traits, use="p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(wgcna_expr))
  wgcna_results <- list(net=net, moduleColors=moduleColors, MEs=MEs,
                        moduleTraitCor=moduleTraitCor, moduleTraitPvalue=moduleTraitPvalue)
  safe_write(data.frame(Gene=colnames(wgcna_expr), Module=moduleColors),
             "results/tables/WGCNA_modules.csv", row.names=FALSE)
  cat(sprintf("  Modules: %d\n", length(unique(moduleColors))-1))
}, error=function(e) cat("  ✗ Skipped\n"))

# Consensus Clustering
cat("\n=== Clustering ===\n")
clustering_results <- NULL
tryCatch({
  ra_samples_clust <- common_samples[meta$Sample_Group=="RA"]
  if(length(ra_samples_clust)>=15) {
    clust_expr <- t(expr_matrix[,ra_samples_clust])
    expr_vars <- apply(clust_expr, 2, var, na.rm=TRUE)
    top_expr <- names(sort(expr_vars, decreasing=TRUE))[1:min(100,sum(!is.na(expr_vars)))]
    clust_data <- cbind(clust_expr[,top_expr], cyto_mat[ra_samples_clust,], facs_mat[ra_samples_clust,])
    clust_data[is.na(clust_data)|is.infinite(clust_data)] <- 0
    feature_vars <- apply(clust_data, 2, var, na.rm=TRUE)
    clust_data <- t(scale(t(clust_data[,!is.na(feature_vars) & feature_vars>0])))
    cc_results <- ConsensusClusterPlus(t(clust_data), maxK=4, reps=1000, pItem=0.8, pFeature=0.8,
                                       clusterAlg="km", distance="euclidean", seed=42, plot="pdf",
                                       title="results/consensus_clustering/consensus")
    cluster_assignments <- data.frame(Sample_ID=ra_samples_clust, Cluster_k2=cc_results[[2]]$consensusClass,
                                      Cluster_k3=cc_results[[3]]$consensusClass, Cluster_k4=cc_results[[4]]$consensusClass)
    safe_write(cluster_assignments, "results/tables/Consensus_clusters.csv", row.names=FALSE)
    clustering_results <- list(consensus_results=cc_results, cluster_assignments=cluster_assignments)
    cat(sprintf("  ✓ %d samples\n", length(ra_samples_clust)))
  }
}, error=function(e) cat("  ✗ Skipped\n"))

# LOOCV
cat("\n=== LOOCV ===\n")
loocv_results <- NULL
tryCatch({
  loocv_samples <- intersect(intersect(rownames(expr_mat),rownames(cyto_mat)),rownames(facs_mat))
  loocv_data <- data.frame(expr_mat[loocv_samples,1:min(50,ncol(expr_mat))], cyto_mat[loocv_samples,],
                           facs_mat[loocv_samples,], Group=meta$Sample_Group[match(loocv_samples,meta$Sample_ID)])
  loocv_data <- loocv_data[!is.na(loocv_data$Group), ]                    # base subset to avoid stats::filter collision
  loocv_data$Group <- factor(loocv_data$Group)                             # ensure factor levels are not NULL
  for(col in setdiff(colnames(loocv_data),"Group"))
    loocv_data[[col]][is.na(loocv_data[[col]])|is.infinite(loocv_data[[col]])] <- median(loocv_data[[col]],na.rm=TRUE)
  feature_vars <- apply(loocv_data[, setdiff(colnames(loocv_data),"Group"), drop=FALSE], 2, var, na.rm=TRUE)  # base indexing to avoid dplyr::select collision
  loocv_data <- loocv_data[,c(names(feature_vars[!is.na(feature_vars) & feature_vars>0]),"Group")]
  if(nrow(loocv_data)>=15 && ncol(loocv_data)>=10) {
    n <- nrow(loocv_data); predictions <- rep(NA,n)
    pred_probs <- matrix(NA, nrow=n, ncol=2, dimnames=list(NULL,levels(loocv_data$Group)))
    for(i in 1:n) {
      features <- setdiff(colnames(loocv_data),"Group")
      rf_temp <- randomForest(x=loocv_data[-i,features], y=loocv_data$Group[-i], ntree=500)
      predictions[i] <- as.character(predict(rf_temp, loocv_data[i,features]))
      pred_probs[i,] <- stats::predict(rf_temp, loocv_data[i,features], type="prob")
    }
    conf_mat <- table(Predicted=predictions, Actual=loocv_data$Group)
    roc_obj <- roc(loocv_data$Group, pred_probs[,"RA"])
    loocv_results <- list(predictions=predictions, probabilities=pred_probs, confusion_matrix=conf_mat,
                          accuracy=sum(diag(conf_mat))/sum(conf_mat), AUC=as.numeric(auc(roc_obj)), roc=roc_obj)
    safe_pdf("results/figures/LOOCV_ROC.pdf", 
             {plot(roc_obj, main=sprintf("LOOCV\nAUC=%.3f",loocv_results$AUC), col="#E41A1C", lwd=2)
               abline(0,1,lty=2,col="gray")}, width=6, height=6)
    cat(sprintf("  AUC: %.3f\n", loocv_results$AUC))
  }
}, error=function(e) cat(sprintf("  ✗ LOOCV error: %s\n", e$message)))

# Pathway Analysis
cat("\n=== Pathways ===\n")
deg_genes <- res_df %>% dplyr::filter(Regulation!="NS", !is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique()
background <- res_df %>% dplyr::filter(!is.na(ENTREZID)) %>% pull(ENTREZID) %>% unique()
pathway_results <- list()
if(length(deg_genes)>=10) {
  tryCatch({
    ego <- enrichGO(gene=deg_genes, universe=background, OrgDb=org.Hs.eg.db, ont="BP", 
                    pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
    if(!is.null(ego) && nrow(ego@result)>0) {
      pathway_results$GO_BP <- ego; safe_write(ego@result, "results/tables/GO_BP.csv", row.names=FALSE)}
  }, error=function(e) NULL)
  tryCatch({
    kegg <- enrichKEGG(gene=deg_genes, organism="hsa", universe=background, pvalueCutoff=0.05)
    if(!is.null(kegg) && nrow(kegg@result)>0) {
      pathway_results$KEGG <- kegg; safe_write(kegg@result, "results/tables/KEGG.csv", row.names=FALSE)}
  }, error=function(e) NULL)
  tryCatch({
    reactome <- enrichPathway(gene=deg_genes, organism="human", universe=background, pvalueCutoff=0.05, readable=TRUE)
    if(!is.null(reactome) && nrow(reactome@result)>0) {
      pathway_results$Reactome <- reactome; safe_write(reactome@result, "results/tables/Reactome.csv", row.names=FALSE)}
  }, error=function(e) NULL)
}
cat("✓ Pathways complete\n")

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
  if(ncol(facs_mat)>0) {
    facs_samples <- intersect(ra_samples, rownames(facs_mat))
    if(length(facs_samples)>=10) {
      for(gs in rownames(gsva_scores)) for(facs_col in colnames(facs_mat)) {
        valid_idx <- !is.na(gsva_scores[gs,facs_samples]) & !is.na(facs_mat[facs_samples,facs_col]) & 
          is.finite(gsva_scores[gs,facs_samples]) & is.finite(facs_mat[facs_samples,facs_col])
        if(sum(valid_idx)>=10) {
          sp_test <- cor.test(gsva_scores[gs,facs_samples][valid_idx], facs_mat[facs_samples,facs_col][valid_idx], 
                              method="spearman")
          external_cor <- rbind(external_cor, data.frame(GeneSet=gs, External_Marker=facs_col, Type="FCM",
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

# Module Enrichment
cat("\n=== Module Enrichment ===\n")
module_enrichment_results <- NULL
if(!is.null(wgcna_results)) {
  tryCatch({
    moduleColors <- wgcna_results$moduleColors
    unique_modules <- setdiff(unique(moduleColors), "grey")
    module_enrichment_results <- list()
    for(module in unique_modules) {
      module_genes <- names(moduleColors)[moduleColors==module]
      module_entrez <- gene_annotations %>% dplyr::filter(Gene_Symbol %in% module_genes, !is.na(ENTREZID)) %>% 
        pull(ENTREZID) %>% unique()
      if(length(module_entrez)>=10) {
        ego <- tryCatch(enrichGO(gene=module_entrez, universe=background, OrgDb=org.Hs.eg.db, ont="BP", 
                                 pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE), error=function(e) NULL)
        if(!is.null(ego) && nrow(ego@result)>0) {
          module_enrichment_results[[module]]$GO_BP <- ego
          safe_write(ego@result, sprintf("results/tables/Module_%s_GO_BP.csv",module), row.names=FALSE)
        }
        kegg <- tryCatch(enrichKEGG(gene=module_entrez, organism="hsa", universe=background, pvalueCutoff=0.05),
                         error=function(e) NULL)
        if(!is.null(kegg) && nrow(kegg@result)>0) {
          module_enrichment_results[[module]]$KEGG <- kegg
          safe_write(kegg@result, sprintf("results/tables/Module_%s_KEGG.csv",module), row.names=FALSE)
        }
      }
    }
    cat(sprintf("  %d modules\n", length(module_enrichment_results)))
  }, error=function(e) cat("  ✗ Skipped\n"))
}

# Microbiome
cat("\n=== Microbiome ===\n")
microbiome_samples <- grep("^(KY|RA[0-9]|Sarcoidosis)", colnames(microbiome_raw), value=TRUE)
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

# Cytokine/FACS
cat("\n=== Cytokine/FACS ===\n")
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
balf_facs_stats <- compare_groups(master_data, balf_facs_cols, "Sample_Group", "Control", "RA")
safe_write(balf_facs_stats, "results/tables/FACS_BALF.csv", row.names=FALSE)
cytokine_facs_results <- list(BALF_cytokines=balf_cyto_stats, Serum_cytokines=serum_cyto_stats, BALF_FACS=balf_facs_stats)
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
            sp_test <- cor.test(gene_vals, th17_matched, method="spearman")
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

# Mito-Immune Axis
cat("\n=== Mito-Immune ===\n")
mito_immune_results <- NULL
tryCatch({
  if(!is.null(gsva_results) && !is.null(gsva_results$scores)) {
    gsva_scores <- gsva_results$scores
    mito_avail <- intersect(c("OXPHOS","Mito_Biogenesis","TCA_Cycle"), rownames(gsva_scores))
    immune_avail <- intersect(c("Innate_Myeloid","IFN_Response","TNFa_NFkB","Adaptive_Tcell","Cytotoxicity"), 
                              rownames(gsva_scores))
    if(length(mito_avail)>0 && length(immune_avail)>0) {
      ra_samples_mito <- colnames(gsva_scores)[colnames(gsva_scores) %in% common_samples[meta$Sample_Group=="RA"]]
      mito_score <- colMeans(gsva_scores[mito_avail, ra_samples_mito, drop=FALSE])
      immune_score <- colMeans(gsva_scores[immune_avail, ra_samples_mito, drop=FALSE])
      cor_test <- cor.test(mito_score, immune_score, method="spearman")
      mito_immune_df <- data.frame(Sample_ID=ra_samples_mito, Mito_Score=mito_score, Immune_Score=immune_score)
      safe_write(mito_immune_df, "results/tables/Mito_Immune_scores.csv", row.names=FALSE)
      safe_pdf("results/figures/Mito_Immune_scatter.pdf", {
        plot(mito_score, immune_score, xlab="Mitochondrial Score", ylab="Immune Score",
             main=sprintf("Mito-Immune Axis\nρ=%.3f, p=%.3e", cor_test$estimate, cor_test$p.value),
             pch=19, col="#E41A1C"); abline(lm(immune_score~mito_score), col="blue", lwd=2, lty=2)
      }, width=6, height=6)
      mito_immune_results <- list(mito_pathways=mito_avail, immune_pathways=immune_avail, 
                                  scores=mito_immune_df, correlation=cor_test)
      cat(sprintf("  ρ=%.3f, p=%.3e\n", cor_test$estimate, cor_test$p.value))
    }
  }
}, error=function(e) cat("✗ Skipped\n"))

# Machine Learning
cat("\n=== ML ===\n")
ml_results <- NULL
tryCatch({
  ml_samples <- intersect(rownames(expr_mat), rownames(cyto_mat))
  ml_data <- data.frame(expr_mat[ml_samples,1:min(50,ncol(expr_mat))], cyto_mat[ml_samples,],
                        Sample_Group=meta$Sample_Group[match(ml_samples,meta$Sample_ID)]) %>% dplyr::filter(!is.na(Sample_Group))
  feature_vars <- apply(ml_data %>% dplyr::select(-Sample_Group), 2, var, na.rm=TRUE)
  ml_data <- ml_data[,c(names(feature_vars[!is.na(feature_vars) & feature_vars>0]),"Sample_Group")]
  if(nrow(ml_data)>=15 && ncol(ml_data)>=10) {
    set.seed(42); train_idx <- createDataPartition(ml_data$Sample_Group, p=0.7, list=FALSE)
    train_data <- ml_data[train_idx,]; test_data <- ml_data[-train_idx,]
    features <- setdiff(colnames(ml_data),"Sample_Group")
    for(col in features) {
      train_data[[col]][is.na(train_data[[col]])] <- median(train_data[[col]], na.rm=TRUE)
      test_data[[col]][is.na(test_data[[col]])] <- median(train_data[[col]], na.rm=TRUE)
    }
    rf_model <- randomForest(x=train_data[,features], y=train_data$Sample_Group, ntree=500, importance=TRUE)
    test_pred_prob <- stats::predict(rf_model, newdata=test_data[,features], type="prob")
    test_pred_class <- predict(rf_model, newdata=test_data[,features])
    roc_obj <- roc(test_data$Sample_Group, test_pred_prob[,"RA"])
    conf_mat <- confusionMatrix(test_pred_class, test_data$Sample_Group)
    ml_results <- list(model=rf_model, AUC=as.numeric(auc(roc_obj)), confusion_matrix=conf_mat, roc=roc_obj)
    cat(sprintf("  AUC=%.3f\n", ml_results$AUC))
  }
}, error=function(e) cat("✗ Skipped\n"))

# Infection Prediction
cat("\n=== Infection Prediction ===\n")
infection_prediction <- NULL
tryCatch({
  if("Infection_Group" %in% colnames(meta)) {
    ra_inf_samples <- common_samples[meta$Sample_Group=="RA" & !is.na(meta$Infection_Group)]
    if(length(ra_inf_samples)>=15) {
      inf_data <- data.frame(cyto_mat[ra_inf_samples,], facs_mat[ra_inf_samples,],
                             Infection=meta$Infection_Group[match(ra_inf_samples,meta$Sample_ID)]) %>% dplyr::filter(!is.na(Infection))
      for(col in setdiff(colnames(inf_data),"Infection"))
        inf_data[[col]][is.na(inf_data[[col]])|is.infinite(inf_data[[col]])] <- median(inf_data[[col]],na.rm=TRUE)
      feature_vars <- apply(inf_data %>% dplyr::select(-Infection), 2, var, na.rm=TRUE)
      inf_data <- inf_data[,c(names(feature_vars[!is.na(feature_vars) & feature_vars>0]),"Infection")]
      if(nrow(inf_data)>=15 && ncol(inf_data)>=5) {
        features <- setdiff(colnames(inf_data),"Infection")
        rf_inf <- randomForest(x=inf_data[,features], y=inf_data$Infection, ntree=500, importance=TRUE)
        importance_df <- data.frame(Feature=rownames(importance(rf_inf)),
                                    MeanDecreaseGini=importance(rf_inf)[,"MeanDecreaseGini"]) %>%
          arrange(desc(MeanDecreaseGini))
        safe_write(importance_df, "results/tables/Infection_prediction_importance.csv", row.names=FALSE)
        n_inf <- nrow(inf_data); inf_predictions <- rep(NA,n_inf)
        inf_probs <- matrix(NA, nrow=n_inf, ncol=2, dimnames=list(NULL,levels(inf_data$Infection)))
        for(i in 1:n_inf) {
          rf_temp_inf <- randomForest(x=inf_data[-i,features], y=inf_data$Infection[-i], ntree=500)
          inf_predictions[i] <- as.character(predict(rf_temp_inf, inf_data[i,features]))
          inf_probs[i,] <- stats::predict(rf_temp_inf, inf_data[i,features], type="prob")
        }
        conf_mat_inf <- table(Predicted=inf_predictions, Actual=inf_data$Infection)
        roc_inf <- roc(inf_data$Infection, inf_probs[,"Infection_Positive"])
        infection_prediction <- list(model=rf_inf, importance=importance_df, 
                                     accuracy=sum(diag(conf_mat_inf))/sum(conf_mat_inf),
                                     AUC=as.numeric(auc(roc_inf)), confusion_matrix=conf_mat_inf, roc=roc_inf)
        safe_pdf("results/figures/Infection_prediction_ROC.pdf", {
          plot(roc_inf, main=sprintf("Infection Prediction\nAUC=%.3f",infection_prediction$AUC), 
               col="#FF7F00", lwd=2); abline(0,1,lty=2,col="gray")
        }, width=6, height=6)
        cat(sprintf("  AUC=%.3f\n", infection_prediction$AUC))
      }
    }
  }
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
      ifelse(grepl("^(KYC|Sarcoidosis)", deconv_common), "Control", "RA"),
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
    save(master_data, gene_annotations, count_matrix, log_cpm, expr_matrix, expr_mat, cyto_mat, facs_mat, meta,
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

# Defense Model
cat("\n=== Defense Model ===\n")
defense_model <- NULL
tryCatch({
  if(!is.null(gsva_results) && "Infection_Group" %in% colnames(meta)) {
    gsva_scores <- gsva_results$scores; defense_axes <- list()
    if("OXPHOS" %in% rownames(gsva_scores)) defense_axes$OXPHOS <- gsva_scores["OXPHOS",]
    if("Surfactant" %in% rownames(gsva_scores)) defense_axes$Surfactant <- gsva_scores["Surfactant",]
    th17_col <- grep("Th17\\.1|Th17_1", colnames(master_data), value=TRUE, ignore.case=TRUE)
    if(length(th17_col)>0) {
      th17_vals <- as.numeric(master_data[[th17_col[1]]][match(colnames(gsva_scores),master_data$Sample_ID)])
      names(th17_vals) <- colnames(gsva_scores); defense_axes$Th17_1 <- th17_vals
    }
    if(length(defense_axes)>=2) {
      ra_samples_def <- colnames(gsva_scores)[colnames(gsva_scores) %in% common_samples[meta$Sample_Group=="RA"]]
      defense_mat <- matrix(NA, nrow=length(defense_axes), ncol=length(ra_samples_def))
      rownames(defense_mat) <- names(defense_axes); colnames(defense_mat) <- ra_samples_def
      for(axis_name in names(defense_axes)) defense_mat[axis_name,] <- scale(defense_axes[[axis_name]][ra_samples_def])[,1]
      composite_score <- colMeans(defense_mat, na.rm=TRUE)
      ra_meta_def <- meta[ra_samples_def,]
      inf_neg_def <- ra_samples_def[!is.na(ra_meta_def$Infection_Group) & 
                                      ra_meta_def$Infection_Group=="Infection_Negative"]
      inf_pos_def <- ra_samples_def[!is.na(ra_meta_def$Infection_Group) & 
                                      ra_meta_def$Infection_Group=="Infection_Positive"]
      if(length(inf_neg_def)>=3 && length(inf_pos_def)>=3) {
        wt_def <- wilcox.test(composite_score[inf_neg_def], composite_score[inf_pos_def], exact=TRUE)
        defense_infection <- factor(ifelse(ra_samples_def %in% inf_pos_def,"Positive","Negative"),
                                    levels=c("Negative","Positive"))
        roc_def <- roc(defense_infection, composite_score)
        defense_scores_df <- data.frame(Sample_ID=ra_samples_def, t(defense_mat), Composite_Score=composite_score)
        safe_write(defense_scores_df, "results/tables/Defense_model_scores.csv", row.names=FALSE)
        defense_model <- list(axes=names(defense_axes), defense_matrix=defense_mat, composite_score=composite_score,
                              p_value=wt_def$p.value, AUC=as.numeric(auc(roc_def)))
        cat(sprintf("  %d axes, AUC=%.3f\n", length(defense_axes), defense_model$AUC))
      }
    }
  }
}, error=function(e) cat("✗ Skipped\n"))

# Network
cat("\n=== Network ===\n")
network_results <- NULL
tryCatch({
  ra_samples_net <- common_samples[meta$Sample_Group=="RA"]
  gene_vars_net <- apply(log_cpm[,ra_samples_net], 1, var, na.rm=TRUE)
  top_genes_net <- names(sort(gene_vars_net, decreasing=TRUE))[1:50]
  cyto_vars <- apply(cyto_mat[ra_samples_net,], 2, var, na.rm=TRUE)
  top_cyto <- names(sort(cyto_vars, decreasing=TRUE))[1:20]
  facs_vars <- apply(facs_mat[ra_samples_net,], 2, var, na.rm=TRUE)
  top_facs <- names(sort(facs_vars, decreasing=TRUE))[1:20]
  network_data <- cbind(t(log_cpm[top_genes_net,ra_samples_net]), cyto_mat[ra_samples_net,top_cyto],
                        facs_mat[ra_samples_net,top_facs])
  network_cor <- cor(network_data, method="spearman", use="pairwise.complete.obs")
  edges <- data.frame()
  for(i in 1:(nrow(network_cor)-1)) for(j in (i+1):ncol(network_cor)) {
    if(abs(network_cor[i,j])>0.5 && !is.na(network_cor[i,j])) {
      node1 <- rownames(network_cor)[i]; node2 <- colnames(network_cor)[j]
      type1 <- ifelse(node1 %in% top_genes_net,"Gene",ifelse(node1 %in% top_cyto,"Cytokine","FCM"))
      type2 <- ifelse(node2 %in% top_genes_net,"Gene",ifelse(node2 %in% top_cyto,"Cytokine","FCM"))
      edges <- rbind(edges, data.frame(Node1=node1, Node2=node2, Type1=type1, Type2=type2,
                                       Correlation=network_cor[i,j], 
                                       Edge_Type=ifelse(type1==type2,"Intra","Inter"), stringsAsFactors=FALSE))
    }
  }
  safe_write(edges, "results/tables/Network_edges.csv", row.names=FALSE)
  network_stats <- edges %>% group_by(Edge_Type) %>% 
    summarise(N_edges=n(), Mean_abs_cor=mean(abs(Correlation)), .groups="drop")
  network_results <- list(correlation_matrix=network_cor, edges=edges, stats=network_stats)
  cat(sprintf("  %d nodes, %d edges\n", nrow(network_cor), nrow(edges)))
}, error=function(e) cat("✗ Skipped\n"))

# Save Results
cat("\n=== Saving ===\n")
results_list <- list(
  master_data=master_data, gene_annotations=gene_annotations, count_matrix=count_matrix, log_cpm=log_cpm,
  expr_matrix=expr_matrix, expr_mat=expr_mat, cyto_mat=cyto_mat, facs_mat=facs_mat, meta=meta,
  common_samples=common_samples, DEG_RA_vs_Control=res_df, DIABLO=diablo_results, WGCNA=wgcna_results,
  Clustering=clustering_results, LOOCV=loocv_results, Pathways=pathway_results, 
  Module_Enrichment=module_enrichment_results, GSVA=gsva_results, Microbiome=microbiome_results,
  Cytokines_FACS=cytokine_facs_results, Th17_Correlations=th17_correlations, 
  Mito_Immune_Axis=mito_immune_results, MachineLearning=ml_results, Infection_Prediction=infection_prediction,
  Deconvolution=deconv_integration_results, Surfactant_Analysis=surfactant_analysis, 
  Defense_Model=defense_model, Network=network_results, analysis_date=Sys.time(), R_version=R.version.string
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
if(!is.null(diablo_results)) cat("DIABLO: ✓\n")
if(!is.null(wgcna_results)) cat(sprintf("WGCNA: %d modules\n", length(unique(wgcna_results$moduleColors))-1))
if(!is.null(ml_results)) cat(sprintf("ML AUC: %.3f\n", ml_results$AUC))
if(deconv_loaded && !is.null(deconv_integration_results)) 
  cat(sprintf("Deconvolution: %d types\n", length(deconv_integration_results$cell_types)))
cat("\n═══════════════════════════════════════════════════════════\n")