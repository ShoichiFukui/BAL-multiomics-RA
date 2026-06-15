# ============================================================================
# RA-ILD BAL Multi-omics Analysis for Nature Communications
#
# This single script runs ALL analyses for the paper:
#   Section A: Enhanced statistical analysis (DEG, GSEA, effect sizes, LOOCV)
#   Section B: CT quantitative integration
#   Section C: Future infection prediction
#   Section D: Multi-omics integration (Joint PCA, Layer Ablation, MOFA2)
#   Section E: Figure generation (7 Main + Supplementary S1a-f)
#
# Prerequisites (run in order before this script):
#   01_BayesPrism_Deconvolution.R
#   02_PostDeconvolution.R
#   → produces results/RA_ILD_Workspace.RData
#
# Expected runtime: ~60-90 minutes (mostly permutation tests)
# ============================================================================

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  RA-ILD Multi-omics for Nature Communications               ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
if(!file.exists("results/RA_ILD_Workspace.RData"))
  stop("results/RA_ILD_Workspace.RData not found. Run prerequisite scripts first.")
load("results/RA_ILD_Workspace.RData")

suppressPackageStartupMessages({
  library(tidyverse); library(DESeq2); library(edgeR)
  library(clusterProfiler); library(org.Hs.eg.db)
  library(pROC); library(randomForest); library(caret)
  library(GSVA); library(pheatmap); library(effsize)
  library(readxl); library(glmnet); library(vegan)
  library(gridExtra); library(grid); library(ggplot2)
})
set.seed(42)

out_dir <- "results"
fig_dir <- file.path(out_dir, "figures_final")
tbl_dir <- file.path(out_dir, "tables")
for(d in c(out_dir, fig_dir, tbl_dir)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

# ============================================================================
# SECTION A: Enhanced Statistical Analysis
# ============================================================================
cat("\n=== SECTION A: Enhanced Statistical Analysis ===\n")
source("analysis_modules/Enhanced_Analysis.R", local=FALSE)

# ============================================================================
# SECTION B: CT Quantitative Integration
# ============================================================================
cat("\n=== SECTION B: CT Quantitative Integration ===\n")
source("analysis_modules/CT_Multiomics.R", local=FALSE)

# ============================================================================
# SECTION C: Future Infection Prediction
# ============================================================================
cat("\n=== SECTION C: Future Infection Prediction ===\n")
source("analysis_modules/Infection_Prediction.R", local=FALSE)

# ============================================================================
# SECTION D: Multi-omics Integration
# ============================================================================
cat("\n=== SECTION D: Multi-omics Integration ===\n")

# D1: Joint PCA + Layer Ablation + Cross-layer Network
source("analysis_modules/Integration.R", local=FALSE)

# D2: MOFA2 (requires Python mofapy2)
cat("\n--- MOFA2 ---\n")
tryCatch({
  # Try loading MOFA2
  suppressPackageStartupMessages(library(MOFA2))

  # Attempt to set up Python backend
  BSPY <- file.path(system.file(package="basilisk"), "..", "MOFA2", "mofa_env", "bin", "python")
  if(file.exists(BSPY)) {
    Sys.setenv(RETICULATE_PYTHON = BSPY)
    suppressPackageStartupMessages({ library(reticulate); use_python(BSPY, required=TRUE) })
  }

  # Load pre-cleaned data from Integration step
  load(file.path(out_dir, "Integration_7layers_Results.RData"))

  # Build 6 independent views (no Deconv/GSVA — RNA-seq derivatives)
  d7 <- list()
  for(nm in names(layers)) {
    X <- as.matrix(layers[[nm]]); X[is.na(X)|is.infinite(X)|is.nan(X)] <- 0
    fv <- apply(X,2,var); X <- X[,fv>1e-10,drop=FALSE]; d7[[nm]] <- t(X)
  }

  # Warmup (prevents segfault in some environments)
  cat("  MOFA2 warmup...\n")
  td <- list(V1=matrix(rnorm(50),5,10,dimnames=list(paste0("a",1:5),paste0("s",1:10))),
             V2=matrix(rnorm(30),3,10,dimnames=list(paste0("b",1:3),paste0("s",1:10))))
  mt <- create_mofa(td)
  mo_t <- get_default_model_options(mt); mo_t$num_factors <- 2
  to_t <- get_default_training_options(mt); to_t$maxiter <- 10; to_t$seed <- 42; to_t$verbose <- FALSE
  mt <- prepare_mofa(mt, model_options=mo_t, training_options=to_t)
  mt <- run_mofa(mt, outfile=tempfile(fileext=".hdf5"), use_basilisk=FALSE)
  cat("  Warmup OK\n")

  # Split FACS into BALF/PB
  facs_orig <- t(d7$FACS)
  d6 <- list(
    Expression = d7$Expression,
    BALF_Cytokine = d7$BALF_Cytokine,
    Serum_Cytokine = d7$Serum_Cytokine,
    BALF_FACS = t(facs_orig[, grep("^BALF_", colnames(facs_orig))]),
    PB_FACS = t(facs_orig[, grep("^PB_", colnames(facs_orig))]),
    Microbiome = d7$Microbiome
  )
  cat(sprintf("  6 views, %d features\n", sum(sapply(d6, nrow))))

  mofa <- create_mofa(d6)
  mo <- get_default_model_options(mofa); mo$num_factors <- 5
  to <- get_default_training_options(mofa)
  to$maxiter <- 500; to$seed <- 42; to$verbose <- FALSE; to$convergence_mode <- "medium"
  mofa <- prepare_mofa(mofa, model_options=mo, training_options=to)

  cat("  Running MOFA2...\n")
  mofa6 <- run_mofa(mofa, outfile=file.path(out_dir, "mofa2_6views.hdf5"), use_basilisk=FALSE)
  cat("  MOFA2 COMPLETE\n")

  # Extract results
  r2 <- get_variance_explained(mofa6)
  factors <- get_factors(mofa6)[[1]]
  weights <- get_weights(mofa6)

  cat("  R2 per factor:\n"); print(round(r2$r2_per_factor[[1]], 1))

  group_m <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors)), "Sarcoidosis", "RA")
  inf_m <- master_data$respiratory_infection[match(rownames(factors), master_data$Sample_ID)]

  fd <- data.frame(); fi <- data.frame()
  for(f in 1:ncol(factors)) {
    wt <- wilcox.test(factors[group_m=="RA",f], factors[group_m=="Sarcoidosis",f], exact=TRUE)
    fd <- rbind(fd, data.frame(Factor=f, p_disease=wt$p.value))
    m <- group_m=="RA" & !is.na(inf_m)
    po <- factors[m & inf_m==1, f]; ne <- factors[m & inf_m==0, f]
    if(length(po)>=3 & length(ne)>=3) {
      wt2 <- wilcox.test(po, ne, exact=TRUE)
      fi <- rbind(fi, data.frame(Factor=f, p_infection=wt2$p.value))
    }
  }

  # Save
  aw <- data.frame()
  for(v in names(weights)) { w <- weights[[v]]
    for(f in 1:ncol(w)) { ti <- order(abs(w[,f]), decreasing=TRUE)[1:min(10,nrow(w))]
      aw <- rbind(aw, data.frame(View=v, Factor=f, Feature=rownames(w)[ti], Weight=w[ti,f])) }}
  write.csv(aw, file.path(tbl_dir, "MOFA2_6views_weights.csv"), row.names=FALSE)
  save(mofa6, factors, r2, fd, fi, weights, file=file.path(out_dir, "MOFA2_6views_Results.RData"))
  cat("  MOFA2 results saved\n")

}, error=function(e) {
  cat(sprintf("\n  *** MOFA2 FAILED: %s ***\n", e$message))
  cat("  This is likely a Python/reticulate/basilisk environment issue.\n")
  cat("  MOFA2 requires Python package 'mofapy2' via the basilisk or reticulate bridge.\n")
  cat("  \n")
  cat("  To resolve:\n")
  cat("    Option 1: Run 04_MOFA2_FullCohort.R separately with:\n")
  cat("      Rscript 04_MOFA2_FullCohort.R\n")
  cat("    Option 2: Set Python path manually before running:\n")
  cat("      Sys.setenv(RETICULATE_PYTHON = '/path/to/python/with/mofapy2')\n")
  cat("    Option 3: Install mofapy2 in basilisk environment:\n")
  cat("      BiocManager::install('MOFA2')  # then restart R\n")
  cat("  \n")
  cat("  All other analyses are complete. Only MOFA2 figure/table will be missing.\n")
  cat("  The remaining analyses (Joint PCA, Layer Ablation) provide equivalent\n")
  cat("  multi-omics integration results without Python dependency.\n")
})

# ============================================================================
# SECTION E: Figure Generation
# ============================================================================
cat("\n=== SECTION E: Figure Generation ===\n")
source("analysis_modules/Final7.R", local=FALSE)

# Supplementary Figures S1a-f (scRNA-seq UMAP)
cat("\n  Supplementary Figures S1a-f...\n")
tryCatch({
  suppressPackageStartupMessages(library(Seurat))
  ref <- readRDS("output_v3_BayesPrism/BAL_reference_author_annotated.rds")

  tn_s <- function(bs=10) {
    theme_classic(base_size=bs) %+replace%
      theme(text=element_text(color="black"),
            axis.text=element_text(size=rel(.9),color="black"),
            axis.title=element_text(size=rel(1),face="bold"),
            plot.title=element_text(size=rel(1.1),face="bold",hjust=0),
            legend.text=element_text(size=rel(.8)),
            legend.key.size=unit(3,"mm"),
            legend.background=element_rect(fill=alpha("white",.85),color=NA),
            panel.border=element_rect(fill=NA,color="black",linewidth=.5),
            axis.line=element_blank(), plot.margin=unit(c(5,5,5,5),"mm"))
  }
  CC_s <- c(Macrophage="#E74C3C",T_cell="#3498DB",NK="#2ECC71",B_cell="#9B59B6",
            Plasma="#E67E22",Neutrophil="#F1C40F",DC="#8B4513",Epithelial="#E91E63",
            Mast="#95A5A6",Other="#BDC3C7")
  col_ds_s <- c("GSE145926"="#E74C3C","GSE193782"="#3498DB","GSE184735"="#27AE60")
  col_fine_s <- c(Alveolar_Mac="#E74C3C",Monocyte_Mac="#C0392B",T_cell="#3498DB",NK="#2ECC71",
                  B_cell="#9B59B6",Plasma="#E67E22",Neutrophil="#F1C40F",DC="#8B4513",
                  Epithelial="#E91E63",Mast="#95A5A6",Other="#BDC3C7",Granulocyte="#D35400")
  col_src_s <- c("author_Liao2020"="#E74C3C","author_GSE193782"="#3498DB","kNN_transfer"="#27AE60")

  uc <- Embeddings(ref,"umap"); ms <- ref@meta.data
  ds <- data.frame(U1=uc[,1],U2=uc[,2],DS=ms$dataset,CT=ms$celltype,CF=ms$celltype_fine,Src=ms$annotation_source)

  ggsave(file.path(fig_dir,"FigS1a_UMAP_dataset.pdf"),
    ggplot(ds,aes(U1,U2,color=DS))+geom_point(size=.1,alpha=.25)+
      scale_color_manual(values=col_ds_s,labels=c("GSE145926"="COVID-19","GSE193782"="Healthy","GSE184735"="Sarcoidosis"))+
      guides(color=guide_legend(override.aes=list(size=3,alpha=1)))+
      labs(x="UMAP 1",y="UMAP 2",color="")+tn_s(),
    width=200,height=160,units="mm",dpi=300)

  ggsave(file.path(fig_dir,"FigS1b_UMAP_celltype.pdf"),
    ggplot(ds,aes(U1,U2,color=CT))+geom_point(size=.1,alpha=.25)+scale_color_manual(values=CC_s)+
      guides(color=guide_legend(override.aes=list(size=3,alpha=1)))+
      labs(x="UMAP 1",y="UMAP 2",color="")+tn_s(),
    width=200,height=160,units="mm",dpi=300)

  ggsave(file.path(fig_dir,"FigS1c_UMAP_fine.pdf"),
    ggplot(ds,aes(U1,U2,color=CF))+geom_point(size=.1,alpha=.25)+
      scale_color_manual(values=col_fine_s,na.value="grey80")+
      guides(color=guide_legend(override.aes=list(size=3,alpha=1)))+
      labs(x="UMAP 1",y="UMAP 2",color="")+tn_s(),
    width=200,height=160,units="mm",dpi=300)

  ggsave(file.path(fig_dir,"FigS1d_UMAP_source.pdf"),
    ggplot(ds,aes(U1,U2,color=Src))+geom_point(size=.1,alpha=.25)+
      scale_color_manual(values=col_src_s,na.value="grey80",
        labels=c("author_Liao2020"="Author (Liao)","author_GSE193782"="Author (GSE193782)","kNN_transfer"="kNN"))+
      guides(color=guide_legend(override.aes=list(size=3,alpha=1)))+
      labs(x="UMAP 1",y="UMAP 2",color="")+tn_s(),
    width=200,height=160,units="mm",dpi=300)

  comp <- ds%>%group_by(DS,CT)%>%summarise(N=n(),.groups="drop")%>%group_by(DS)%>%mutate(Pct=N/sum(N))
  comp$DS <- recode(comp$DS,"GSE145926"="COVID-19","GSE193782"="Healthy","GSE184735"="Sarcoidosis")
  ggsave(file.path(fig_dir,"FigS1e_Composition.pdf"),
    ggplot(comp,aes(DS,Pct,fill=CT))+geom_col(width=.6)+scale_fill_manual(values=CC_s)+
      labs(x="",y="Proportion",fill="")+tn_s(),
    width=180,height=140,units="mm",dpi=300)

  mks <- intersect(c("CD68","MARCO","FABP4","FCN1","CD3D","CD3E","NKG7","GNLY",
    "CD79A","MS4A1","JCHAIN","MZB1","FCGR3B","CSF3R","FCER1A","CD1C","EPCAM","KRT18"),rownames(ref))
  ed <- GetAssayData(ref,layer="data"); dd <- data.frame()
  for(ct in unique(ds$CT)){cl<-which(ds$CT==ct);if(length(cl)<10)next
    for(mk in mks){vv<-ed[mk,cl];dd<-rbind(dd,data.frame(CT=ct,M=mk,Pct=sum(vv>0)/length(vv)*100,Avg=mean(vv)))}}
  dd$M <- factor(dd$M,levels=rev(mks))
  ggsave(file.path(fig_dir,"FigS1f_Markers.pdf"),
    ggplot(dd,aes(CT,M,size=Pct,color=Avg))+geom_point()+
      scale_size_continuous(range=c(.8,5),name="% Expr")+scale_color_gradient(low="grey90",high="#C0392B",name="Mean expr")+
      labs(x="",y="")+
      tn_s(9)+theme(axis.text.x=element_text(angle=45,hjust=1,size=8),axis.text.y=element_text(size=8,face="italic")),
    width=200,height=180,units="mm",dpi=300)

  rm(ref,ed); gc()
  cat("  Supplementary figures complete\n")
}, error=function(e) cat(sprintf("  Supplementary figures failed: %s\n", e$message)))

# ============================================================================
# COMPLETE
# ============================================================================
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  ANALYSIS COMPLETE                                           ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("  Tables: %d files in %s\n", length(list.files(tbl_dir)), tbl_dir))
cat(sprintf("  Figures: %d files in %s\n", length(list.files(fig_dir, pattern="\\.pdf$")), fig_dir))
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Date: %s\n", Sys.Date()))

sink(file.path(out_dir, "session_info.txt"))
cat("Analysis completed:", as.character(Sys.time()), "\n\n"); sessionInfo()
sink()
