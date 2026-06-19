# ============================================================================
# 20260604_regen_Fig5_panels_descriptive.R
# Figure 5 panels with descriptive "Factor 1 / Factor 2" labels, and Fig 5b
# showing Factor 1 only (RA vs Sarcoidosis).
#   5a R2 heatmap  : Factor labels
#   5b boxplot     : Factor 1 only (RA vs Sarcoidosis)
#   5c disease ROC : Factor 1 vs best serum biomarker
# Outputs: Fig5a_R2_noacronym, Fig5b_Factor1_only, Fig5c_diseaseROC_noacronym
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({ library(tidyverse); library(pROC); library(patchwork) })
load("results/RA_ILD_Workspace.RData")
e <- new.env(); load("results/MOFA2_6views_Results.RData", envir = e)
factors <- get("factors", e); r2 <- get("r2", e)
panel_dir <- "results/panels"

view_labels <- c(Expression="Expression", BALF_Cytokine="BAL Cytokine",
                 Serum_Cytokine="Serum Cytokine", BALF_FACS="BAL FCM",
                 PB_FACS="PB FCM", Microbiome="Microbiome")
view_cols <- c(Expression="#E67E22", BALF_Cytokine="#3498DB", Serum_Cytokine="#E74C3C",
               BALF_FACS="#27AE60", PB_FACS="#9B59B6", Microbiome="#95A5A6")

tn <- function(bs = 8) {
  theme_classic(base_size = bs, base_family = "Helvetica") %+replace%
    theme(text = element_text(color = "black"),
          axis.text = element_text(size = rel(0.9), color = "black"),
          axis.title = element_text(size = rel(0.95), face = "bold"),
          plot.title = element_text(size = rel(1.05), face = "bold", hjust = 0,
                                    margin = ggplot2::margin(0,0,3,0)),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 0.33),
          axis.line = element_blank(), plot.margin = unit(c(2,3,2,3),"mm"))
}
sv <- function(p, name, w, h) {
  ggsave(file.path(panel_dir, paste0(name,".pdf")), p, width=w, height=h, units="mm", dpi=300, device=cairo_pdf)
  ggsave(file.path(panel_dir, paste0(name,".png")), p, width=w, height=h, units="mm", dpi=300, bg="white")
  cat("  saved", name, "\n")
}
group_m <- ifelse(grepl("^(KYC|Sarcoidosis)", rownames(factors)), "Sarcoidosis", "RA")

# ---------- 5a: R2 heatmap with Factor 1..5 (no IDD/IIV) ----------
fnames <- c(Factor1="Factor 1", Factor2="Factor 2", Factor3="Factor 3",
            Factor4="Factor 4", Factor5="Factor 5")
r2m <- r2$r2_per_factor[[1]]
r2l <- as.data.frame(r2m) %>% rownames_to_column("Factor") %>%
  pivot_longer(-Factor, names_to="View", values_to="R2")
r2l$Factor <- factor(r2l$Factor, levels = rev(rownames(r2m)))
r2l$View <- factor(r2l$View, levels = colnames(r2m))
r2l$FLabel <- factor(fnames[as.character(r2l$Factor)],
                     levels = rev(fnames[rownames(r2m)]))
p_hm <- ggplot(r2l, aes(View, FLabel, fill = R2)) +
  geom_tile(color="white", linewidth=0.8) +
  scale_fill_gradient(low="white", high="#2471A3", name="R² (%)") +
  geom_text(aes(label=ifelse(R2>0.5, sprintf("%.1f",R2), ""),
                color=ifelse(R2>25,"white","black")), size=2.5, fontface="bold", show.legend=FALSE) +
  scale_color_identity() + scale_x_discrete(labels=view_labels) +
  labs(x="", y="") +
  tn(8) + theme(axis.text.x=element_text(angle=40, hjust=1, size=7,
                  color=view_cols[levels(r2l$View)]),
                legend.position="left", legend.key.height=unit(10,"mm"))
total_r2 <- data.frame(View=factor(names(r2$r2_total[[1]]), levels=colnames(r2m)),
                       Total=as.numeric(r2$r2_total[[1]]))
p_bar <- ggplot(total_r2, aes(Total, View, fill=View)) +
  geom_col(width=0.65) + scale_fill_manual(values=view_cols) +
  scale_y_discrete(labels=view_labels) + labs(x="Total R² (%)", y="") +
  tn(7) + theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank())
sv(p_hm + p_bar + plot_layout(widths=c(4,1.2)), "Fig5a_R2_noacronym", 130, 75)

# ---------- 5b: Factor 1 only, RA vs Sarcoidosis ----------
pval <- suppressWarnings(wilcox.test(factors[group_m=="RA",1], factors[group_m=="Sarcoidosis",1],
                                     exact=TRUE)$p.value)
plab <- if (pval < 0.0001) "p < 0.0001" else sprintf("p = %.4f", pval)
cat(sprintf("  Factor 1 disease wilcox exact p = %.6f -> %s\n", pval, plab))
df5b <- data.frame(Group=factor(group_m, levels=c("RA","Sarcoidosis")), Value=factors[,1])
sv(ggplot(df5b, aes(Group, Value, fill=Group)) +
     geom_boxplot(outlier.shape=NA, width=0.6, linewidth=0.3) +
     geom_jitter(width=0.15, size=0.9, alpha=0.6, color="grey25") +
     annotate("text", x=1.5, y=max(df5b$Value), label=plab, size=2.8, fontface="italic", vjust=-0.2) +
     scale_fill_manual(values=c(RA="#C0392B", Sarcoidosis="#2980B9")) +
     scale_y_continuous(expand=expansion(mult=c(0.05,0.15))) +
     labs(x="", y="Factor 1 value") + ggtitle("Factor 1") +
     tn(8) + theme(legend.position="none"),
   "Fig5b_Factor1_only", 62, 75)

# ---------- 5c: disease ROC with "Factor 1" (no IDD) ----------
y_disease <- ifelse(group_m=="RA",1,0)
roc_f1 <- roc(y_disease, factors[,1], quiet=TRUE, direction="auto")
serum <- master_data[match(rownames(factors), master_data$Sample_ID), ]
best<-NULL; ba<-0
for(cc in grep("^Serum_", colnames(serum), value=TRUE)){
  v<-as.numeric(serum[[cc]])
  if(sum(!is.na(v))>=30){ r<-tryCatch(roc(y_disease,v,quiet=TRUE,direction="auto"),error=function(x)NULL)
    if(!is.null(r)&&as.numeric(auc(r))>ba){ba<-as.numeric(auc(r));best<-cc}}}
roc_best<-roc(y_disease, as.numeric(serum[[best]]), quiet=TRUE, direction="auto")
cf<-gsub("^Serum_","",best); cf<-gsub("^GCSF$","G-CSF",cf); cf<-gsub("^IL(\\d)","IL-\\1",cf)
roc_df<-rbind(
  data.frame(Sens=roc_f1$sensitivities, Spec=1-roc_f1$specificities,
             Model=sprintf("Factor 1 (AUC = %.3f)", auc(roc_f1))),
  data.frame(Sens=roc_best$sensitivities, Spec=1-roc_best$specificities,
             Model=sprintf("%s (AUC = %.3f)", cf, auc(roc_best))))
cat(sprintf("  Factor 1 disease AUC=%.3f, %s AUC=%.3f\n", auc(roc_f1), cf, auc(roc_best)))
sv(ggplot(roc_df, aes(Spec, Sens, color=Model)) +
     geom_step(linewidth=0.6, direction="vh") +
     geom_abline(slope=1, intercept=0, linetype="dotted", color="grey60", linewidth=0.3) +
     annotate("segment", x=0, xend=1, y=1, yend=1, linetype="dashed", color="#27AE60", linewidth=0.4) +
     annotate("text", x=0.98, y=0.05, size=2.3, hjust=1, color="#27AE60",
              label="Nested LOOCV (AUC = 0.962)", fontface="italic") +
     scale_color_manual(values=c("#E67E22","#C0392B")) +
     labs(x="1 – Specificity", y="Sensitivity", color="") + ggtitle("RA vs Sarcoidosis") +
     tn(8) + theme(legend.position=c(0.62,0.22),
                   legend.background=element_rect(fill=alpha("white",0.95), color=NA),
                   legend.text=element_text(size=6.5), legend.key.height=unit(3,"mm")),
   "Fig5c_diseaseROC_noacronym", 80, 75)
