# ============================================================================
# 20260605_rebuild_supp_tables_renumbered.R
# Rebuild the 9-sheet Supplementary Tables workbook in the CORRECTED first-mention
# order (Nature Comms), with the same content/format fixes as the non-renumbered
# build plus: (a) renumbering, (b) ST2 (healthy serum) healthy-vs-X test columns
# removed (study protocol forbids hypothesis testing of the healthy comparison;
# only KW omnibus + RA-vs-Sarc kept).
# New order (display ST1..ST9) follows body first-mention order 1,8,2,3,4,5,6,7,9:
#   ST1 = clinical RA vs sarcoidosis            (old ST1, verbatim)
#   ST2 = serum cytokine RA/sarc/healthy        (old ST8; drop RA-vs-Healthy & Sarc-vs-Healthy P)
#   ST3 = clinical by infection status          (old ST2)
#   ST4 = single biomarker ROC (infection)      (old ST3, 23 rows desc AUC)
#   ST5 = BAL microbiome by infection           (old ST4, all, asc P)
#   ST6 = clinical vs multi-omics               (old ST5, clean 17)
#   ST7 = BAL/serum profiles by ILD             (old ST6, all 108, asc P)
#   ST8 = CT progression correlations           (old ST7, all 432, asc p)
#   ST9 = medication confounding                (old ST9; IDD->Factor 1, IIV->Factor 2)
# Output: manuscript&figures/20260605_Supplementary_Tables.xlsx (overwrite)
# ============================================================================
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages(library(openxlsx))
TBL <- "results/tables"
SRC_XLSX <- "manuscript&figures/20260601_Supplementary_Tables.xlsx"
OUT <- "manuscript&figures/20260605_Supplementary_Tables.xlsx"
MINUS<-"−"; SUB2<-"₂"; RHO<-"ρ"; ENDASH<-"–"
rnd<-function(x,d) ifelse(is.na(x),NA,round(as.numeric(x),d))
sgf<-function(x,d=3) ifelse(is.na(x),NA,signif(as.numeric(x),d))

title_style  <- createStyle(fontSize=11, textDecoration="bold")
header_style <- createStyle(fontSize=10, textDecoration="bold", border="TopBottom",
                            borderStyle="thin", fgFill="#F2F2F2", halign="left",
                            valign="center", wrapText=TRUE)
data_style   <- createStyle(fontSize=10, halign="left", valign="top")
section_style<- createStyle(fontSize=10, textDecoration="bold")
wb <- createWorkbook()
add_sheet <- function(sheet,title,headers,dat){
  addWorksheet(wb,sheet)
  writeData(wb,sheet,title,startRow=1,startCol=1); addStyle(wb,sheet,title_style,rows=1,cols=1)
  colnames(dat)<-headers
  writeData(wb,sheet,dat,startRow=2,startCol=1,headerStyle=header_style)
  if(nrow(dat)>0) addStyle(wb,sheet,data_style,rows=3:(2+nrow(dat)),cols=1:ncol(dat),gridExpand=TRUE,stack=TRUE)
  setColWidths(wb,sheet,cols=1:ncol(dat),widths="auto")
  cat(sprintf("  %s: %d x %d  | %s\n",sheet,nrow(dat),ncol(dat),substr(title,1,45)))
}

## ST1 = clinical RA vs sarc (verbatim from existing workbook)
# readWorkbook can return numeric character refs as literal "&#NNNN;" strings; decode them.
dehtml<-function(x){
  if(is.na(x)) return(x)
  repeat{ g<-regmatches(x,regexpr("&#[0-9]+;",x)); if(length(g)==0) break
    x<-sub("&#[0-9]+;", intToUtf8(as.integer(sub("&#([0-9]+);","\\1",g))), x) }
  x<-gsub("&gt;",">",x,fixed=TRUE); x<-gsub("&lt;","<",x,fixed=TRUE); x<-gsub("&amp;","&",x,fixed=TRUE); x }
wb_old<-loadWorkbook(SRC_XLSX); st1<-readWorkbook(wb_old,sheet="ST1",colNames=FALSE,skipEmptyRows=FALSE)[,1:3]
st1[]<-lapply(st1,function(col) vapply(col,dehtml,character(1)))
st1[]<-lapply(st1,function(col) gsub("†","",col))   # drop dagger footnote marker -> plain N/A
addWorksheet(wb,"ST1"); writeData(wb,"ST1",st1,startRow=1,startCol=1,colNames=FALSE)
addStyle(wb,"ST1",title_style,rows=1,cols=1); addStyle(wb,"ST1",header_style,rows=2,cols=1:3,gridExpand=TRUE)
for(r in 3:nrow(st1)){ c2<-st1[r,2];c3<-st1[r,3]
  sec<-!is.na(st1[r,1])&&(is.na(c2)||c2=="")&&(is.na(c3)||c3=="")
  addStyle(wb,"ST1",if(sec)section_style else data_style,rows=r,cols=1:3,gridExpand=TRUE,stack=TRUE)}
setColWidths(wb,"ST1",cols=1:3,widths="auto"); cat("  ST1: verbatim\n")

## ST2 = serum cytokine (old ST8); drop RA-vs-Healthy & Sarc-vs-Healthy P
d<-read.csv(file.path(TBL,"S3_Healthy_Serum_Comparison.csv"),check.names=FALSE)
d[["RA_vs_Sarc_p"]]<-sgf(d[["RA_vs_Sarc_p"]],3)
add_sheet("ST2","Supplementary Table 2. Serum cytokine comparison: RA, sarcoidosis, and healthy controls.",
  c("Cytokine","Healthy\nN","Healthy\nmedian [95% CI]","Sarcoidosis\nN","Sarcoidosis\nmedian [IQR]",
    "RA\nN","RA\nmedian [IQR]","RA vs Sarc\nP"),
  d[,c("Cytokine","Healthy_n","Healthy_median_95CI","Sarcoidosis_n","Sarcoidosis_median_IQR","RA_n","RA_median_IQR","RA_vs_Sarc_p")])

## ST3 = clinical by infection (old ST2)
d<-read.csv(file.path(TBL,"Infection_Clinical_Characteristics.csv"),check.names=FALSE,colClasses="character")
add_sheet("ST3","Supplementary Table 3. Clinical characteristics by prospective infection status (RA, n = 24).",
  c("Variable","Infection (+)\n(n = 6)",paste0("Infection (",MINUS,")\n(n = 18)"),"P value"),
  d[,c("Variable","InfPos","InfNeg","P")])

## ST4 = single biomarker ROC (old ST3, 23 rows desc AUC)
d<-read.csv(file.path(TBL,"Infection_SingleBiomarker_ROC.csv"),check.names=FALSE); d<-d[order(-d$AUC),]
d$AUC<-rnd(d$AUC,3); d$CI_lo<-rnd(d$CI_lo,3); d$CI_hi<-rnd(d$CI_hi,3)
add_sheet("ST4","Supplementary Table 4. Single biomarker ROC analysis for infection prediction.",
  c("Biomarker","AUC","95% CI\nlower","95% CI\nupper","Direction","N"),
  d[,c("Biomarker","AUC","CI_lo","CI_hi","Direction","N")])

## ST5 = microbiome (old ST4, all, asc P)
d<-read.csv(file.path(TBL,"Microbiome_Genus_Infection.csv"),check.names=FALSE); d<-d[order(d$P_value),]
d$Mean_InfPos<-rnd(d$Mean_InfPos,4); d$Mean_InfNeg<-rnd(d$Mean_InfNeg,4); d$Log2FC<-rnd(d$Log2FC,2)
d$P_value<-sgf(d$P_value,3); d$P_adjusted<-sgf(d$P_adjusted,3)
add_sheet("ST5","Supplementary Table 5. BAL microbiome genus-level comparison by prospective infection status (RA, n = 24).",
  c("Genus","Mean (Infection +)",paste0("Mean (Infection ",MINUS,")"),paste0("Log",SUB2,"FC"),"P value","FDR-adjusted P"),
  d[,c("Genus","Mean_InfPos","Mean_InfNeg","Log2FC","P_value","P_adjusted")])

## ST6 = clinical vs omics (old ST5, clean 17)
d<-read.csv(file.path(TBL,"Clinical_vs_Multiomics_Comparison.csv"),check.names=FALSE)
d$Inf_AUC<-rnd(d$Inf_AUC,3); d$Inf_P<-rnd(d$Inf_P,3); d$Prog_AUC<-rnd(d$Prog_AUC,3); d$Prog_P<-rnd(d$Prog_P,3)
add_sheet("ST6","Supplementary Table 6. Clinical markers vs multi-omics biomarkers comparison.",
  c("Category","Variable","Infection\nAUC","Infection\nP value","Progression\nAUC","Progression\nP value"),
  d[,c("Category","Variable","Inf_AUC","Inf_P","Prog_AUC","Prog_P")])

## ST7 = ILD strat (old ST6, all 108, asc P)
d<-read.csv(file.path(TBL,"ST4_ILD_stratification_updated.csv"),check.names=FALSE); d<-d[order(d$P_value),]
d$Mean_ILD_pos<-rnd(d$Mean_ILD_pos,2); d$Mean_ILD_neg<-rnd(d$Mean_ILD_neg,2)
d$P_value<-sgf(d$P_value,3); d$P_adjusted<-sgf(d$P_adjusted,3)
add_sheet("ST7","Supplementary Table 7. BAL cellular and serum cytokine profiles by ILD status (ILD+, n = 11; ILD−, n = 13).",
  c("Variable","Mean\nILD (+)",paste0("Mean\nILD (",MINUS,")"),"N\nILD (+)",paste0("N\nILD (",MINUS,")"),"P value","FDR-adjusted\nP value"),
  d[,c("Variable","Mean_ILD_pos","Mean_ILD_neg","N_ILD_pos","N_ILD_neg","P_value","P_adjusted")])

## ST8 = CT progression correlations (old ST7, all 432, asc p)
d<-read.csv(file.path(TBL,"CT_Delta_vs_NonGene_RAonly.csv"),check.names=FALSE); d<-d[order(d$p),]
d$rho<-rnd(d$rho,3); d$p<-sgf(d$p,3)
add_sheet("ST8","Supplementary Table 8. CT progression correlations with non-gene biomarkers (RA only, n = 24).",
  c("CT variable","Biomarker",paste0("Spearman ",RHO),"P value","N"),
  d[,c("CT_Var","Biomarker","rho","p","n")])

## ST9 = medication (old ST9; IDD->Factor 1, IIV->Factor 2)
d<-read.csv(file.path(TBL,"Medication_Confounding_RAfactor.csv"),check.names=FALSE)
d$Median_on<-rnd(d$Median_on,2); d$Median_off<-rnd(d$Median_off,2); d$P<-rnd(d$P,4)
add_sheet("ST9","Supplementary Table 9. Medication confounding analysis for RA-specific outcome predictors (RA Factor 1 and BAL Th17.1).",
  c("Medication","Variable","N (on)","N (off)","Median (on)","Median (off)","P value"),
  d[,c("Medication","Variable","N_on","N_off","Median_on","Median_off","P")])

saveWorkbook(wb,OUT,overwrite=TRUE); cat(sprintf("\nSaved: %s\n",OUT))
