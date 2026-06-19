# Regenerate Supplementary Fig 4 (medication confounding) with the MOFA factor
# axis labels as "Factor 1" / "Factor 2".
suppressPackageStartupMessages({library(tidyverse);library(ggplot2);library(patchwork)})
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR); panel_dir<-"results/panels"
EXCL<-character(0)
tn<-function(bs=7){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.0),face="bold",hjust=0.5,margin=ggplot2::margin(0,0,2,0)),
        legend.text=element_text(size=rel(0.8)),legend.key.size=unit(3,"mm"),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,2,2,2),"mm"))}
sv<-function(p,name,w,h){ggsave(file.path(panel_dir,paste0(name,".pdf")),p,width=w,height=h,units="mm",dpi=300,device=cairo_pdf)
  ggsave(file.path(panel_dir,paste0(name,".png")),p,width=w,height=h,units="mm",dpi=300,bg="white");cat("Saved",name,"\n")}
fmt_p<-function(p){if(p<0.001)return("p < 0.001");sprintf("p = %.2f",p)}
col_ra<-"#C0392B"; col_sarc<-"#2980B9"
load("results/RA_ILD_Workspace.RData"); load("results/MOFA2_6views_Results.RData")
ct_data<-read.csv("results/tables/master_data_with_CT_all.csv"); ct_unique<-ct_data[!duplicated(ct_data$Sample_ID),]
ra_samples<-setdiff(rownames(factors)[grepl("^(KY[0-9]|RA[0-9])",rownames(factors))],EXCL)
med_cols<-c("steroid.1.0.","MTX","NSAIDS"); med_labels<-c("GC","MTX","NSAIDs")
factor_names<-c("Factor 1","Factor 2")
plots_s4<-list(); pi<-1
for(fi in 1:2){ fn<-factor_names[fi]; fac_vals<-factors[ra_samples,fi]
  for(mi in seq_along(med_cols)){ mc<-med_cols[mi]; ml<-med_labels[mi]
    med_vals<-ct_unique[[mc]][match(ra_samples,ct_unique$Sample_ID)]; ok<-!is.na(fac_vals)&!is.na(med_vals)
    on_vals<-fac_vals[ok&med_vals==1]; off_vals<-fac_vals[ok&med_vals==0]
    if(length(on_vals)<2||length(off_vals)<2) next
    wt<-wilcox.test(on_vals,off_vals,exact=TRUE); p_lab<-fmt_p(wt$p.value)
    df<-data.frame(Value=c(on_vals,off_vals),
      Med=factor(c(rep(paste0(ml," (+)"),length(on_vals)),rep(paste0(ml," (−)"),length(off_vals))),
                 levels=c(paste0(ml," (−)"),paste0(ml," (+)"))))
    yr<-range(df$Value); ymax<-yr[2]+diff(yr)*0.25
    p<-ggplot(df,aes(Med,Value,fill=Med))+geom_boxplot(outlier.shape=NA,width=0.55,linewidth=0.3,alpha=0.7)+
      geom_jitter(width=0.12,size=0.8,alpha=0.6,show.legend=FALSE)+scale_fill_manual(values=c(col_sarc,col_ra))+
      labs(x="",y=fn,title=ml)+annotate("text",x=1.5,y=ymax,label=p_lab,size=2.2,fontface="italic")+
      coord_cartesian(ylim=c(min(yr[1],0),ymax*1.05),clip="off")+tn(7)+
      theme(legend.position="none",axis.text.x=element_text(size=6))
    plots_s4[[pi]]<-p; pi<-pi+1 }}
p_s4<-wrap_plots(plots_s4,ncol=3)+plot_annotation(theme=theme(plot.margin=ggplot2::margin(2,2,2,2)))
sv(p_s4,"FigS4_medication",160,100)
