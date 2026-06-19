# Medication-confounding for the RA-specific outcome predictors used in Fig 5:
# RA-only MOFA Factor 1 (ra_factors[,1]) and BAL Th17.1, vs GC/MTX/NSAID.
# Writes Medication_Confounding_RAfactor.csv and regenerates FigS4_medication
# (2 rows x 3 medications).
suppressPackageStartupMessages({library(tidyverse);library(ggplot2);library(patchwork)})
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR); panel_dir<-"results/panels"
EXCL<-character(0)
tn<-function(bs=7){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.0),face="bold",hjust=0.5,margin=ggplot2::margin(0,0,2,0)),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,2,2,2),"mm"))}
sv<-function(p,name,w,h){ggsave(file.path(panel_dir,paste0(name,".pdf")),p,width=w,height=h,units="mm",dpi=300,device=cairo_pdf)
  ggsave(file.path(panel_dir,paste0(name,".png")),p,width=w,height=h,units="mm",dpi=300,bg="white");cat("Saved",name,"\n")}
fmt_p<-function(p){if(p<0.001)return("p < 0.001");sprintf("p = %.2f",p)}
col_ra<-"#C0392B"; col_sarc<-"#2980B9"
load("results/RA_ILD_Workspace.RData")
e<-new.env(); load("results/MOFA2_RA_only_Results.RData",envir=e); ra_factors<-get("ra_factors",e)
ct<-read.csv("results/tables/master_data_with_CT_all.csv"); ct<-ct[!duplicated(ct$Sample_ID),]
ids<-rownames(ra_factors)
md<-master_data[!master_data$Sample_ID %in% EXCL,]
vars<-list(list(nm="RA Factor 1", lab="RA Factor 1", v=as.numeric(ra_factors[,1])),
           list(nm="BAL Th17.1", lab="BAL Th17.1 (%)", v=as.numeric(md$BALF_Th17.1_CCR6pos_CXCR3pos[match(ids,md$Sample_ID)])))
med_cols<-c("steroid.1.0.","MTX","NSAIDS"); med_lab<-c("GC","MTX","NSAIDs")
rows<-list(); plots<-list(); pi<-1
for(vv in vars){ for(k in seq_along(med_cols)){
  mv<-ct[[med_cols[k]]][match(ids,ct$Sample_ID)]; ok<-!is.na(vv$v)&!is.na(mv)
  on<-vv$v[ok&mv==1]; off<-vv$v[ok&mv==0]
  w<-suppressWarnings(wilcox.test(on,off,exact=TRUE))
  rows[[length(rows)+1]]<-data.frame(Medication=med_lab[k],Variable=vv$nm,N_on=length(on),N_off=length(off),
                                     Median_on=round(median(on),2),Median_off=round(median(off),2),P=round(w$p.value,4))
  df<-data.frame(Value=c(on,off),Med=factor(c(rep(paste0(med_lab[k]," (+)"),length(on)),rep(paste0(med_lab[k]," (−)"),length(off))),
                 levels=c(paste0(med_lab[k]," (−)"),paste0(med_lab[k]," (+)"))))
  yr<-range(df$Value); ymax<-yr[2]+diff(yr)*0.25
  plots[[pi]]<-ggplot(df,aes(Med,Value,fill=Med))+geom_boxplot(outlier.shape=NA,width=0.55,linewidth=0.3,alpha=0.7)+
    geom_jitter(width=0.12,size=0.8,alpha=0.6,show.legend=FALSE)+scale_fill_manual(values=c(col_sarc,col_ra))+
    labs(x="",y=vv$lab,title=med_lab[k])+annotate("text",x=1.5,y=ymax,label=fmt_p(w$p.value),size=2.2,fontface="italic")+
    coord_cartesian(ylim=c(min(yr[1],0),ymax*1.05),clip="off")+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=6))
  pi<-pi+1 }}
res<-do.call(rbind,rows); write.csv(res,"results/tables/Medication_Confounding_RAfactor.csv",row.names=FALSE)
print(res)
p_s4<-wrap_plots(plots,ncol=3)+plot_annotation(theme=theme(plot.margin=ggplot2::margin(2,2,2,2)))
sv(p_s4,"FigS4_medication",160,100)
