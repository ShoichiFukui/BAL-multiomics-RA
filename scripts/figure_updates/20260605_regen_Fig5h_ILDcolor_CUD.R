# Fig5h: color points by baseline ILD status (CUD: ILD− blue, ILD+ orange) + legend.
# Single overall regression line + rho/p (overall Spearman) retained.
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({library(tidyverse);library(ggplot2)})
EXCL<-character(0); panel_dir<-"results/panels"
CUD_BLUE<-"#0041FF"; CUD_ORANGE<-"#FF9900"
load("results/RA_ILD_Workspace.RData")
e<-new.env(); load("results/MOFA2_RA_only_Results.RData",envir=e); ra_factors<-get("ra_factors",e)
ct<-read.csv("results/tables/master_data_with_CT_all.csv"); ct<-ct[!duplicated(ct$Sample_ID),]
ids<-rownames(ra_factors)
dhl<-as.numeric(ct$CT_Delta_Healthy_Lung_pct[match(ids,ct$Sample_ID)])
ild<-master_data$ILD_Group[match(ids,master_data$Sample_ID)]
ild<-ifelse(ild=="ILD_Positive","ILD+",ifelse(ild=="ILD_Negative","ILD−",NA))
f1<-as.numeric(ra_factors[,1]); ok<-!is.na(f1)&!is.na(dhl)
sp<-suppressWarnings(cor.test(f1[ok],dhl[ok],method="spearman",exact=TRUE))
cat(sprintf("rho=%.3f p=%.4f | ILD+ n=%d ILD- n=%d\n",sp$estimate,sp$p.value,sum(ild[ok]=="ILD+",na.rm=T),sum(ild[ok]=="ILD−",na.rm=T)))
df<-data.frame(Fval=f1[ok],DHL=dhl[ok],ILD=factor(ild[ok],levels=c("ILD−","ILD+")))
tn<-function(bs=8){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        legend.text=element_text(size=rel(0.8)),legend.title=element_text(size=rel(0.85),face="bold"),
        legend.key.size=unit(3,"mm"),legend.background=element_rect(fill=alpha("white",0.95),color=NA),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,3,2,3),"mm"))}
p<-ggplot(df,aes(Fval,DHL))+
  geom_smooth(method="lm",se=TRUE,color="grey30",linewidth=0.4,fill="grey90",formula=y~x)+
  geom_hline(yintercept=0,linetype="dotted",color="grey60",linewidth=0.2)+
  geom_point(aes(color=ILD),size=1.8,alpha=0.85)+
  scale_color_manual(values=c("ILD−"=CUD_BLUE,"ILD+"=CUD_ORANGE),name="ILD status")+
  annotate("text",x=-Inf,y=Inf,hjust=-0.05,vjust=1.3,size=3.3,
           label=paste0("ρ=",sprintf("%.3f",sp$estimate),"\np=",sprintf("%.4f",sp$p.value)))+
  labs(x="RA Factor 1",y="Δ Healthy Lung (%)")+
  tn(8)+theme(legend.position=c(0.83,0.22))
ggsave(file.path(panel_dir,"Fig5h_RA_factor_CT.pdf"),p,width=80,height=70,units="mm",dpi=300,device=cairo_pdf)
ggsave(file.path(panel_dir,"Fig5h_RA_factor_CT.png"),p,width=80,height=70,units="mm",dpi=300,bg="white")
cat("Saved Fig5h_RA_factor_CT (ILD-colored, CUD)\n")
