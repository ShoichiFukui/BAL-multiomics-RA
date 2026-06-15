# Remove the green "Nested LOOCV (AUC = 0.962)" dashed line + text from Fig5c
# (redundant with Fig2l; overlapped the legend). CUD colors retained.
# Adjust BASEDIR to your local clone of the repository.
BASEDIR <- Sys.getenv("HOSOGAYA_BAL_DIR", unset = getwd())
setwd(BASEDIR)
suppressPackageStartupMessages({library(tidyverse);library(pROC)})
CUD_ORANGE<-"#FF9900"; CUD_BLUE<-"#0041FF"; panel_dir<-"results/panels"
load("results/RA_ILD_Workspace.RData"); master_all<-master_data
tn<-function(bs=8){theme_classic(base_size=bs,base_family="Helvetica")%+replace%
  theme(text=element_text(color="black"),axis.text=element_text(size=rel(0.9),color="black"),
        axis.title=element_text(size=rel(0.95),face="bold"),
        plot.title=element_text(size=rel(1.05),face="bold",hjust=0,margin=ggplot2::margin(0,0,3,0)),
        panel.border=element_rect(fill=NA,color="black",linewidth=0.33),axis.line=element_blank(),
        plot.margin=unit(c(2,3,2,3),"mm"))}
add_ep<-function(spec,sens){ord<-order(spec,sens);spec<-spec[ord];sens<-sens[ord]
  if(spec[1]!=0||sens[1]!=0){spec<-c(0,spec);sens<-c(0,sens)}
  if(tail(spec,1)!=1||tail(sens,1)!=1){spec<-c(spec,1);sens<-c(sens,1)};list(spec=spec,sens=sens)}
e6<-new.env(); load("results/MOFA2_6views_Results.RData",envir=e6); factors6<-get("factors",e6)
y<-ifelse(grepl("^(KYC|Sarcoidosis)",rownames(factors6)),0,1)
roc_f1<-roc(y,factors6[,1],quiet=TRUE,direction="auto")
serum<-master_all[match(rownames(factors6),master_all$Sample_ID),]
best<-NULL;ba<-0
for(cc in grep("^Serum_",colnames(serum),value=TRUE)){v<-as.numeric(serum[[cc]])
  if(sum(!is.na(v))>=30){r<-tryCatch(roc(y,v,quiet=TRUE,direction="auto"),error=function(x)NULL)
    if(!is.null(r)&&as.numeric(auc(r))>ba){ba<-as.numeric(auc(r));best<-cc}}}
roc_best<-roc(y,as.numeric(serum[[best]]),quiet=TRUE,direction="auto")
cf<-gsub("^Serum_","",best);cf<-gsub("^GCSF$","G-CSF",cf);cf<-gsub("^IL(\\d)","IL-\\1",cf)
ep1<-add_ep(1-roc_f1$specificities,roc_f1$sensitivities);ep2<-add_ep(1-roc_best$specificities,roc_best$sensitivities)
df<-rbind(data.frame(Spec=ep1$spec,Sens=ep1$sens,Model=sprintf("Factor 1 (AUC = %.3f)",auc(roc_f1))),
          data.frame(Spec=ep2$spec,Sens=ep2$sens,Model=sprintf("%s (AUC = %.3f)",cf,auc(roc_best))))
df$Model<-factor(df$Model,levels=unique(df$Model))
p<-ggplot(df,aes(Spec,Sens,color=Model))+geom_step(linewidth=0.6,direction="vh")+
  geom_abline(slope=1,intercept=0,linetype="dotted",color="grey60",linewidth=0.3)+
  scale_color_manual(values=c(CUD_BLUE,CUD_ORANGE))+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="on")+
  labs(x="1 – Specificity",y="Sensitivity",color="")+ggtitle("RA vs Sarcoidosis")+
  tn(8)+theme(legend.position=c(0.62,0.18),legend.background=element_rect(fill=alpha("white",0.95),color=NA),
              legend.text=element_text(size=6.5),legend.key.height=unit(3,"mm"))
ggsave(file.path(panel_dir,"Fig5c_diseaseROC_noacronym.pdf"),p,width=80,height=75,units="mm",dpi=300,device=cairo_pdf)
ggsave(file.path(panel_dir,"Fig5c_diseaseROC_noacronym.png"),p,width=80,height=75,units="mm",dpi=300,bg="white")
cat(sprintf("5c regenerated WITHOUT LOOCV annot: F1=%.3f, %s=%.3f\n",auc(roc_f1),cf,auc(roc_best)))
