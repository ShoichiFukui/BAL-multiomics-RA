# ============================================================================
# Nature Communications — 7 Main Figures (Final with MOFA2 Integration)
# ============================================================================
if (!exists("BASEDIR")) BASEDIR <- getwd()
setwd(BASEDIR)
load("results/RA_ILD_Workspace.RData")
load("results/CT_Multiomics_Results.RData")

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(gridExtra)
  library(readxl); library(pROC); library(effsize); library(randomForest)
  library(grid); library(vegan); library(Seurat)
})
set.seed(42)
fig_dir <- "results/figures_final"
dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)

# n=35 enforcement
EXCL <- character(0)
master_data <- master_data[!master_data$Sample_ID %in% EXCL, ]
if(exists("master_ext")) master_ext <- master_ext[!master_ext$Sample_ID %in% EXCL, ]
stopifnot(nrow(master_data) == 35)

tn <- function(bs=7.5) {
  theme_classic(base_size=bs) %+replace%
    theme(text=element_text(color="black"),
          axis.text=element_text(size=rel(.9),color="black"),
          axis.title=element_text(size=rel(.95),face="bold"),
          plot.title=element_text(size=rel(1.05),face="bold",hjust=0),
          legend.text=element_text(size=rel(.78)),
          legend.key.size=unit(2.5,"mm"),
          legend.background=element_rect(fill=alpha("white",.92),color=NA),
          panel.border=element_rect(fill=NA,color="black",linewidth=.4),
          axis.line=element_blank(),
          plot.margin=unit(c(3,5,3,4),"mm"))
}
CG <- c("RA"="#C0392B","Sarcoidosis"="#2980B9")
CS <- c("Sarcoidosis"="#2980B9","RA-nonILD"="#27AE60","RA-ILD"="#C0392B")
CI <- c("Infection (-)"="#3498DB","Infection (+)"="#E74C3C")
CC <- c(Macrophage="#E74C3C",T_cell="#3498DB",NK="#2ECC71",B_cell="#9B59B6",
        Plasma="#E67E22",Neutrophil="#F1C40F",DC="#8B4513",Epithelial="#E91E63",
        Mast="#95A5A6",Other="#BDC3C7")
col_ds <- c("GSE145926"="#E74C3C","GSE193782"="#3498DB","GSE184735"="#27AE60")

# Helpers
sct <- function(x,y,xl,yl,ttl,col="#C0392B",grp=NULL,gc=NULL,bs=7) {
  v<-!is.na(x)&!is.na(y); sp<-cor.test(x[v],y[v],method="spearman",exact=TRUE)
  if(is.null(grp)){df<-data.frame(X=x[v],Y=y[v]);p<-ggplot(df,aes(X,Y))+geom_point(size=1.5,color=col,alpha=.7)
  }else{df<-data.frame(X=x[v],Y=y[v],G=grp[v]);p<-ggplot(df,aes(X,Y,color=G))+geom_point(size=1.5,alpha=.8)+scale_color_manual(values=gc)}
  p+geom_smooth(method="lm",se=T,color="grey30",linewidth=.35,fill="grey90",formula=y~x,inherit.aes=F,aes(x=X,y=Y),data=df)+
    geom_hline(yintercept=0,linetype="dotted",color="grey60",linewidth=.2)+
    annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.3,size=2,
             label=sprintf("rho=%+.3f\np=%.4f (n=%d)",sp$estimate,sp$p.value,sum(v)))+
    labs(x=xl,y=yl,title=ttl)+tn(bs)+theme(legend.position="none")
}
bxi <- function(col,yl,ttl,rd,bs=7) {
  v<-as.numeric(rd[[col]]);ig<-rd$IG;ok<-!is.na(v)&!is.na(ig);df<-data.frame(V=v[ok],IG=ig[ok])
  wt<-wilcox.test(df$V[df$IG=="Infection (+)"],df$V[df$IG=="Infection (-)"], exact=TRUE)
  ymx<-max(df$V,na.rm=T);ymn<-min(df$V,na.rm=T);ysp<-(ymx-ymn)*0.12
  ggplot(df,aes(IG,V,fill=IG))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+
    geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+scale_fill_manual(values=CI)+
    labs(x="",y=yl,title=ttl)+
    annotate("text",x=1.5,y=ymx+ysp,label=sprintf("p=%.3f",wt$p.value),size=2)+
    coord_cartesian(ylim=c(ymn-ysp*.3,ymx+ysp*1.5),clip="off")+
    tn(bs)+theme(legend.position="none",axis.text.x=element_text(size=5.5))
}
bvp <- function(bc,pc,nm,ttl) {
  ba<-as.numeric(master_data[[bc]]);pb<-as.numeric(master_data[[pc]])
  gg<-ifelse(master_data$Sample_Group=="Control","Sarcoidosis","RA")
  v<-!is.na(ba)&!is.na(pb);sp<-cor.test(ba[v],pb[v],method="spearman",exact=TRUE)
  df<-data.frame(X=pb[v],Y=ba[v],G=gg[v])
  ggplot(df,aes(X,Y,color=G))+geom_point(size=1.5,alpha=.8)+
    geom_abline(slope=1,intercept=0,linetype="dotted",color="grey60")+scale_color_manual(values=CG)+
    annotate("text",x=Inf,y=-Inf,hjust=1.05,vjust=-.3,size=1.9,label=sprintf("rho=%.3f\np=%.4f",sp$estimate,sp$p.value))+
    labs(x=paste("PB",nm,"(%)"),y=paste("BAL",nm,"(%)"),title=ttl,color="")+
    tn(7)+theme(legend.position="bottom",legend.direction="horizontal")
}

samples <- common_samples
pca <- prcomp(t(expr_matrix),scale.=T); ve <- summary(pca)$importance[2,1:5]*100
pd <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2],SID=rownames(pca$x))%>%
  mutate(Group=ifelse(grepl("^(KYC|Sarcoidosis)",SID),"Sarcoidosis","RA"),
         Sub=recode(as.character(meta$Subgroup[match(SID,meta$Sample_ID)]),"Control"="Sarcoidosis","RA_nonILD"="RA-nonILD","RA_ILD"="RA-ILD"),
         InfSt=recode(as.character(meta$Infection_Group[match(SID,meta$Sample_ID)]),"Infection_Negative"="Infection (-)","Infection_Positive"="Infection (+)"))
xl<-sprintf("PC1 (%.1f%%)",ve[1]); yl<-sprintf("PC2 (%.1f%%)",ve[2])
deg <- read.csv("results/tables/DEG_ShrunkLFC.csv")
ggo <- read.csv("results/tables/GSEA_GO_BP.csv")
gkk <- read.csv("results/tables/GSEA_KEGG.csv")
es <- read.csv("results/tables/Effect_sizes_RA_vs_Control.csv")
inf_comp <- read.csv("results/tables/Infection_AllOmics_Comparison.csv")
bio_roc <- read.csv("results/tables/Infection_SingleBiomarker_ROC.csv")
deg_inf <- read.csv("results/tables/DEG_Infection_pos_vs_neg.csv")
bp_cell <- read.csv("output_v3_BayesPrism/celltype_proportions_cellFraction.csv",row.names=1)
bp_theta <- read.csv("output_v3_BayesPrism/celltype_proportions_BayesPrism.csv",row.names=1)
fcm <- read_excel("FCM_integrated_data_transformed.xlsx",sheet="BALF_analysis");colnames(fcm)[1]<-"Sample_ID"
ra_md <- master_data[master_data$Sample_Group=="RA"&!is.na(master_data$respiratory_infection),]
ra_md$IG <- ifelse(ra_md$respiratory_infection==1,"Infection (+)","Infection (-)")
grp_all <- recode(as.character(master_ext$Subgroup),"Control"="Sarcoidosis","RA_nonILD"="RA-nonILD","RA_ILD"="RA-ILD")
rate_hl <- as.numeric(master_ext$CT_Rate_Healthy_Lung_pct)

# ============================================================================
# Fig 1: Study Design + scRNA-seq Reference (8 panels)
# ============================================================================
cat("Fig 1...\n")
f1a<-ggplot(pd,aes(PC1,PC2,color=Group))+geom_point(size=1.8,alpha=.8)+scale_color_manual(values=CG)+
  labs(x=xl,y=yl,title="a",color="")+tn()+theme(legend.position="bottom",legend.direction="horizontal")
f1b<-ggplot(pd%>%filter(!is.na(Sub)),aes(PC1,PC2,color=Sub))+geom_point(size=1.8,alpha=.8)+
  scale_color_manual(values=CS)+labs(x=xl,y=yl,title="b",color="")+tn()+theme(legend.position="bottom",legend.direction="horizontal")

ref<-readRDS("output_v3_BayesPrism/BAL_reference_author_annotated.rds")
uc<-Embeddings(ref,"umap");ms<-ref@meta.data
ds_orig<-data.frame(U1=uc[,1],U2=uc[,2],DS=ms$dataset,CT=ms$celltype)

# Compute dotplot BEFORE shuffling (indices must match expression matrix)
mks<-intersect(c("CD68","MARCO","FABP4","FCN1","CD3D","CD3E","NKG7","GNLY","CD79A","MS4A1","JCHAIN","MZB1","FCGR3B","FCER1A","EPCAM","KRT18"),rownames(ref))
ed<-GetAssayData(ref,layer="data");dd<-data.frame()
for(ct in unique(ds_orig$CT)){cl<-which(ds_orig$CT==ct);if(length(cl)<10)next;for(mk in mks){vv<-ed[mk,cl];dd<-rbind(dd,data.frame(CT=ct,M=mk,Pct=sum(vv>0)/length(vv)*100,Avg=mean(vv)))}}

# Shuffle for UMAP plotting only
set.seed(42);ds<-ds_orig[sample(nrow(ds_orig)),]

f1c<-ggplot(ds,aes(U1,U2,color=DS))+geom_point(size=.06,alpha=.2)+
  scale_color_manual(values=col_ds,labels=c("GSE145926"="COVID-19","GSE193782"="Healthy","GSE184735"="Sarcoidosis"))+
  guides(color=guide_legend(override.aes=list(size=2,alpha=1)))+
  labs(x="UMAP 1",y="UMAP 2",title="c",color="")+tn(7)+theme(legend.position="bottom",legend.direction="horizontal")
f1d<-ggplot(ds,aes(U1,U2,color=CT))+geom_point(size=.06,alpha=.2)+scale_color_manual(values=CC)+
  guides(color=guide_legend(override.aes=list(size=1.5,alpha=1),ncol=5))+
  labs(x="UMAP 1",y="UMAP 2",title="d",color="")+tn(7)+theme(legend.position="bottom")
dd$M<-factor(dd$M,levels=rev(mks))
f1e<-ggplot(dd,aes(CT,M,size=Pct,color=Avg))+geom_point()+
  scale_size_continuous(range=c(.3,3),name="%Exp")+scale_color_gradient(low="grey90",high="#C0392B",name="Avg")+
  labs(x="",y="",title="e")+tn(6)+theme(axis.text.x=element_text(angle=45,hjust=1,size=5),axis.text.y=element_text(size=5,face="italic"))

cm<-ds%>%group_by(DS,CT)%>%summarise(N=n(),.groups="drop")%>%group_by(DS)%>%mutate(P=N/sum(N))
cm$DS<-recode(cm$DS,"GSE145926"="COVID-19","GSE193782"="Healthy","GSE184735"="Sarcoidosis")
f1f<-ggplot(cm,aes(DS,P,fill=CT))+geom_col(width=.6)+scale_fill_manual(values=CC)+
  labs(x="",y="Proportion",title="f",fill="")+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=7))
rm(ref,ed);gc()

pd_inf<-pd%>%filter(!is.na(InfSt))
cen<-pd_inf%>%group_by(InfSt)%>%summarise(PC1m=mean(PC1),PC2m=mean(PC2),.groups="drop")
expr_inf_m<-t(expr_matrix[,pd_inf$SID]);dist_inf<-dist(expr_inf_m)
set.seed(42)
suppressMessages(perm_inf<-adonis2(dist_inf~pd_inf$InfSt,permutations=999))
f1g<-ggplot(pd_inf,aes(PC1,PC2,color=InfSt))+stat_ellipse(level=.95,linewidth=.5)+
  geom_point(size=1.8,alpha=.8)+geom_point(data=cen,aes(PC1m,PC2m),shape=4,size=2.4,stroke=1.0,show.legend=F)+
  scale_color_manual(values=c("Infection (-)"="#3498DB","Infection (+)"="#E74C3C"))+
  annotate("text",x=Inf,y=Inf,hjust=1.05,vjust=1.3,size=2,
           label=sprintf("PERMANOVA\nR2=%.3f, p=%.3f",perm_inf$R2[1],perm_inf$`Pr(>F)`[1]))+
  labs(x=xl,y=yl,title="g",color="")+tn()+theme(legend.position="bottom",legend.direction="horizontal")
f1h<-ggplot(data.frame(PC=factor(paste0("PC",1:5),levels=paste0("PC",1:5)),V=ve[1:5]),aes(PC,V))+
  geom_col(fill="#34495E",width=.5)+labs(x="",y="Variance(%)",title="h")+tn()

ggsave(file.path(fig_dir,"Fig1_StudyDesign.pdf"),
       arrangeGrob(f1a,f1b,f1c,f1d,f1e,f1f,f1g,f1h,
                   layout_matrix=rbind(c(1,2,3,4),c(5,5,6,7),c(NA,NA,NA,8)),heights=c(1,1.1,.55)),
       width=210,height=200,units="mm",dpi=300)
cat("  done\n")

# ============================================================================
# Fig 2: RA vs Sarcoidosis (10 panels)
# ============================================================================
cat("Fig 2...\n")
deg$sig<-case_when(deg$padj<.05&deg$log2FoldChange>.585~"Up",deg$padj<.05&deg$log2FoldChange< -.585~"Down",deg$padj<.05~"Sig",TRUE~"NS")
deg$lab<-ifelse(!is.na(deg$padj)&deg$padj<.05&!grepl("^ENSG",deg$Gene_Symbol),deg$Gene_Symbol,NA)
f2a<-ggplot(deg,aes(log2FoldChange,-log10(padj),color=sig))+
  geom_point(data=deg%>%filter(sig=="NS"),size=.2,alpha=.15)+geom_point(data=deg%>%filter(sig!="NS"),size=1.5,alpha=.9)+
  geom_text(data=deg%>%filter(!is.na(lab)),aes(label=lab),hjust=-.1,size=1.9,fontface="italic",color="black",show.legend=F)+
  scale_color_manual(values=c(Up="#C0392B",Down="#2980B9",Sig="#E67E22",NS="grey85"))+
  geom_hline(yintercept=-log10(.05),linetype="dashed",color="grey50",linewidth=.2)+
  geom_vline(xintercept=c(-.585,.585),linetype="dashed",color="grey50",linewidth=.2)+
  xlim(-2,2)+labs(x="Shrunken log2FC",y="-log10(padj)",title="a",color="")+tn(7)+theme(legend.position="right")

tgo<-head(ggo[order(ggo$p.adjust),],10)%>%mutate(D=str_trunc(Description,44));tgo<-tgo[order(tgo$NES),];tgo$D<-fct_inorder(tgo$D)
f2b<-ggplot(tgo,aes(NES,D,fill=NES>0))+geom_col(width=.6,show.legend=F)+scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2980B9"))+
  geom_vline(xintercept=0,linewidth=.2)+labs(x="NES",y="",title="b  GO BP")+tn(6)+theme(axis.text.y=element_text(size=5.5))

tkk<-head(gkk[order(gkk$p.adjust),],10)%>%mutate(D=str_trunc(Description,44));tkk<-tkk[order(tkk$NES),];tkk$D<-fct_inorder(tkk$D)
f2c<-ggplot(tkk,aes(NES,D,fill=NES>0))+geom_col(width=.6,show.legend=F)+scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2980B9"))+
  geom_vline(xintercept=0,linewidth=.2)+labs(x="NES",y="",title="c  KEGG")+tn(6)+theme(axis.text.y=element_text(size=5.5))

gm<-gsva_scores[,common_samples];so<-common_samples[order(meta$Sample_Group[match(common_samples,meta$Sample_ID)])]
glg<-as.data.frame(gm)%>%rownames_to_column("PW")%>%pivot_longer(-PW,names_to="S",values_to="Sc");glg$S<-factor(glg$S,levels=so)
f2d<-ggplot(glg,aes(S,PW,fill=Sc))+geom_tile()+scale_fill_gradient2(low="#2980B9",mid="white",high="#C0392B",midpoint=0,name="Score")+
  labs(x="",y="",title="d  GSVA")+tn(5)+theme(axis.text.x=element_text(angle=55,hjust=1,size=3),axis.text.y=element_text(size=5))

# ROC
lo_s<-intersect(intersect(rownames(expr_mat),rownames(cyto_mat)),rownames(facs_mat))
lo_d<-data.frame(expr_mat[lo_s,1:min(50,ncol(expr_mat))],cyto_mat[lo_s,],facs_mat[lo_s,],Group=meta$Sample_Group[match(lo_s,meta$Sample_ID)])
lo_d<-lo_d[!is.na(lo_d$Group),];lo_d$Group<-factor(lo_d$Group)
for(cc in setdiff(colnames(lo_d),"Group"))lo_d[[cc]][is.na(lo_d[[cc]])|is.infinite(lo_d[[cc]])]<-median(lo_d[[cc]],na.rm=T)
fvv<-apply(lo_d[,setdiff(colnames(lo_d),"Group"),drop=F],2,var,na.rm=T);lo_d<-lo_d[,c(names(fvv[!is.na(fvv)&fvv>0]),"Group")]
nl<-nrow(lo_d);np<-matrix(NA,nl,2,dimnames=list(NULL,levels(lo_d$Group)))
for(i in 1:nl){tr<-lo_d[-i,];ft<-setdiff(colnames(tr),"Group");pv<-apply(tr[,ft],2,function(x)tryCatch(wilcox.test(x~tr$Group, exact=TRUE)$p.value,error=function(e)1))
  tp<-names(sort(pv))[1:min(50,length(pv))];rf<-randomForest(x=tr[,tp],y=tr$Group,ntree=500);np[i,]<-stats::predict(rf,lo_d[i,tp],type="prob")}
roc1<-roc(lo_d$Group,np[,"RA"],quiet=T);ci1<-ci.auc(roc1,method="bootstrap",boot.n=2000,quiet=T)
rd1<-data.frame(S=roc1$sensitivities,Sp=1-roc1$specificities)
rd1<-rd1[order(rd1$Sp,rd1$S),]
f2e<-ggplot(rd1,aes(Sp,S))+geom_step(color="#C0392B",linewidth=.8,direction="vh")+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60")+
  annotate("text",x=.55,y=.1,size=2.2,hjust=0,label=sprintf("AUC=%.3f\n95%%CI:%.3f-%.3f",as.numeric(auc(roc1)),ci1[1],ci1[3]))+
  labs(x="1-Specificity",y="Sensitivity",title="e")+coord_equal()+tn(7)

ses<-es%>%filter(!is.na(p_adjusted),p_adjusted<.05)%>%mutate(L=gsub("^BALF_|^Serum_","",Variable),Src=ifelse(grepl("Serum",Category),"Serum","BALF"))%>%
  arrange(cliff_delta)%>%tail(15);ses$L<-fct_inorder(ses$L)
f2f<-ggplot(ses,aes(cliff_delta,L,fill=Src))+geom_col(width=.5)+scale_fill_manual(values=c(Serum="#E74C3C",BALF="#3498DB"))+
  geom_vline(xintercept=0,linewidth=.2)+labs(x="Cliff's delta",y="",title="f",fill="")+tn(6.5)+theme(legend.position="top",axis.text.y=element_text(size=5.5))

ga_pl<-ifelse(grepl("^(KYC|Sarcoidosis)",rownames(bp_theta)),"Sarcoidosis","RA")
dp<-data.frame(P=bp_theta[,"Plasma"]*100,B=bp_theta[,"B_cell"]*100,G=ga_pl);dp$Ratio<-dp$P/(dp$B+dp$P+1e-10)*100
wtp<-wilcox.test(dp$P[dp$G=="RA"],dp$P[dp$G=="Sarcoidosis"], exact=TRUE)
f2g<-ggplot(dp,aes(G,P+1e-6,fill=G))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+
  scale_y_log10()+scale_fill_manual(values=CG)+labs(x="",y="Plasma (%, log10)",title="g")+
  annotate("text",x=1.5,y=max(dp$P)*3,label=sprintf("p=%.3f",wtp$p.value),size=2)+tn(7)+theme(legend.position="none")

ig<-deg[grepl("^IGHG|^IGHM$|^IGKC$|^IGHA[12]$|^JCHAIN$",deg$Gene_Symbol)&!is.na(deg$padj),]%>%arrange(desc(abs(log2FoldChange)))%>%head(7)
ig$Gene_Symbol<-fct_reorder(ig$Gene_Symbol,ig$log2FoldChange);ig$sg<-ifelse(ig$padj<.05,"*","NS")
f2h<-ggplot(ig,aes(log2FoldChange,Gene_Symbol,fill=sg))+geom_col(width=.5)+scale_fill_manual(values=c("*"="#C0392B","NS"="#95A5A6"))+
  geom_vline(xintercept=0,linewidth=.2)+labs(x="log2FC",y="",title="h  Ig genes",fill="")+tn(7)+theme(legend.position="bottom",legend.direction="horizontal")

ggsave(file.path(fig_dir,"Fig2_Molecular.pdf"),
       arrangeGrob(f2a,f2b,f2c,f2d,f2e,f2f,f2g,f2h,layout_matrix=rbind(c(1,1,2,2),c(3,3,4,4),c(5,6,6,7),c(NA,NA,NA,8)),heights=c(1,1,1,.7)),
       width=200,height=260,units="mm",dpi=300)
cat("  done\n")

# ============================================================================
# Fig 3: Deconvolution + BAL vs PB (10 panels)
# ============================================================================
cat("Fig 3...\n")
bp_l<-bp_cell%>%rownames_to_column("S")%>%mutate(G=ifelse(grepl("^(KYC|Sarcoidosis)",S),"Sarcoidosis","RA"))%>%
  mutate(S_display=ifelse(grepl("^(KYC|Sarcoidosis)",S),paste0("Sarcoidosis",gsub("^(KYC0*|Sarcoidosis)","",S)),paste0("RA",gsub("^(KY0*|RA)","",S))))%>%
  pivot_longer(-c(S,G,S_display),names_to="CT",values_to="Fr")%>%mutate(S_display=fct_reorder(S_display,as.numeric(factor(G))))
f3a<-ggplot(bp_l,aes(S_display,Fr,fill=CT))+geom_col(width=.85)+scale_fill_manual(values=CC)+
  labs(x="",y="Proportion",title="a",fill="")+tn(6)+theme(axis.text.x=element_text(angle=55,hjust=1,size=3.5),legend.position="bottom",legend.direction="horizontal",legend.text=element_text(size=5))

cf<-intersect(rownames(bp_cell),fcm$Sample_ID);fi<-match(cf,fcm$Sample_ID);lc<-intersect(c("T_cell","NK","B_cell","Plasma"),colnames(bp_cell))
mkf<-function(bv,fv,lb,tt){v<-!is.na(fv)&!is.na(bv);sp<-cor.test(bv[v],fv[v],method="spearman",exact=TRUE)
  df<-data.frame(F=fv[v],B=bv[v],G=ifelse(grepl("^(KYC|Sarcoidosis)",cf[v]),"Sarcoidosis","RA"))
  ggplot(df,aes(F,B,color=G))+geom_point(size=1.2,alpha=.8)+geom_smooth(method="lm",se=F,color="grey30",linewidth=.3,linetype="dashed",formula=y~x,inherit.aes=F,aes(x=F,y=B))+
    geom_abline(slope=1,intercept=0,linetype="dotted",color="grey60")+scale_color_manual(values=CG)+
    annotate("text",x=Inf,y=-Inf,hjust=1.05,vjust=-.3,size=1.8,label=sprintf("rho=%.3f\np=%.4f",sp$estimate,sp$p.value))+
    labs(x=paste("FCM",lb,"(%)"),y=paste("BP",lb,"(%)"),title=tt)+tn(6.5)+theme(legend.position="none")}
f3b<-mkf(bp_cell[cf,"Macrophage"]*100,as.numeric(fcm$BALF_Macrophage_Percent[fi]),"Mac","b")
f3c<-mkf(rowSums(bp_cell[cf,lc])*100,as.numeric(fcm$BALF_Lymphocyte_Percent[fi]),"Lym","c")
f3d<-mkf(bp_cell[cf,"Neutrophil"]*100,as.numeric(fcm$BALF_Neutrophil_Percent[fi]),"Neu","d")

tc<-c("BALF_Th17_CCR6pos_CXCR3neg","BALF_Th17.1_CCR6pos_CXCR3pos","BALF_Th1_CCR6neg_CXCR3pos","BALF_Th2_CCR6neg_CXCR3neg","BALF_Treg_CD127neg_CD25pos")
tn2<-c("Th17","Th17.1","Th1","Th2","Treg");td<-data.frame()
for(i in seq_along(tc))if(tc[i]%in%colnames(master_data)){v<-as.numeric(master_data[[tc[i]]]);g<-master_data$Sample_Group;ok<-!is.na(v)&!is.na(g)
  td<-rbind(td,data.frame(Sub=tn2[i],V=v[ok],G=ifelse(g[ok]=="Control","Sarcoidosis","RA")))}
td$Sub<-factor(td$Sub,levels=tn2)
# Add p-values for each subset
pv_td<-sapply(tn2,function(s){d<-td[td$Sub==s,];tryCatch(wilcox.test(d$V[d$G=="RA"],d$V[d$G=="Sarcoidosis"], exact=TRUE)$p.value,error=function(e)NA)})
pv_lab<-data.frame(Sub=factor(tn2,levels=tn2),Y=sapply(tn2,function(s)max(td$V[td$Sub==s],na.rm=T)*1.12),
                   Lab=ifelse(pv_td<0.05,sprintf("*p=%.3f",pv_td),"NS"))
f3e<-ggplot(td,aes(Sub,V,fill=G))+geom_boxplot(outlier.shape=NA,width=.5,linewidth=.2)+
  geom_jitter(aes(group=G),position=position_jitterdodge(jitter.width=.08,dodge.width=.5),size=.6,alpha=.4,show.legend=F)+
  scale_fill_manual(values=CG)+
  geom_text(data=pv_lab,aes(Sub,Y,label=Lab),inherit.aes=F,size=1.8,color="grey30")+
  labs(x="",y="% of CD4+",title="e",fill="")+tn(7)+theme(legend.position="right",axis.text.x=element_text(size=6.5))

f3f<-bvp("BALF_Th17.1_CCR6pos_CXCR3pos","PB_Th17.1_CCR6pos_CXCR3pos","Th17.1","f")
f3g<-bvp("BALF_Treg_CD127neg_CD25pos","PB_Treg_CD127neg_CD25pos","Treg","g")
f3h<-bvp("BALF_Th17_CCR6pos_CXCR3neg","PB_Th17_CCR6pos_CXCR3neg","Th17","h")
f3i<-bvp("BALF_Activated_Macrophage_CD86pos_CD14pos","PB_Activated_Macrophage_CD86pos_CD14pos","Act.Mac","i")

ggsave(file.path(fig_dir,"Fig3_Deconv_Compartment.pdf"),
       arrangeGrob(f3a,f3b,f3c,f3d,f3e,f3f,f3g,f3h,f3i,
                   layout_matrix=rbind(c(1,1,2,3,4),c(5,5,6,7,8),c(NA,NA,NA,9,NA)),heights=c(1,1,.7)),
       width=210,height=210,units="mm",dpi=300)
cat("  done\n")

# ============================================================================
# Fig 4: Infection Prediction (10 panels)
# ============================================================================
cat("Fig 4...\n")
top4<-head(bio_roc%>%filter(!grepl("Gene_",Biomarker)),4)
clrs<-c("#C0392B","#2980B9","#27AE60","#8E44AD")
ra_c2<-master_data[master_data$Sample_Group=="RA"&!is.na(master_data$respiratory_infection)&master_data$Sample_ID%in%common_samples,]
y_inf<-factor(ifelse(ra_c2$respiratory_infection==1,"Inf_pos","Inf_neg"),levels=c("Inf_neg","Inf_pos"))
rl<-data.frame()
for(i in 1:nrow(top4)){nm<-top4$Biomarker[i];vv<-NULL
  if(nm%in%colnames(ra_c2))vv<-as.numeric(ra_c2[[nm]])
  else if(grepl("^Ratio_Th17.1",nm))vv<-as.numeric(ra_c2[["BALF_Th17.1_CCR6pos_CXCR3pos"]])/(as.numeric(ra_c2[["PB_Th17.1_CCR6pos_CXCR3pos"]])+.01)
  else if(grepl("^GSVA_",nm)&exists("gsva_scores")){gs<-gsub("^GSVA_","",nm);vv<-gsva_scores[gs,ra_c2$Sample_ID]}
  else if(grepl("^Deconv_",nm)&exists("deconv_mat")){ct<-gsub("^Deconv_","",nm);vv<-deconv_mat[ra_c2$Sample_ID,ct]}
  if(!is.null(vv)){v<-!is.na(vv)&!is.na(y_inf);if(sum(v)>=10){ro<-roc(y_inf[v],vv[v],quiet=T,direction="auto")
    short<-gsub("BALF_|PB_|GSVA_|Deconv_|Ratio_|_CCR6.*|_CD127.*|_BAL_PB","",nm)
    d<-data.frame(S=ro$sensitivities,Sp=1-ro$specificities)
    d<-d[order(d$Sp,d$S),]
    d$Label<-sprintf("%s(%.2f)",short,as.numeric(auc(ro)))
    d$Rank<-i
    rl<-rbind(rl,d)}}}
rl$Label<-factor(rl$Label,levels=unique(rl$Label))
f4a<-ggplot(rl,aes(Sp,S,color=Label))+geom_step(linewidth=.6,direction="vh")+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey60")+
  scale_color_manual(values=clrs)+labs(x="1-Specificity",y="Sensitivity",title="a  Biomarker ROC",color="")+
  coord_equal()+tn(7)+theme(legend.position=c(.62,.25),legend.text=element_text(size=5.5),legend.key.size=unit(2,"mm"),
                            legend.background=element_rect(fill=alpha("white",.95),color=NA))

ti<-inf_comp%>%filter(P_value<.05)%>%mutate(L=gsub("BALF_|PB_|GSVA_|Deconv_|_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14.*","",Variable))%>%arrange(Cliff_delta)
ti$L<-fct_inorder(ti$L)
f4b<-ggplot(ti,aes(Cliff_delta,L,fill=Category))+geom_col(width=.5)+scale_fill_brewer(palette="Set2")+geom_vline(xintercept=0,linewidth=.2)+
  labs(x="Cliff's delta",y="",title="b",fill="")+tn(6)+theme(legend.position="bottom",legend.key.size=unit(2,"mm"),axis.text.y=element_text(size=5),legend.text=element_text(size=5))

sf<-c("SFTPB","ETV5","LPCAT1","PGC","SFTA3","CLDN18","SFTPD","SFTPA1","SFTPA2")
sfd<-deg_inf[deg_inf$Gene_Symbol%in%sf,]%>%arrange(log2FoldChange);sfd$Gene_Symbol<-fct_inorder(sfd$Gene_Symbol)
sfd$sg<-ifelse(!is.na(sfd$padj)&sfd$padj<.05,"padj<0.05",ifelse(!is.na(sfd$padj)&sfd$padj<.1,"padj<0.1","NS"))
f4c<-ggplot(sfd,aes(log2FoldChange,Gene_Symbol,fill=sg))+geom_col(width=.5)+
  scale_fill_manual(values=c("padj<0.05"="#C0392B","padj<0.1"="#E67E22","NS"="#95A5A6"))+geom_vline(xintercept=0,linewidth=.2)+
  labs(x="log2FC",y="",title="c  Surfactant",fill="")+tn(7)+theme(legend.position="bottom",legend.direction="horizontal",
       legend.text=element_text(size=5),legend.key.size=unit(2,"mm"),plot.margin=unit(c(3,6,3,4),"mm"))

f4d<-bxi("BALF_Th17.1_CCR6pos_CXCR3pos","BALF Th17.1 (%)","d",ra_md)
f4e<-bxi("BALF_Treg_CD127neg_CD25pos","BALF Treg (%)","e",ra_md)

ra_md$Ratio_T171<-as.numeric(ra_md[["BALF_Th17.1_CCR6pos_CXCR3pos"]])/(as.numeric(ra_md[["PB_Th17.1_CCR6pos_CXCR3pos"]])+.01)
dfr<-data.frame(V=ra_md$Ratio_T171,IG=ra_md$IG);dfr<-dfr[!is.na(dfr$V)&!is.na(dfr$IG),]
wtr2<-wilcox.test(dfr$V[dfr$IG=="Infection (+)"],dfr$V[dfr$IG=="Infection (-)"], exact=TRUE)
f4f<-ggplot(dfr,aes(IG,V,fill=IG))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+
  scale_fill_manual(values=CI)+scale_y_log10()+labs(x="",y="Th17.1 BAL/PB\n(log10)",title="f")+
  annotate("text",x=1.5,y=max(dfr$V,na.rm=T)*2,label=sprintf("p=%.3f",wtr2$p.value),size=2)+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=5.5))

ra_ids_inf<-intersect(ra_md$Sample_ID,colnames(gsva_scores))
dfo<-data.frame(V=gsva_scores["OXPHOS",ra_ids_inf],IG=ra_md$IG[match(ra_ids_inf,ra_md$Sample_ID)])
dfo<-dfo[!is.na(dfo$V)&!is.na(dfo$IG),]
wto<-wilcox.test(dfo$V[dfo$IG=="Infection (+)"],dfo$V[dfo$IG=="Infection (-)"], exact=TRUE)
f4g<-ggplot(dfo,aes(IG,V,fill=IG))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+
  scale_fill_manual(values=CI)+labs(x="",y="OXPHOS score",title="g")+annotate("text",x=1.5,y=max(dfo$V)*1.08,label=sprintf("p=%.3f",wto$p.value),size=2)+
  tn(7)+theme(legend.position="none",axis.text.x=element_text(size=5.5),plot.margin=unit(c(5,5,3,4),"mm"))

ginf<-tryCatch(read.csv("results/tables/GSEA_Infection_GO_BP.csv"),error=function(e)NULL)
tgi<-head(ginf[order(ginf$p.adjust),],8)%>%mutate(D=str_trunc(Description,42));tgi<-tgi[order(tgi$NES),];tgi$D<-fct_inorder(tgi$D)
f4h<-ggplot(tgi,aes(NES,D,fill=NES>0))+geom_col(width=.55,show.legend=F)+scale_fill_manual(values=c("TRUE"="#C0392B","FALSE"="#2980B9"))+
  geom_vline(xintercept=0,linewidth=.2)+labs(x="NES",y="",title="h  GSEA Infection")+tn(6)+theme(axis.text.y=element_text(size=5.5))

mi<-tryCatch(read.csv("results/tables/Microbiome_Genus_Infection.csv"),error=function(e){
  tryCatch({md<-read.csv("results/tables/Microbiome_Infection_DiffAbundance.csv")
    md$Genus<-sapply(strsplit(as.character(md$Taxon),";"),function(x) trimws(tail(x,2)[1]))
    md%>%group_by(Genus)%>%dplyr::slice(1)%>%ungroup()%>%dplyr::select(Genus,Log2FC,P_value,P_adjusted)%>%as.data.frame()},error=function(e)NULL)})
if(!is.null(mi)){tmi<-head(mi%>%arrange(P_value),8)%>%arrange(Log2FC);tmi$Genus<-fct_inorder(tmi$Genus);tmi$sg<-ifelse(tmi$P_value<.05,"*","NS")
  f4i<-ggplot(tmi,aes(Log2FC,Genus,fill=sg))+geom_col(width=.5)+scale_fill_manual(values=c("*"="#C0392B","NS"="#95A5A6"))+
    geom_vline(xintercept=0,linewidth=.2)+labs(x="log2FC(Infection+ vs −)",y="",title="i  Microbiome",fill="")+
    tn(6.5)+theme(axis.text.y=element_text(size=5.5,face="italic"),legend.position="bottom")
} else {f4i<-ggplot()+geom_blank()+labs(title="i")+tn()}

ggsave(file.path(fig_dir,"Fig4_Infection.pdf"),
       arrangeGrob(f4a,f4b,f4c,f4d,f4e,f4f,f4g,f4h,f4i,
                   layout_matrix=rbind(c(1,1,2,3),c(4,5,6,7),c(8,8,9,9))),
       width=220,height=230,units="mm",dpi=300)
cat("  done\n")

# ============================================================================
# Fig 5: CT Progression (8 panels)
# ============================================================================
cat("Fig 5...\n")
ra_ext2<-master_ext[master_ext$Sample_Group=="RA" & master_ext$Sample_ID %in% common_samples,]  # n=24
ctd<-data.frame(HU=as.numeric(master_ext$CT_T1_lung_mean_HU),SG=grp_all)%>%filter(!is.na(HU),!is.na(SG))
ctd$SG<-factor(ctd$SG,levels=c("Sarcoidosis","RA-nonILD","RA-ILD"))
ctd$SG_short<-recode(ctd$SG,"Sarcoidosis"="Sarcoidosis","RA-nonILD"="RA-\nnon-ILD","RA-ILD"="RA-\nILD")
ctd$SG_short<-factor(ctd$SG_short,levels=c("Sarcoidosis","RA-\nnon-ILD","RA-\nILD"))
f5a<-ggplot(ctd,aes(SG_short,HU,fill=SG))+geom_boxplot(outlier.shape=NA,width=.45,linewidth=.2)+geom_jitter(width=.07,size=1.2,alpha=.5,show.legend=F)+
  scale_fill_manual(values=CS)+labs(x="",y="Lung mean HU",title="a")+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=6))

f5b<-sct(as.numeric(ra_ext2$ILD_Score),as.numeric(ra_ext2$CT_Delta_GGO_pct),"ILD score","Delta GGO (%)","b")
# RA only for all CT progression scatter plots
f5c<-sct(as.numeric(ra_ext2$`PB_Th17.1_CCR6pos_CXCR3pos`),as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),"PB Th17.1 (%)","Delta Healthy Lung (%)","c")
f5d<-sct(as.numeric(ra_ext2$`PB_Th17_CCR6pos_CXCR3neg`),as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),"PB Th17 (%)","Delta Healthy Lung (%)","d")
f5e<-sct(as.numeric(ra_ext2$BALF_Neutrophil_Percent),as.numeric(ra_ext2$CT_Delta_Healthy_Lung_pct),"BALF Neutrophil (%)","Delta Healthy Lung (%)","e")

# BAL/PB Th17.1 ratio vs infection (dual role)
f5g<-bxi("BALF_Th17.1_CCR6pos_CXCR3pos","BALF Th17.1 (%)","g  BAL Th17.1",ra_md)

ra_md$Ratio_T171_fig<-as.numeric(ra_md[["BALF_Th17.1_CCR6pos_CXCR3pos"]])/(as.numeric(ra_md[["PB_Th17.1_CCR6pos_CXCR3pos"]])+.01)
dfr_fig<-data.frame(V=ra_md$Ratio_T171_fig,IG=ra_md$IG);dfr_fig<-dfr_fig[!is.na(dfr_fig$V)&!is.na(dfr_fig$IG),]
wtr_fig<-wilcox.test(dfr_fig$V[dfr_fig$IG=="Infection (+)"],dfr_fig$V[dfr_fig$IG=="Infection (-)"], exact=TRUE)
f5h<-ggplot(dfr_fig,aes(IG,V,fill=IG))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+
  scale_fill_manual(values=CI)+scale_y_log10()+labs(x="",y="BAL/PB Th17.1 ratio",title="h  BAL/PB ratio")+
  annotate("text",x=1.5,y=max(dfr_fig$V,na.rm=T)*2,label=sprintf("p=%.3f",wtr_fig$p.value),size=2)+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=5.5))

ggsave(file.path(fig_dir,"Fig5_CT_Th171.pdf"),arrangeGrob(f5a,f5b,f5c,f5d,f5e,f5g,f5h,layout_matrix=rbind(c(1,2,3,4),c(5,6,7,NA))),
       width=220,height=170,units="mm",dpi=300)
cat("  done\n")

# ============================================================================
# Fig 6: Multi-omics Integration (MOFA2 + Ablation) — 8 panels
# ============================================================================
cat("Fig 6...\n")

# Load MOFA2 6-view results
if(file.exists("results/MOFA2_6views_Results.RData")) {
  e6 <- new.env(); load("results/MOFA2_6views_Results.RData", envir=e6)

  # a: R2 heatmap
  r2m<-e6$r2$r2_per_factor[[1]];r2l<-as.data.frame(r2m)%>%rownames_to_column("Factor")%>%pivot_longer(-Factor,names_to="View",values_to="R2")
  r2l$Factor<-factor(r2l$Factor,levels=rownames(r2m))
  f6a<-ggplot(r2l,aes(View,Factor,fill=R2))+geom_tile(color="white")+scale_fill_gradient(low="white",high="#C0392B",name="R2(%)")+
    geom_text(aes(label=sprintf("%.1f",R2)),size=1.8)+labs(x="",y="",title="a  MOFA2 (6 independent views)")+
    tn(7)+theme(axis.text.x=element_text(angle=35,hjust=1,size=6))

  # b: Factor scatter
  group_m<-ifelse(grepl("^(KYC|Sarcoidosis)",rownames(e6$factors)),"Sarcoidosis","RA")
  df_m<-data.frame(F1=e6$factors[,1],F2=e6$factors[,2],Group=group_m)
  f6b<-ggplot(df_m,aes(F1,F2,color=Group))+geom_point(size=2.5,alpha=.8)+stat_ellipse(level=.95,linewidth=.5)+
    scale_color_manual(values=CG)+labs(x="Factor 1",y="Factor 2",title="b",color="")+tn()+theme(legend.position="bottom")

  # c: Disease factor boxplot
  sf6<-e6$fd$Factor[e6$fd$p_disease<.1];if(length(sf6)==0)sf6<-1:2
  fb6<-data.frame();for(f in sf6)fb6<-rbind(fb6,data.frame(Factor=sprintf("F%d",f),Value=e6$factors[,f],Group=group_m))
  fb6$Group<-recode(fb6$Group,"Sarcoidosis"="Sarcoidosis")
  f6c<-ggplot(fb6,aes(Group,Value,fill=Group))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+
    geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+scale_fill_manual(values=c("RA"="#C0392B","Sarcoidosis"="#2980B9"))+facet_wrap(~Factor,scales="free_y")+
    labs(x="",y="Factor value",title="c  Disease")+tn(7)+theme(legend.position="none")

  # d: Infection factor boxplot
  inf_m<-master_data$respiratory_infection[match(rownames(e6$factors),master_data$Sample_ID)]
  sfi6<-e6$fi$Factor[e6$fi$p_infection<.1];if(length(sfi6)==0)sfi6<-1:2
  ib6<-data.frame();for(f in sfi6){m<-group_m=="RA"&!is.na(inf_m)
    ib6<-rbind(ib6,data.frame(Factor=sprintf("F%d",f),Value=e6$factors[m,f],IG=ifelse(inf_m[m]==1,"Infection (+)","Infection (-)")))}
  f6d<-ggplot(ib6,aes(IG,Value,fill=IG))+geom_boxplot(outlier.shape=NA,width=.4,linewidth=.2)+
    geom_jitter(width=.07,size=1,alpha=.5,show.legend=F)+scale_fill_manual(values=CI)+facet_wrap(~Factor,scales="free_y")+
    labs(x="",y="Factor value",title="d  Infection (p=0.077)")+tn(7)+theme(legend.position="none",axis.text.x=element_text(size=5))
}

# e: Layer Ablation AUC
if(file.exists("results/tables/LayerAblation_7layers_AUC.csv")) {
  abl<-read.csv("results/tables/LayerAblation_7layers_AUC.csv")
  abl$Model<-fct_reorder(abl$Model,abl$AUC)
  abl$Type<-recode(abl$Type,"Single"="1-layer","Combined"="Full","Ablation"="Leave-out")
  abl$Type<-factor(abl$Type,levels=c("1-layer","Full","Leave-out"))
  f6e<-ggplot(abl,aes(AUC,Model,fill=Type))+geom_col(width=.55)+
    scale_fill_manual(values=c("1-layer"="#3498DB","Full"="#C0392B","Leave-out"="#E67E22"))+
    geom_vline(xintercept=.5,linetype="dashed",color="grey60",linewidth=.2)+
    labs(x="LOOCV AUC",y="",title="e  Layer ablation",fill="")+xlim(0,1)+
    tn(7)+theme(legend.position="bottom",axis.text.y=element_text(size=5.5),legend.text=element_text(size=6))
}

# f: Cross-layer network
if(file.exists("results/Integration_Results.RData")) {
  load("results/Integration_Results.RData")
  edge_sum<-edges%>%mutate(Pair=paste(pmin(Layer1,Layer2),"-",pmax(Layer1,Layer2)))%>%
    group_by(Pair)%>%summarise(N=n(),.groups="drop")%>%arrange(desc(N))
  edge_sum$Pair<-fct_reorder(edge_sum$Pair,edge_sum$N)
  f6f<-ggplot(edge_sum,aes(N,Pair))+geom_col(fill="#34495E",width=.55)+
    labs(x="Edges (|rho|>0.5)",y="",title="f  Network")+tn(7)+theme(axis.text.y=element_text(size=5.5))
}

# g: Integrated heatmap (simplified)
gm_all<-gsva_scores[,common_samples]
so_all<-common_samples[order(meta$Sample_Group[match(common_samples,meta$Sample_ID)])]
glg_all<-as.data.frame(gm_all)%>%rownames_to_column("PW")%>%pivot_longer(-PW,names_to="S",values_to="Sc")
glg_all$S_d<-ifelse(grepl("^(KYC|Sarcoidosis)",glg_all$S),paste0("Sarcoidosis",gsub("^(KYC0*|Sarcoidosis)","",glg_all$S)),paste0("RA",gsub("^(KY0*|RA)","",glg_all$S)))
so_d<-ifelse(grepl("^(KYC|Sarcoidosis)",so_all),paste0("Sarcoidosis",gsub("^(KYC0*|Sarcoidosis)","",so_all)),paste0("RA",gsub("^(KY0*|RA)","",so_all)))
glg_all$S_d<-factor(glg_all$S_d,levels=so_d)
f6g<-ggplot(glg_all,aes(S_d,PW,fill=Sc))+geom_tile()+scale_fill_gradient2(low="#2980B9",mid="white",high="#C0392B",midpoint=0,name="Score")+
  labs(x="",y="",title="g  GSVA overview")+tn(5)+theme(axis.text.x=element_text(angle=55,hjust=1,size=3),axis.text.y=element_text(size=5))

# h: AUC summary comparison
auc_comp<-data.frame(
  Model=c("RA vs Sarcoidosis\n(multi-omics)","BALF Th17.1\n(infection)","BAL/PB Th17.1\nratio","GSVA OXPHOS","BALF MCP3","PB Th17.1\nvs CT Delta"),
  AUC_or_rho=c(0.962,0.870,0.861,0.861,0.866,0.484),
  Type=c("AUC","AUC","AUC","AUC","AUC","rho"))
auc_comp$Model<-fct_inorder(auc_comp$Model)
f6h<-ggplot(auc_comp%>%filter(Type=="AUC"),aes(AUC_or_rho,Model))+geom_col(fill="#C0392B",width=.5)+
  geom_vline(xintercept=.5,linetype="dashed",color="grey60",linewidth=.2)+
  labs(x="AUC",y="",title="h  Performance")+xlim(0,1)+tn(7)+theme(axis.text.y=element_text(size=6))

# i: MOFA2 Factor 1 feature weights
if(file.exists("results/MOFA2_6views_Results.RData")) {
  e7 <- new.env(); load("results/MOFA2_6views_Results.RData", envir=e7)
  fw1 <- data.frame()
  for(v in names(e7$weights)) {
    wv <- e7$weights[[v]][, 1]
    df_v <- data.frame(Feature=names(wv), Weight=as.numeric(wv), View=v, stringsAsFactors=FALSE)
    fw1 <- rbind(fw1, df_v)
  }
  fw1$AbsW <- abs(fw1$Weight)
  fw1_top <- fw1 %>% arrange(desc(AbsW)) %>% head(15)
  fw1_top$Feature <- gsub("^BALF_|^Serum_|^PB_", "", fw1_top$Feature)
  fw1_top$Feature <- gsub("_CCR6.*|_CD127.*|_CD86.*|_CD66.*|_CD14pos.*", "", fw1_top$Feature)
  fw1_top$Feature <- fct_reorder(fw1_top$Feature, fw1_top$AbsW)
  fw1_top$View <- recode(fw1_top$View, "Serum_Cytokine"="Serum", "BALF_Cytokine"="BALF Cyto",
                          "BALF_FACS"="BALF FACS", "PB_FACS"="PB FACS",
                          "Expression"="Expr", "Microbiome"="Micro")
  view_cols <- c("Serum"="#E74C3C", "BALF Cyto"="#3498DB", "BALF FACS"="#27AE60",
                 "PB FACS"="#9B59B6", "Expr"="#E67E22", "Micro"="#34495E")
  f6i <- ggplot(fw1_top, aes(Weight, Feature, fill=View)) +
    geom_col(width=0.6) + scale_fill_manual(values=view_cols) +
    geom_vline(xintercept=0, linewidth=0.2) +
    labs(x="Weight", y="", title="i  Factor 1 weights", fill="") +
    tn(7) + theme(legend.position="bottom", legend.key.size=unit(2.5,"mm"),
                  axis.text.y=element_text(size=6), legend.text=element_text(size=5.5))
} else {
  f6i <- ggplot() + geom_blank() + labs(title="i") + tn()
}

ggsave(file.path(fig_dir,"Fig6_Integration.pdf"),
       arrangeGrob(f6a,f6b,f6c,f6d,f6e,f6f,f6g,f6h,f6i,
                   layout_matrix=rbind(c(1,1,1,2,2,2),c(3,3,4,4,5,5),c(6,6,7,7,8,9))),
       width=220,height=240,units="mm",dpi=300)
cat("  done\n")

# Figure 7 removed — all content integrated into Figures 1-6
# Key findings → discussed in text; MOFA2 weights → Fig 6i; Th17.1 dual role → Fig 5g-h

cat("\n=== All 6 figures complete ===\n")
