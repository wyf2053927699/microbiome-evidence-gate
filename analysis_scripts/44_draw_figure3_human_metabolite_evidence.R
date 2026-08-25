options(stringsAsFactors=FALSE)
root<-Sys.getenv("HCC_PROJECT_ROOT",unset="");if(!nzchar(root))root<-normalizePath(".",winslash="/",mustWork=TRUE)
outdir<-Sys.getenv("HCC_OUTPUT_DIR",unset="");if(!nzchar(outdir))outdir<-file.path(root,"07_figures");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,"05_scripts","figure_theme_times_q2.R"))
plasma<-read.csv(file.path(root,"06_results","integration","figure3_mtbls8764_promoted_longitudinal_metabolites.csv"),check.names=FALSE)
tissue<-read.csv(file.path(root,"06_results","ST001152","paired_metabolomics","st001152_promoted_identified_metabolite_results.csv"),check.names=FALSE)
overlap<-read.csv(file.path(root,"06_results","integration","figure3_exact_metabolite_cross_specimen_overlap.csv"),check.names=FALSE)
plasma$label<-paste0(plasma$metabolite," [",ifelse(plasma$window=="within_12m","<=12 m","all pre-HCC"),"]")
plasma$effect_min<-apply(plasma[,grep("^beta_group_",names(plasma))],1,min,na.rm=TRUE);plasma$effect_max<-apply(plasma[,grep("^beta_group_",names(plasma))],1,max,na.rm=TRUE);plasma$effect_mid<-(plasma$effect_min+plasma$effect_max)/2
tissue$short_tier<-ifelse(grepl("level_1",tissue$identification_tier),"Level 1 authentic standard","Putative level 2")
decision<-data.frame(metric=c("Longitudinal plasma exploratory candidates","Paired tissue identified promotions","Exact cross-specimen rows","Exact metabolites passing FDR in both specimens"),n=c(nrow(plasma),nrow(tissue),nrow(overlap),sum(overlap$both_fdr_promoted)),interpretation=c("direction-consistent; both GEE q<0.10","paired all-adjacent contrast","five unique names across seven platform rows","No-Go for biochemical continuity"))
write.csv(plasma,file.path(outdir,"Figure3_source_data_longitudinal_plasma_promotions.csv"),row.names=FALSE)
write.csv(tissue,file.path(outdir,"Figure3_source_data_paired_tissue_promotions.csv"),row.names=FALSE)
write.csv(overlap,file.path(outdir,"Figure3_source_data_exact_cross_specimen_overlap.csv"),row.names=FALSE)
write.csv(decision,file.path(outdir,"Figure3_source_data_promotion_decisions.csv"),row.names=FALSE)

pos=Q2_COL[["teal"]];neg=Q2_COL[["red"]];level1=Q2_COL[["blue"]];level2=Q2_COL[["amber"]]
short<-function(z){z<-gsub("gamma-glutamyltryptophan","gamma-Glu-Trp",z,fixed=TRUE);z<-gsub("1-palmitoleoyl-GPC (16:1)*","palmitoleoyl-GPC",z,fixed=TRUE);gsub("taurocholenate sulfate*","T-cholenate SO4",z,fixed=TRUE)}
draw<-function(){
 par(oma=c(.15,.15,.15,.15));layout(matrix(c(1,1,2,3,4,5),3,2,byrow=TRUE),heights=c(.48,1.08,1.00),widths=c(1.02,.98));par(family=Q2_FONT,fg=Q2_COL[["ink"]])
 # a
 par(mar=c(.45,.45,1.1,.35));plot.new();plot.window(c(0,1),c(0,1));q2_panel("A","Specimen- and design-specific human metabolite evidence",title_adj=.055)
 bx<-c(.20,.52,.82);labs<-c("MTBLS8764 plasma\n203 participants | 612 samples\nlongitudinal pre-HCC","ST001152 liver\n40 tumour-adjacent pairs\n10 cirrhosis-background pairs","Exact-name audit\n5 unique metabolites\n7 platform-specific rows")
for(i in 1:3){q2_card(bx[i]-.14,.25,bx[i]+.14,.80,c(Q2_COL[["blue_light"]],Q2_COL[["teal_light"]],Q2_COL[["amber_light"]])[i],c(level1,pos,level2)[i],1.0);text(bx[i],.53,labs[i],font=2,cex=.57,col=c(level1,pos,level2)[i])};q2_arrow(.34,.53,.38,.53);q2_arrow(.66,.53,.68,.53);text(.50,.11,"Specimen layers are not participant-linked",cex=.55,col=neg,font=3)
 # b plasma
par(mar=c(5.0,8.4,1.1,.55));plasma<-plasma[order(plasma$effect_mid),];yy<-seq_len(nrow(plasma));cc<-ifelse(plasma$effect_mid>0,pos,neg);plot(plasma$effect_mid,yy,xlim=range(c(plasma$effect_min,plasma$effect_max))+c(-.15,.15),ylim=c(.5,nrow(plasma)+.5),yaxt="n",pch=19,col=cc,xlab="",ylab="",cex=.70);segments(plasma$effect_min,yy,plasma$effect_max,yy,lwd=1.8,col=cc);abline(v=0,lty=2,col=Q2_COL[["pale"]]);axis(2,at=yy,labels=short(plasma$label),las=1,tick=FALSE,cex.axis=.56);q2_panel("B","Six stable exploratory plasma candidates",title_adj=.105);mtext("Longitudinal group effect across four models",side=1,line=2.0,cex=.61);mtext("Both GEE models BH q<0.10; not FDR<0.05 confirmation",side=1,line=3.45,cex=.51,col=level2)
 # c tissue
 par(mar=c(3.4,8.6,1.1,.55));tissue<-tissue[order(tissue$median_log2_difference),];yy<-seq_len(nrow(tissue));cc<-ifelse(tissue$short_tier=="Level 1 authentic standard",level1,level2);plot(tissue$median_log2_difference,yy,xlim=range(tissue$median_log2_difference)+c(-.2,.2),ylim=c(.5,nrow(tissue)+.5),yaxt="n",pch=19,col=cc,xlab="Median paired log2 difference (HCC - adjacent)",ylab="",cex=.70);abline(v=0,lty=2,col=Q2_COL[["pale"]]);segments(0,yy,tissue$median_log2_difference,yy,col=cc,lwd=1.1);axis(2,at=yy,labels=tissue$metabolite,las=1,tick=FALSE,cex.axis=.50);legend("bottomright",c("Level 1 authentic standard","Putative level 2"),pch=19,col=c(level1,level2),bty="n",cex=.52);q2_panel("C","Tissue features passing paired FDR",title_adj=.105)
 # d overlap
 par(mar=c(3.4,4.0,1.1,.55));cc<-ifelse(overlap$direction_concordant,pos,neg);xr<-range(overlap$beta_group_provided_area_gee)+c(-.04,.08);yr<-range(overlap$median_log2_difference)+c(-.18,.18);plot(overlap$beta_group_provided_area_gee,overlap$median_log2_difference,pch=19,col=cc,xlim=xr,ylim=yr,xlab="Plasma longitudinal GEE effect",ylab="Paired tissue effect",cex=.70);abline(h=0,v=0,lty=2,col=Q2_COL[["pale"]]);labdat<-aggregate(cbind(beta_group_provided_area_gee,median_log2_difference)~name_key,overlap,mean);labdat$lab<-c(glucose="glucose",glycerol3phosphate="glycerol 3-phosphate",hypoxanthine="hypoxanthine",propionylcarnitine="propionylcarnitine",xanthosine="xanthosine")[labdat$name_key];text(labdat$beta_group_provided_area_gee,labdat$median_log2_difference,labels=labdat$lab,pos=c(4,4,2,2,2),offset=.40,cex=.49,col=Q2_COL[["ink"]]);q2_panel("D","No dual-specimen FDR support",title_adj=.105);legend("topleft",c("direction concordant","discordant"),pch=19,col=c(pos,neg),bty="n",cex=.52)
 # e decision and identification
par(mar=c(.6,.5,1.1,.4));plot.new();plot.window(c(0,1),c(0,1));q2_panel("E","Identification confidence and evidence boundary",title_adj=.105)
 n1<-sum(tissue$short_tier=="Level 1 authentic standard");n2<-nrow(tissue)-n1;rect(.05,.62,.44,.80,col="#EAF0F7",border=level1);text(.245,.71,paste0(n1," tissue features\nLevel 1 authentic standard"),font=2,cex=.60,col="#365A6C");rect(.56,.62,.95,.80,col="#F8EFE3",border=level2);text(.755,.71,paste0(n2," tissue features\nputative level 2"),font=2,cex=.60,col="#7A4D18")
 yy<-c(.48,.36,.24,.12);for(i in 1:4){text(.05,yy[i],decision$metric[i],adj=0,font=if(i==4)2 else 1,cex=.48);text(.68,yy[i],decision$n[i],font=2,cex=.60,col=if(i==4)neg else level1);text(.95,yy[i],decision$interpretation[i],adj=1,cex=.43,col=if(i==4)neg else Q2_COL[["neutral"]])}
}
q2_save_base(draw,file.path(outdir,"Figure3_human_metabolite_evidence"),width_mm=180,height_mm=140)
