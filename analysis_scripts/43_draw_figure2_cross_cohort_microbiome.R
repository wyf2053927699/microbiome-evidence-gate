options(stringsAsFactors=FALSE)
if(.Platform$OS.type=="windows")windowsFonts(Arial=windowsFont("Arial"))
root<-Sys.getenv("HCC_PROJECT_ROOT",unset="");if(!nzchar(root))root<-normalizePath(".",winslash="/",mustWork=TRUE)
outdir<-Sys.getenv("HCC_OUTPUT_DIR",unset="");if(!nzchar(outdir))outdir<-file.path(root,"07_figures");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
x<-read.csv(file.path(root,"06_results","microbiome_cross_cohort","three_cohort_harmonised_genus_effects.csv"),check.names=FALSE)
loo<-read.csv(file.path(root,"06_results","microbiome_cross_cohort","three_cohort_leave_one_cohort_out.csv"),check.names=FALSE)
sp<-split(loo,loo$genus);stable_names<-names(Filter(function(z)all(z$random_p<.05)&&length(unique(sign(z$random_estimate)))==1,sp))
stable<-x[x$genus%in%stable_names,];stable<-stable[order(stable$random_estimate),]
loo_stable<-loo[loo$genus%in%stable_names,]
flow<-data.frame(stage=c("Harmonised genera","Discovery q<0.10","Locked validation pass","External age-adjusted pass"),n=c(70,44,30,0))
write.csv(stable,file.path(outdir,"Figure2_source_data_stable_genera.csv"),row.names=FALSE)
write.csv(loo_stable,file.path(outdir,"Figure2_source_data_stable_genera_leave_one_out.csv"),row.names=FALSE)
write.csv(flow,file.path(outdir,"Figure2_source_data_promotion_flow.csv"),row.names=FALSE)

cols<-c(discovery="#4477AA",validation="#D08B35",external="#3A7D6B",pooled="#222222")
short<-function(z)gsub(" group$"," grp",sub("Lachnospiraceae ","Lachno. ",sub("[Eubacterium] ","Eub. ",z,fixed=TRUE),fixed=TRUE))
draw<-function(){
 layout(matrix(c(1,1,2,3,4,5),3,2,byrow=TRUE),heights=c(.5,1.18,1),widths=c(1.08,.92));par(family="Arial",fg="#222222")
 # a design
par(mar=c(.8,.8,2,.8));plot.new();plot.window(c(0,1),c(0,1));title("Discovery, locked validation and external sensitivity are kept separate",adj=.03,cex.main=.91,font.main=2,line=.55)
bx<-c(.18,.50,.82);labs<-c("PRJNA784025 cohort 1\n30 HCC | 59 cirrhosis\nDISCOVERY","PRJNA784025 cohort 2\n25 HCC | 40 cirrhosis\nLOCKED VALIDATION","PRJEB54571\n95 HCC | 10 cirrhosis\nEXTERNAL SENSITIVITY")
 for(i in 1:3){rect(bx[i]-.13,.25,bx[i]+.13,.80,col=paste0(cols[i],"18"),border=cols[i],lwd=1.2);text(bx[i],.54,labs[i],font=2,cex=.60,col="#365A6C")};arrows(.31,.53,.37,.53,length=.07,col="#777777");arrows(.63,.53,.69,.53,length=.07,col="#777777");text(.82,.12,"Comparator imbalance requires qualified interpretation",cex=.47,col="#8A3333");mtext("A",3,line=.7,at=0,adj=0,font=2,cex=1.05)
 # b conventional versus modified Hartung-Knapp uncertainty for the DL-stable subset
 par(mar=c(3.8,7.1,2,1));yy<-seq_len(nrow(stable));lim<-range(c(stable$hksj_ci_low_modified,stable$hksj_ci_high_modified));plot(stable$random_estimate,yy,xlim=lim+ c(-.2,.2),ylim=c(.5,nrow(stable)+.5),yaxt="n",pch=19,col=cols['pooled'],xlab="Random-effects CLR difference (HCC - cirrhosis)",ylab="",cex=.7);segments(stable$hksj_ci_low_modified,yy,stable$hksj_ci_high_modified,yy,lwd=2,col="#888888");segments(stable$random_ci_low,yy,stable$random_ci_high,yy,lwd=4,col="#3A7D6B88");points(stable$random_estimate,yy,pch=19,col=cols['pooled'],cex=.65);abline(v=0,lty=2,col="#999999");axis(2,at=yy,labels=short(stable$genus),las=1,tick=FALSE,cex.axis=.56);title("Small-k intervals weaken the conventional DL-stable subset",adj=.02,cex.main=.79,font.main=2,line=.55);legend("bottomright",c("Modified Hartung-Knapp 95% CI","Conventional DL 95% CI"),lwd=c(2,4),col=c("#888888","#3A7D6B88"),bty="n",cex=.44);mtext("B",3,line=1.45,at=par("usr")[1],adj=0,font=2,cex=1.05)
 # c all genus map
par(mar=c(4,4.2,2,1));cc<-ifelse(x$three_cohort_direction_concordant,"#3A7D6B88","#BBBBBB88");plot(x$random_estimate,x$i2_percent,pch=19,col=cc,xlab="Random-effects CLR difference",ylab=expression(paste("Heterogeneity ",I^2," (%)")),cex=.65);abline(v=0,lty=2,col="#999999");abline(h=75,lty=3,col="#C95A5A");title("Direction concordance does not remove heterogeneity",adj=.05,cex.main=.86,font.main=2,line=.55);text(par("usr")[2],78,expression(I^2==75*"%"),adj=1,cex=.48,col="#8A3333");legend("topleft",c("3-cohort direction concordant","other"),pch=19,col=c("#3A7D6B","#BBBBBB"),bty="n",cex=.50);mtext("C",3,line=1.45,at=par("usr")[1],adj=0,font=2,cex=1.05)
 # d promotion
 par(mar=c(1.5,1,2,1));plot.new();plot.window(c(0,1),c(0,1));title("A universal genus claim fails the external gate",adj=.03,cex.main=.87,font.main=2,line=.55);xx<-c(.12,.38,.64,.88);ww<-c(.20,.17,.14,.10);for(i in 1:4){rect(xx[i]-ww[i]/2,.47,xx[i]+ww[i]/2,.70,col=if(i<4)"#EAF0F2" else "#FFF0EF",border=if(i<4)"#6F8794" else "#C95A5A",lwd=1.2);text(xx[i],.61,flow$n[i],font=2,cex=.86,col=if(i<4)"#365A6C" else "#8A3333");text(xx[i],.38,flow$stage[i],cex=.48);if(i<4)arrows(xx[i]+ww[i]/2,.585,xx[i+1]-ww[i+1]/2,.585,length=.06,col="#777777")};rect(.18,.10,.82,.22,col="#FFF5F0",border="#C95A5A");text(.50,.16,"25/70 are direction-concordant, but 0/70 pass the full three-stage gate",font=2,cex=.51,col="#8A3333");mtext("D",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
 # e leave-one-out range
par(mar=c(5.4,6.7,2,1));rng<-do.call(rbind,lapply(sp[stable_names],function(z)data.frame(genus=z$genus[1],lo=min(z$random_estimate),hi=max(z$random_estimate),mid=mean(z$random_estimate))));rng<-rng[order(rng$mid),];yy<-seq_len(nrow(rng));plot(rng$mid,yy,xlim=range(c(rng$lo,rng$hi)),ylim=c(.5,nrow(rng)+.5),yaxt="n",pch=19,col="#3A7D6B",xlab="",ylab="",cex=.65);segments(rng$lo,yy,rng$hi,yy,lwd=2,col="#83AFA2");abline(v=0,lty=2,col="#999999");axis(2,at=yy,labels=short(rng$genus),las=1,tick=FALSE,cex.axis=.50);title("Direction persists, but small-k significance does not",adj=.02,cex.main=.82,font.main=2,line=.55);mtext("Leave-one-cohort-out pooled CLR difference",side=1,line=2.2,cex=.58);mtext("0/210 leave-one-out rows have modified Hartung-Knapp P<0.05",side=1,line=3.8,cex=.43,col="#8A3333");mtext("E",3,line=1.45,at=par("usr")[1],adj=0,font=2,cex=1.05)
}
base<-file.path(outdir,"Figure2_cross_cohort_microbiome");wmm<-180;hmm<-150;win<-wmm/25.4;hin<-hmm/25.4
if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(paste0(base,".svg"),width=win,height=hin,system_fonts=list(Arial="Arial"),pointsize=8)else svg(paste0(base,".svg"),width=win,height=hin,family="Arial",pointsize=8);par(oma=c(.2,.2,.4,.2));draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=win,height=hin,family="Arial",pointsize=8);par(oma=c(.2,.2,.4,.2));draw();dev.off()
png(paste0(base,".png"),width=wmm,height=hmm,units="mm",res=600,pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
tiff(paste0(base,".tiff"),width=wmm,height=hmm,units="mm",res=600,compression="lzw",pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
