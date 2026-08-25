options(stringsAsFactors=FALSE)
if(.Platform$OS.type=="windows") windowsFonts(Arial=windowsFont("Arial"))
root<-Sys.getenv("HCC_PROJECT_ROOT",unset="");if(!nzchar(root))root<-normalizePath(".",winslash="/",mustWork=TRUE)
outdir<-Sys.getenv("HCC_OUTPUT_DIR",unset="");if(!nzchar(outdir))outdir<-file.path(root,"07_figures");dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

registry<-data.frame(
 dataset=c("PRJNA784025 cohort 1","PRJNA784025 cohort 2","PRJEB54571","MTBLS8764","ST001152","PRJEB54571 RNA-seq","TCGA-LIHC","GSE63898","GSE149614"),
 layer=c(rep("Microbiome",3),rep("Metabolome",2),rep("Bulk host",3),"Single cell"),
 design=c("HCC vs cirrhosis","HCC vs cirrhosis","HCC vs cirrhosis","Longitudinal pre-HCC vs controls","Paired tumour-adjacent liver","HCC vs cirrhosis","Paired tumour-adjacent","HCC vs cirrhosis","Paired tumour-non-tumour"),
 independent_unit=c("89 participants","65 participants","105 participants","203 participants; 612 samples","40 pairs; 10 cirrhosis-background pairs","22 participants (17 HCC; 5 cirrhosis)","50 patient pairs","396 participants","10 donors; 8 paired"),
  role=c("discovery","locked validation","external sensitivity","longitudinal candidate discovery","cross-specimen tissue audit","bulk discovery","paired validation","direct external validation","cell-type localisation"),
 status="included"
)
excluded<-data.frame(
 dataset=c("PRJNA647523","PRJNA540574","ST004733"),
 layer=c("Microbiome","Microbiome","Metabolome"),
 disposition=c("summary evidence only","excluded from patient-level pooling","formally unavailable/excluded"),
 reason=c("published summary lacks harmonisable participant-level matrix","disease labels unresolved after metadata audit","required sample-level files unavailable at frozen search date")
)
gates<-data.frame(
 gate=c("Metadata","Independent unit","Harmonisation","Multiplicity","Replication","Claim"),
 rule=c("disease, specimen and pairing verified","participant/donor, never library or cell","common taxon/identifier and contrast only","family-wide BH or frozen gate","independent cohort/model/donor required","no causality without linked mediation/intervention")
)
write.csv(registry,file.path(outdir,"Figure1_source_data_cohort_registry.csv"),row.names=FALSE)
write.csv(excluded,file.path(outdir,"Figure1_source_data_exclusion_registry.csv"),row.names=FALSE)
write.csv(gates,file.path(outdir,"Figure1_source_data_analysis_gates.csv"),row.names=FALSE)

layer_cols<-c(Microbiome="#4477AA",Metabolome="#D08B35","Bulk.host"="#3A7D6B","Single.cell"="#8B6FA8")
lc<-function(x)layer_cols[gsub(" ",".",x)]
draw<-function(){
 layout(matrix(c(1,2,3,4),2,2,byrow=TRUE),heights=c(1.18,.98),widths=c(1.12,.88));par(family="Arial",fg="#222222")
 # a included registry
 par(mar=c(1,1,2.6,1));plot.new();plot.window(c(0,1),c(0,1));title("Nine analyzable cohorts span four evidence layers",adj=.03,cex.main=.91,font.main=2,line=.10)
 yy<-seq(.90,.10,length.out=nrow(registry));for(i in seq_len(nrow(registry))){rect(.02,yy[i]-.038,.98,yy[i]+.038,col=if(i%%2)"#F8F8F8" else "white",border="#E1E1E1");rect(.02,yy[i]-.038,.035,yy[i]+.038,col=lc(registry$layer[i]),border=NA);text(.05,yy[i],registry$dataset[i],adj=0,font=2,cex=.48);text(.30,yy[i],registry$design[i],adj=0,cex=.43);text(.63,yy[i],registry$independent_unit[i],adj=0,cex=.42);text(.96,yy[i],registry$role[i],adj=1,cex=.40,col="#555555")};mtext("A",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
 # b layer map
 par(mar=c(1,1,2.6,1));plot.new();plot.window(c(0,1),c(0,1));title("Evidence layers are integrated without pooling unlike units",adj=.03,cex.main=.90,font.main=2,line=.10)
 bx<-c(.16,.50,.84);by<-c(.74,.74,.74);lab<-c("Microbiome\n3 cohorts | 259 participants","Metabolome\nlongitudinal + paired tissue","Bulk host\n3 independent datasets")
 for(i in 1:3){rect(bx[i]-.14,by[i]-.11,bx[i]+.14,by[i]+.11,col=paste0(lc(c("Microbiome","Metabolome","Bulk host"))[i],"22"),border=lc(c("Microbiome","Metabolome","Bulk host"))[i],lwd=1.2);text(bx[i],by[i],lab[i],font=2,cex=.57,col="#365A6C")}
 arrows(.30,.74,.35,.74,length=.07,col="#777777",lty=2);arrows(.64,.74,.69,.74,length=.07,col="#777777",lty=2)
 rect(.34,.35,.66,.53,col="#F3EDF7",border=lc("Single cell"),lwd=1.2);text(.50,.44,"Single cell\n71,915 cells | 10 donors | 8 paired",font=2,cex=.60,col="#5D4770")
 arrows(.84,.63,.63,.53,length=.07,col="#777777",lty=2)
 text(.50,.18,"Dashed connections denote evidence triangulation, not participant-linked mediation",cex=.52,col="#8A3333");mtext("B",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
 # c exclusions
 par(mar=c(1,1,2.6,1));plot.new();plot.window(c(0,1),c(0,1));title("Three datasets are downgraded or excluded with traceable reasons",adj=.03,cex.main=.90,font.main=2,line=.10)
 yy<-c(.76,.50,.24);for(i in 1:3){rect(.03,yy[i]-.085,.97,yy[i]+.085,col="#FFF5F0",border="#C95A5A");text(.06,yy[i]+.028,excluded$dataset[i],adj=0,font=2,cex=.58);text(.94,yy[i]+.028,excluded$disposition[i],adj=1,font=2,cex=.48,col="#8A3333");text(.06,yy[i]-.034,excluded$reason[i],adj=0,cex=.46,col="#555555")};mtext("C",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
 # d gates
 par(mar=c(1,1,2.6,1));plot.new();plot.window(c(0,1),c(0,1));title("Every promoted claim passes six workflow safeguards",adj=.03,cex.main=.90,font.main=2,line=.10)
 yy<-seq(.84,.16,length.out=6);for(i in 1:6){rect(.04,yy[i]-.045,.25,yy[i]+.045,col="#EAF0F2",border="#6F8794");text(.145,yy[i],paste0(i,"  ",gates$gate[i]),font=2,cex=.50,col="#365A6C");text(.29,yy[i],gates$rule[i],adj=0,cex=.43);if(i<6)arrows(.145,yy[i]-.05,.145,yy[i+1]+.05,length=.05,col="#777777")};mtext("D",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
}
base<-file.path(outdir,"Figure1_cohort_registry_harmonisation");wmm<-180;hmm<-150;win<-wmm/25.4;hin<-hmm/25.4
if(requireNamespace("svglite",quietly=TRUE))svglite::svglite(paste0(base,".svg"),width=win,height=hin,system_fonts=list(Arial="Arial"),pointsize=8)else svg(paste0(base,".svg"),width=win,height=hin,family="Arial",pointsize=8);par(oma=c(.2,.2,.4,.2));draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=win,height=hin,family="Arial",pointsize=8);par(oma=c(.2,.2,.4,.2));draw();dev.off()
png(paste0(base,".png"),width=wmm,height=hmm,units="mm",res=600,pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
tiff(paste0(base,".tiff"),width=wmm,height=hmm,units="mm",res=600,compression="lzw",pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
