options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset = "")
if (!nzchar(outdir)) outdir <- file.path(root, "07_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loo <- read.csv(file.path(root,"06_results","microbiome_cross_cohort","three_cohort_leave_one_cohort_out.csv"))
met <- read.csv(file.path(root,"06_results","integration","figure3_mtbls8764_four_model_robustness.csv"))
axis <- read.csv(file.path(root,"06_results","axis_prioritisation","figure4_rank_stability_summary.csv"))
bulk <- read.csv(file.path(root,"06_results","integration","figure5_host_program_contradiction_ledger.csv"))
sc <- read.csv(file.path(root,"06_results","single_cell","paired_programs","prespecified_program_paired_results.csv"))

genera <- split(loo, loo$genus)
stable_genera <- names(Filter(function(z) all(z$random_p < 0.05) && length(unique(sign(z$random_estimate))) == 1, genera))
robust <- data.frame(
  layer=c("Cross-cohort microbiome","Longitudinal metabolome","Mechanism-axis ranking","Independent bulk host","Single-cell host"),
  primary=c("13 DL omission-stable; 0 mHKSJ",paste0(sum(met$longitudinal_promotion)," exploratory plasma rows"),"4 ranked candidates","3 concordant programmes","3 paired-t BH signals"),
  perturbation=c("Leave one cohort out + mHKSJ","4 scaling/model combinations","4 prespecified weight schemes","3 independent liver datasets","Wilcoxon + exact sign-flip + donor omission"),
  outcome=c("Directional only; external gate fails","Both GEE q<0.10; not FDR<0.05","Ranks stable; all essential gates fail","Three concordant; one contradiction","Wilcoxon/exact q=0.086; 0/22 exact-FDR passes"),
  status=c("qualified","qualified","no_go","qualified","qualified")
)
write.csv(robust,file.path(outdir,"Figure7_source_data_robustness_registry.csv"),row.names=FALSE)

negative <- data.frame(
  evidence_item=c("Cross-cohort taxon heterogeneity","Exact cross-specimen metabolite overlap","Mechanism-axis eligibility","Respiratory electron transport","Axis-specific single-cell programmes","Causal/mediation link"),
  observed=c(
    paste0(sum(vapply(genera,function(z) max(z$i2_percent,na.rm=TRUE)>=75,logical(1)))," genera reach I2 >=75% in at least one omission"),
    "No exact metabolite passes FDR in both specimens",
    "0/4 candidates pass all essential gates",
    "Significant in three datasets but directionally discordant",
    "0 proposed axis programmes pass global FDR; three broad paired-t signals have Wilcoxon and exact sign-flip q=0.086",
    "No linked participant-level mediation or intervention"
  ),
  disposition=c("retain heterogeneity qualifier","do not infer biochemical continuity","No-Go; exploratory ranking only","report contradiction, not replication","retain explicit null result","association-only working model")
)
write.csv(negative,file.path(outdir,"Figure7_source_data_negative_evidence_ledger.csv"),row.names=FALSE)

claim <- data.frame(
  claim=c("Qualified microbial directional differences","Longitudinal human metabolite candidates","Broad host-response associations","Specific microbiota-metabolite-host axis","Causal disease mechanism"),
  evidence=c("three cohorts; external and small-k limits retained","four-model directional stability; both GEE q<0.10","three bulk datasets plus donor-aware localisation","essential gates and cell-specific tests fail","no mediation or intervention"),
  allowed=c(TRUE,TRUE,TRUE,FALSE,FALSE),
  wording=c("directional association only","exploratory candidate","qualified association","not supported","not tested")
)
write.csv(claim,file.path(outdir,"Figure7_source_data_claim_boundary.csv"),row.names=FALSE)

cols <- c(pass="#3A7D6B",qualified="#D1A65A",no_go="#C95A5A")
draw <- function(){
  layout(matrix(c(1,1,2,3,4,5),3,2,byrow=TRUE),heights=c(.55,1.15,1.05),widths=c(1.08,.92))
  par(family="Arial",fg="#222222",col.axis="#222222",col.lab="#222222")
  # a registry
  par(mar=c(.8,.8,2,.8));plot.new();plot.window(c(0,1),c(0,1));
  title("Prespecified perturbations are retained across every evidence layer",adj=.03,cex.main=.92,font.main=2,line=.55)
  xs <- seq(.10,.90,length.out=5)
  for(i in 1:5){rect(xs[i]-.085,.19,xs[i]+.085,.84,col="#F5F7F8",border="#738995",lwd=1.1);text(xs[i],.65,robust$primary[i],font=2,cex=.64,col="#365A6C");text(xs[i],.42,robust$perturbation[i],cex=.50);rect(xs[i]-.06,.24,xs[i]+.06,.32,col=cols[robust$status[i]],border=NA);text(xs[i],.28,toupper(gsub("_","-",robust$status[i])),cex=.45,col="white",font=2)}
  mtext("A",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
  # b matrix
  par(mar=c(6.2,7.8,2,1));plot(NA,xlim=c(.5,3.5),ylim=c(.5,5.5),axes=FALSE,xlab="",ylab="")
  title("Sensitivity and decision matrix",adj=.02,cex.main=.9,font.main=2,line=.55)
  labs <- c("Primary result","Perturbation","Decision")
  axis(3,at=1:3,labels=labs,tick=FALSE,line=-1.45,cex.axis=.62)
  axis(2,at=5:1,labels=robust$layer,las=1,tick=FALSE,cex.axis=.62)
  for(i in 1:5){y=6-i;rect(.6,y-.35,3.4,y+.35,col="#FAFAFA",border="#DDDDDD");text(1,y,robust$primary[i],cex=.54);text(2,y,robust$perturbation[i],cex=.52);rect(2.72,y-.18,3.28,y+.18,col=cols[robust$status[i]],border=NA);text(3,y,if(robust$status[i]=="pass")"PASS" else if(robust$status[i]=="qualified")"QUALIFIED" else "NO-GO",col="white",font=2,cex=.49)}
  mtext("B",3,line=1.45,at=par("usr")[1],adj=0,font=2,cex=1.05)
  # c ledger
  par(mar=c(2,1,2,1));plot.new();plot.window(c(0,1),c(0,1));title("Material negative and contradictory evidence",adj=.03,cex.main=.9,font.main=2,line=.55)
  yy<-seq(.87,.12,length.out=nrow(negative));for(i in seq_len(nrow(negative))){rect(.02,yy[i]-.06,.98,yy[i]+.06,col=if(i%%2)"#F7F7F7" else "#FFFFFF",border="#DDDDDD");text(.04,yy[i]+.024,negative$evidence_item[i],adj=0,font=2,cex=.51);text(.04,yy[i]-.026,negative$observed[i],adj=0,cex=.43,col="#444444");text(.96,yy[i]-.026,negative$disposition[i],adj=1,cex=.43,col="#8A3333")};mtext("C",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
  # d bounded model
  par(mar=c(2.5,1,2,1));plot.new();plot.window(c(0,1),c(0,1));title("Bounded working model",adj=.03,cex.main=.9,font.main=2,line=.55)
  bx<-c(.12,.38,.64,.88);bl<-c("Microbiome\nassociations","Human metabolite\ncandidates","Broad host\nprogrammes","HCC context")
  for(i in 1:4){rect(bx[i]-.09,.54,bx[i]+.09,.76,col="#EEF3F4",border="#587889",lwd=1.2);text(bx[i],.65,bl[i],font=2,cex=.61,col="#365A6C")}
  arrows(bx[1]+.09,.65,bx[2]-.09,.65,lty=2,col="#888888",length=.08);arrows(bx[2]+.09,.65,bx[3]-.09,.65,lty=2,col="#888888",length=.08);arrows(bx[3]+.09,.65,bx[4]-.09,.65,lty=2,col="#888888",length=.08)
  text(.5,.39,"Dashed links are hypotheses: cohorts are not participant-linked and no mediation was tested",cex=.55,col="#8A3333")
  rect(.15,.14,.85,.29,col="#FFF0EF",border="#C95A5A",lwd=1.2);text(.5,.215,"Supported: layer-specific human associations   |   Unsupported: a continuous causal mechanism axis",font=2,cex=.56,col="#8A3333");mtext("D",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
  # e claim boundary
  par(mar=c(3.8,1,2,1));plot.new();plot.window(c(0,1),c(0,1));title("Claim boundary",adj=.03,cex.main=.9,font.main=2,line=.55)
  yy<-seq(.84,.18,length.out=nrow(claim));for(i in seq_len(nrow(claim))){cc<-if(claim$allowed[i])"#E7F2EE" else "#FFF0EF";bc<-if(claim$allowed[i])"#3A7D6B" else "#C95A5A";rect(.03,yy[i]-.055,.97,yy[i]+.055,col=cc,border=bc);text(.05,yy[i],claim$claim[i],adj=0,font=2,cex=.53);text(.95,yy[i],toupper(claim$wording[i]),adj=1,cex=.47,col=bc)}
  mtext("Association-level evidence can motivate validation but cannot substitute for linked mediation or intervention",1,line=2.2,cex=.55,col="#666666");mtext("E",3,line=1.45,at=0,adj=0,font=2,cex=1.05)
}

base <- file.path(outdir,"Figure7_robustness_bounded_model")
width_mm<-180; height_mm<-150; width_in<-width_mm/25.4; height_in<-height_mm/25.4
if(requireNamespace("svglite",quietly=TRUE)) svglite::svglite(paste0(base,".svg"),width=width_in,height=height_in,system_fonts=list(Arial="Arial"),pointsize=8) else svg(paste0(base,".svg"),width=width_in,height=height_in,family="Arial",pointsize=8)
par(oma=c(.2,.2,.4,.2));draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=width_in,height=height_in,family="Arial",pointsize=8);par(oma=c(.2,.2,.4,.2));draw();dev.off()
png(paste0(base,".png"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
tiff(paste0(base,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,compression="lzw",pointsize=8,family="Arial",bg="white");par(oma=c(.2,.2,.4,.2));draw();dev.off()
