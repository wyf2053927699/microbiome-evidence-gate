options(stringsAsFactors = FALSE)
root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
fd <- file.path(root, "07_figures")
rd <- file.path(root, "06_results")
out <- file.path(fd, "core_results_redesign_approval_2026-08-20")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

font <- "Times New Roman"
if (.Platform$OS.type == "windows") windowsFonts(`Times New Roman` = windowsFont("Times New Roman"))
P <- c(ink="#18242B", muted="#62727A", rule="#D5DEE2", pale="#F5F8F9",
       blue="#3578A8", teal="#2F8B7E", amber="#C18124", violet="#735EA3",
       red="#C45151", green="#3D805B")
alpha_col <- function(col, a=.18) adjustcolor(col, alpha.f=a)
short <- function(x,n=28) ifelse(nchar(x)>n,paste0(substr(x,1,n-3),"..."),x)
wrap <- function(x, width=24) vapply(as.character(x), function(s) paste(strwrap(gsub("_", " ", s), width=width), collapse="\n"), character(1))
truth <- function(x) tolower(as.character(x)) %in% c("true","yes","1")

tag <- function(x) {
  xx <- grconvertX(par("fig")[1]+.012,"ndc","user")
  yy <- grconvertY(par("fig")[4]-.012,"ndc","user")
  text(xx,yy,x,adj=c(0,1),font=2,cex=1.05,xpd=NA,family=font,col=P["ink"])
}
base_par <- function(mar) par(family=font,mar=mar,mgp=c(2.1,.65,0),tcl=-.22,
                              col.axis=P["ink"],col.lab=P["ink"],fg=P["ink"])
save_pub <- function(draw, stem, height_mm, pointsize=11.5) {
  b <- file.path(out,stem); w <- 183
  png(paste0(b,".png"),w,height_mm,units="mm",res=300,pointsize=pointsize,family=font,bg="white");draw();dev.off()
  tiff(paste0(b,".tiff"),w,height_mm,units="mm",res=600,compression="lzw",pointsize=pointsize,family=font,bg="white");draw();dev.off()
  cairo_pdf(paste0(b,".pdf"),w/25.4,height_mm/25.4,pointsize=pointsize,family=font);draw();dev.off()
  svg(paste0(b,".svg"),w/25.4,height_mm/25.4,pointsize=pointsize,family=font);draw();dev.off()
}
heat <- function(mat, rows, cols, lim=max(abs(mat),na.rm=TRUE), values=TRUE) {
  nr<-nrow(mat); nc<-ncol(mat); plot.new();plot.window(c(0,nc),c(0,nr))
  for(i in seq_len(nr)) for(j in seq_len(nc)) {
    v<-mat[i,j]; frac<-if(is.na(v))0 else min(1,abs(v)/lim)
    cc<-if(is.na(v))P["pale"] else if(v>=0)alpha_col(P["red"],.18+.68*frac) else alpha_col(P["blue"],.18+.68*frac)
    rect(j-1,nr-i,j,nr-i+1,col=cc,border="white")
    if(values && !is.na(v)) text(j-.5,nr-i+.5,sprintf("%.2f",v),cex=.56,col=if(frac>.65)"white" else P["ink"])
  }
  axis(2,nr-seq_len(nr)+.5,rows,las=1,tick=FALSE,cex.axis=.68)
  axis(3,seq_len(nc)-.5,wrap(cols,18),tick=FALSE,cex.axis=.68)
  box(col=P["rule"])
}

# Figure 1: data-rich study landscape
cr <- read.csv(file.path(fd,"Figure1_source_data_cohort_roles.csv"),check.names=FALSE)
cr$role_group <- ifelse(grepl("[Dd]iscovery",cr$role),"Discovery",
                 ifelse(grepl("[Ll]ocked|[Pp]aired",cr$role),"Locked validation",
                 ifelse(grepl("[Ee]xternal",cr$role),"External sensitivity","Localisation")))
cr$role <- factor(cr$role_group,levels=c("Discovery","Locked validation","External sensitivity","Localisation"))
layer_cols <- c(Microbiome=unname(P["blue"]),Metabolome=unname(P["amber"]),
                `Bulk host`=unname(P["teal"]),`Single cell`=unname(P["violet"]))
draw1 <- function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(1:2,2,1),heights=c(1.3,1))
  base_par(c(3.3,8.2,1.1,1.0)); y<-rev(seq_len(nrow(cr))); x<-as.numeric(cr$role)
  plot(x,y,xlim=c(.55,4.45),ylim=c(.4,nrow(cr)+.6),xaxt="n",yaxt="n",xlab="",ylab="",pch=21,
       bg=layer_cols[cr$layer],col="white",cex=1.2+3.1*sqrt(cr$n/max(cr$n)))
  axis(1,1:4,levels(cr$role),cex.axis=.78);axis(2,y,cr$dataset,las=1,tick=FALSE,cex.axis=.72)
  abline(v=1:4,col=P["rule"],lty=3);tag("A")
  legend("bottomright",legend=names(layer_cols),pt.bg=layer_cols,pch=21,pt.cex=1.5,bty="n",cex=.70,ncol=2)
  text(x,y,cr$n,cex=.57,font=2,col=P["ink"])
  base_par(c(3.0,8.2,1.0,1.0)); roles<-levels(cr$role); layers<-names(layer_cols)
  m<-matrix(0,length(layers),length(roles),dimnames=list(layers,roles))
  for(i in seq_len(nrow(cr)))m[cr$layer[i],as.character(cr$role[i])]<-m[cr$layer[i],as.character(cr$role[i])]+cr$n[i]
  plot.new();plot.window(c(0,length(roles)),c(0,length(layers)))
  for(i in seq_along(layers))for(j in seq_along(roles)){
    v<-m[i,j];rect(j-1,length(layers)-i,j,length(layers)-i+1,col=if(v>0)alpha_col(layer_cols[layers[i]],.18+.65*v/max(m)) else "white",border=P["rule"])
    text(j-.5,length(layers)-i+.5,if(v>0) as.character(v) else "",font=if(v>0)2 else 1,cex=.78,col=if(v>0)P["ink"] else P["muted"])
  }
  axis(2,length(layers)-seq_along(layers)+.5,layers,las=1,tick=FALSE,cex.axis=.82);axis(3,seq_along(roles)-.5,roles,tick=FALSE,cex.axis=.80);tag("B")
}
save_pub(draw1,"Figure_1_CORE_RESULTS",125,12.5)

# Figure 2: full microbiome replication landscape
x <- read.csv(file.path(rd,"microbiome_cross_cohort","three_cohort_harmonised_genus_effects.csv"),check.names=FALSE)
sel <- truth(x$selected_in_discovery_q_lt_0_10); val <- truth(x$locked_validation_pass); ext <- truth(x$external_age_adjusted_pass)
draw2 <- function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(1,1),heights=c(1,1.05))
  base_par(c(3.2,3.5,1.1,1)); yy<--log10(pmax(x$PRJNA784025_discovery_p,1e-12));
  plot(x$PRJNA784025_discovery_estimate,yy,pch=21,bg=ifelse(sel,P["amber"],alpha_col(P["muted"],.45)),col="white",cex=ifelse(sel,1.1,.72),xlab="Discovery CLR effect",ylab="-log10(P)")
  abline(v=0,lty=2,col=P["rule"]);abline(h=-log10(.05),lty=3,col=P["rule"]);tag("A")
  lab<-order(yy,decreasing=TRUE)[1:4];text(x$PRJNA784025_discovery_estimate[lab],yy[lab],short(x$genus[lab],18),pos=ifelse(x$PRJNA784025_discovery_estimate[lab]>0,2,4),cex=.52)
  base_par(c(3.2,3.5,1.1,1)); cc<-ifelse(ext,P["green"],ifelse(val,P["teal"],ifelse(sel,P["amber"],alpha_col(P["muted"],.35))))
  plot(x$PRJNA784025_discovery_estimate,x$PRJNA784025_validation_estimate,pch=21,bg=cc,col="white",cex=.85,xlab="Discovery effect",ylab="Locked-validation effect")
  abline(h=0,v=0,lty=2,col=P["rule"]);abline(0,1,col=P["muted"],lty=3);tag("B")
  lab<-head(which(sel | val),8);text(x$PRJNA784025_discovery_estimate[lab],x$PRJNA784025_validation_estimate[lab],short(x$genus[lab],16),pos=ifelse(x$PRJNA784025_discovery_estimate[lab]>0,2,4),cex=.46)
  o<-order(abs(x$random_estimate),decreasing=TRUE)[1:18];z<-x[o,];y<-rev(seq_along(o))
  base_par(c(3.4,10.5,1.1,1));plot(z$random_estimate,y,xlim=range(c(z$hksj_ci_low_modified,z$hksj_ci_high_modified)),yaxt="n",pch=21,bg=P["teal"],col="white",xlab="HKSJ pooled CLR effect",ylab="")
  segments(z$hksj_ci_low_modified,y,z$hksj_ci_high_modified,y,col=P["muted"],lwd=1.35);abline(v=0,lty=2,col=P["rule"]);axis(2,y,wrap(z$genus,28),las=1,tick=FALSE,cex.axis=.70);tag("C")
}
save_pub(draw2,"Figure_2_CORE_RESULTS",170,12.5)

# Figure 3: metabolite effect and cross-specimen concordance
pl<-read.csv(file.path(fd,"Figure3_source_data_longitudinal_plasma_promotions.csv"),check.names=FALSE)
ti<-read.csv(file.path(fd,"Figure3_source_data_paired_tissue_promotions.csv"),check.names=FALSE)
ov<-read.csv(file.path(fd,"Figure3_source_data_exact_cross_specimen_overlap.csv"),check.names=FALSE)
draw3<-function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(1:4,2,2,byrow=TRUE))
  y<-rev(seq_len(nrow(pl)));base_par(c(3.1,8.0,1.1,1));plot(pl$effect_mid,y,xlim=range(c(pl$effect_min,pl$effect_max)),yaxt="n",pch=21,bg=P["amber"],xlab="Longitudinal coefficient",ylab="")
  segments(pl$effect_min,y,pl$effect_max,y,col=P["amber"],lwd=4);axis(2,y,wrap(pl$metabolite,26),las=1,tick=FALSE,cex.axis=.72);abline(v=0,lty=2,col=P["rule"]);tag("A")
  o<-order(ti$median_log2_difference);q<-ti[o,];y<-seq_len(nrow(q));base_par(c(3.1,8.2,1.1,1));plot(q$median_log2_difference,y,yaxt="n",pch=ifelse(truth(q$authentic_standard),23,21),bg=ifelse(truth(q$authentic_standard),P["green"],P["violet"]),xlab="Paired liver median log2 difference",ylab="")
  axis(2,y,wrap(q$metabolite,27),las=1,tick=FALSE,cex.axis=.69);abline(v=0,lty=2,col=P["rule"]);tag("B");legend("bottomright",c("Authentic standard","Putative level 2"),pch=c(23,21),pt.bg=c(P["green"],P["violet"]),bty="n",cex=.66)
  base_par(c(3.2,3.8,1.1,1));cc<-ifelse(truth(ov$direction_concordant),P["teal"],P["red"]);plot(ov$beta_group_provided_area_gee,ov$median_log2_difference,pch=21,bg=cc,col="white",cex=1.15,xlab="Plasma longitudinal coefficient",ylab="Paired liver effect")
  abline(h=0,v=0,lty=2,col=P["rule"]);dx<-c(-.012,.012,-.016,.016,-.012,.014,-.014);dy<-c(.02,-.03,.035,-.025,.02,.035,-.035);text(ov$beta_group_provided_area_gee+dx,ov$median_log2_difference+dy,seq_len(nrow(ov)),cex=.64,font=2);legend("bottomright",paste(seq_len(nrow(ov)),wrap(ov$name_key,20)),bty="n",cex=.49);tag("C")
  base_par(c(3.2,7.2,1.1,1));states<-cbind(Plasma=truth(ov$longitudinal_promotion),Tissue=truth(ov$analysis_pass),Direction=truth(ov$direction_concordant),Both_FDR=truth(ov$both_fdr_promoted));plot.new();plot.window(c(0,4),c(0,nrow(states)))
  for(i in seq_len(nrow(states)))for(j in 1:4)rect(j-1,nrow(states)-i,j,nrow(states)-i+1,col=if(states[i,j])P["teal"] else "white",border=P["rule"])
  axis(1,(1:4)-.5,c("Plasma","Tissue","Direction","Both FDR"),tick=FALSE,cex.axis=.66);axis(2,nrow(states)-(1:nrow(states))+.5,short(ov$name_key,18),las=1,tick=FALSE,cex.axis=.58);box(col=P["rule"]);tag("D")
}
save_pub(draw3,"Figure_3_CORE_RESULTS",155,11.5)

# Figure 4: interpretable axis-prioritisation evidence synthesis
ep<-read.csv(file.path(fd,"Figure4_source_data_edge_provenance.csv"),check.names=FALSE)
gs<-read.csv(file.path(fd,"Figure4_source_data_gate_scores.csv"),check.names=FALSE)
ws<-read.csv(file.path(fd,"Figure4_source_data_weight_sensitivity.csv"),check.names=FALSE)
draw4<-function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(c(1,1,2,3),2,2,byrow=TRUE),heights=c(1.05,1))
  gm<-as.matrix(gs[,2:9]);base_par(c(4.2,12.0,1.5,1));heat(gm,wrap(gs$candidate_family,36),c("Microbial","Metabolite","Microbe-\nmetabolite","Metabolite-\ntarget","Bulk host","Single cell","Contradiction","Robustness"),lim=1,values=TRUE);tag("A")
  base_par(c(4.0,7.8,1.2,1));meas<-aggregate(truth(ep$measured_in_project),list(Axis=ep$candidate_family),sum);total<-aggregate(ep$measured_in_project,list(Axis=ep$candidate_family),length);pv<-merge(data.frame(Axis=gs$candidate_family),merge(meas,total,by="Axis",all=TRUE),by="Axis",all.x=TRUE);names(pv)[2:3]<-c("Measured","Total");pv[is.na(pv)]<-0;pv$Curated<-pv$Total-pv$Measured;yy<-rev(seq_len(nrow(pv)));plot(pv$Total,yy,type="n",xlim=c(0,max(pv$Total)*1.15),yaxt="n",xlab="Number of traceable evidence edges",ylab="");segments(0,yy,pv$Measured,yy,col=P["teal"],lwd=12);segments(pv$Measured,yy,pv$Total,yy,col=P["amber"],lwd=12);points(pv$Total,yy,pch=21,bg=ifelse(pv$Total>0,P["amber"],"white"),cex=.8);axis(2,yy,wrap(pv$Axis,28),las=1,tick=FALSE,cex.axis=.66);legend("topright",c("Measured in project","Curated prior only"),col=c(P["teal"],P["amber"]),lwd=6,bty="n",cex=.64);tag("B")
  base_par(c(3.8,7.8,1.2,1));schemes<-unique(ws$weighting_scheme);cands<-unique(ws$candidate_family);plot(NA,xlim=c(1,length(schemes)),ylim=range(ws$weighted_score),xaxt="n",xlab="Weighting scheme",ylab="Exploratory score")
  for(i in seq_along(cands)){q<-ws[ws$candidate_family==cands[i],];q<-q[match(schemes,q$weighting_scheme),];lines(seq_along(schemes),q$weighted_score,col=c(P["blue"],P["teal"],P["amber"],P["violet"])[i],lwd=2);points(seq_along(schemes),q$weighted_score,pch=16,col=c(P["blue"],P["teal"],P["amber"],P["violet"])[i])}
  axis(1,seq_along(schemes),wrap(schemes,18),cex.axis=.69);for(i in seq_along(cands)){q<-ws[ws$candidate_family==cands[i],];q<-q[match(schemes,q$weighting_scheme),];text(length(schemes)-.03,tail(q$weighted_score,1),wrap(cands[i],25),pos=2,cex=.50,col=c(P["blue"],P["teal"],P["amber"],P["violet"])[i])};tag("C")
}
save_pub(draw4,"Figure_4_CORE_RESULTS",160,12.5)

# Figure 5: host replication with participant-level distributions
s5<-read.csv(file.path(rd,"integration","figure5_host_program_cross_dataset_source.csv"),check.names=FALSE)
tc<-read.csv(file.path(fd,"Figure5_source_data_TCGA_paired_differences.csv"),check.names=FALSE)
gse<-read.csv(file.path(rd,"GSE63898","host_programs","gse63898_program_validation.csv"),check.names=FALSE)
progs<-unique(s5$program);dsets<-unique(s5$dataset);m5<-matrix(NA,length(progs),length(dsets),dimnames=list(progs,dsets));for(i in seq_len(nrow(s5)))m5[s5$program[i],s5$dataset[i]]<-s5$effect[i]
draw5<-function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(1:4,2,2,byrow=TRUE),heights=c(1,1))
  base_par(c(4.2,10.0,1.1,.8));heat(m5,wrap(rownames(m5),30),colnames(m5),values=TRUE);tag("A")
  base_par(c(4.2,4.2,1.1,.8));lev<-unique(tc$program);cols<-c(P["blue"],P["teal"],P["amber"],P["violet"]);plot(NA,xlim=c(.5,length(lev)+.5),ylim=range(tc$paired_difference),xaxt="n",xlab="",ylab="TCGA paired difference")
  for(i in seq_along(lev)){v<-tc$paired_difference[tc$program==lev[i]];j<-seq(-.16,.16,length.out=length(v));points(i+j,v,pch=16,col=alpha_col(cols[i],.45),cex=.48);segments(i-.22,median(v),i+.22,median(v),col=cols[i],lwd=3)};axis(1,seq_along(lev),short(lev,14),las=2,cex.axis=.60);abline(h=0,lty=2,col=P["rule"]);tag("B")
  q<-gse[gse$family=="primary_PRJEB54571_GO_programs",];y<-rev(seq_len(nrow(q)));base_par(c(3.4,10.0,1.1,.8));plot(q$plate_adjusted_beta,y,xlim=range(c(q$plate_adjusted_ci_low,q$plate_adjusted_ci_high)),yaxt="n",pch=21,bg=ifelse(truth(q$strict_three_dataset_replication),P["green"],P["red"]),xlab="GSE63898 adjusted effect",ylab="");segments(q$plate_adjusted_ci_low,y,q$plate_adjusted_ci_high,y,col=P["muted"],lwd=1.4);axis(2,y,wrap(q$program,30),las=1,tick=FALSE,cex.axis=.72);abline(v=0,lty=2,col=P["rule"]);tag("C")
  base_par(c(3.2,3.8,1.1,.8));xr<-range(q$tcga_median_paired_difference);yr<-range(q$plate_adjusted_beta);plot(q$tcga_median_paired_difference,q$plate_adjusted_beta,xlim=xr+c(-.06,.06)*diff(xr),ylim=yr+c(-.10,.10)*diff(yr),pch=21,bg=ifelse(truth(q$strict_three_dataset_replication),P["green"],P["red"]),col="white",cex=1.25,xlab="TCGA paired effect",ylab="GSE63898 adjusted effect");abline(h=0,v=0,lty=2,col=P["rule"]);text(q$tcga_median_paired_difference,q$plate_adjusted_beta,wrap(q$program,18),pos=ifelse(q$tcga_median_paired_difference>mean(xr),2,4),cex=.52);tag("D")
}
save_pub(draw5,"Figure_5_CORE_RESULTS",165,13)

# Figure 7: leakage-free CRC component validation and framework behaviour
sp<-read.csv(file.path(rd,"second_scenario_crc_v2","species_leakage_free_results.csv"),check.names=FALSE)
repd<-truth(sp$heldout_reproduced);cand<-truth(sp$discovery_candidate)
bulk<-read.csv(file.path(rd,"second_scenario_crc","host_validation","bulk_program_replication.csv"),check.names=FALSE)
sc<-read.csv(file.path(rd,"second_scenario_crc_v2","host_validation","single_cell_whole_atlas_exact_signflip.csv"),check.names=FALSE)
err<-read.csv(file.path(rd,"framework_benchmark_v2","coding_error_summary.csv"),check.names=FALSE)
trans<-read.csv(file.path(rd,"framework_benchmark_v2","bidirectional_transition_summary.csv"),check.names=FALSE)
names(err)[1] <- "per_dimension_coding_error"
names(trans)[1] <- "workflow"
draw7<-function(){
  par(family=font,oma=c(.2,.2,.2,.2));layout(matrix(1:6,3,2,byrow=TRUE),widths=c(1,1),heights=c(1,1,1))
  base_par(c(3.2,3.6,1.1,.8));cc<-ifelse(repd,P["green"],ifelse(cand,P["amber"],alpha_col(P["muted"],.30)));plot(sp$rank_biserial_discovery,sp$rank_biserial_validation,pch=21,bg=cc,col="white",cex=ifelse(repd,1.25,.60),xlab="Discovery rank-biserial effect",ylab="Held-out validation effect");abline(h=0,v=0,lty=2,col=P["rule"]);abline(0,1,lty=3,col=P["muted"]);text(sp$rank_biserial_discovery[repd],sp$rank_biserial_validation[repd],short(sp$feature[repd],18),pos=3,cex=.52);tag("A")
  q<-sp[repd,];y<-rev(seq_len(nrow(q)));base_par(c(3.4,10.0,1.1,.8));plot(q$rank_biserial_validation,y,xlim=range(c(q$rank_biserial_validation_ci_low,q$rank_biserial_validation_ci_high)),yaxt="n",pch=22,bg=P["teal"],xlab="Held-out effect (95% bootstrap CI)",ylab="");segments(q$rank_biserial_validation_ci_low,y,q$rank_biserial_validation_ci_high,y,col=P["teal"],lwd=2);axis(2,y,wrap(q$feature,30),las=1,tick=FALSE,cex.axis=.72);tag("B")
  base_par(c(2.6,2.2,1.1,.8));labs<-c("Species","Orthologues","Faecal metab.","Bulk host","Whole-atlas scRNA");tested<-c(219,5141,299,3,3);passed<-c(4,0,0,1,0);plot(seq_along(labs),passed,pch=21,bg=ifelse(passed>0,P["green"],P["red"]),col="white",cex=2.2,xaxt="n",xlab="",ylab="Passed components",ylim=c(0,4.8));axis(1,seq_along(labs),labs,las=2,cex.axis=.58);text(seq_along(labs),passed,paste0(passed,"/",tested),cex=.58,font=2);tag("C")
  programs<-c("TLR4_NFKB_inflammatory_response","epithelial_EMT_barrier_response","autophagy_stress_response");mat<-matrix(NA,3,3,dimnames=list(c("TLR4/NF-kB","EMT/barrier","Autophagy/stress"),c("GSE74602","GSE39582","Whole-atlas scRNA")));for(i in 1:3){br<-bulk[bulk[[1]]==programs[i],];sr<-sc[sc[[1]]==programs[i],];mat[i,]<-c(br$GSE74602_effect,br$GSE39582_effect,sr$mean_tumor_minus_normal)};base_par(c(3.8,6.6,1.1,.8));heat(mat,rownames(mat),colnames(mat),values=TRUE);tag("D")
  base_par(c(3.2,3.7,1.1,.8));wf<-unique(err$workflow);cols<-c(P["red"],P["violet"],P["blue"],P["green"]);plot(NA,xlim=range(err$per_dimension_coding_error),ylim=c(0,1.03),xlab="Per-dimension coding error",ylab="False-promotion rate");for(i in seq_along(wf)){z<-err[err$workflow==wf[i],];lines(z$per_dimension_coding_error,z$mean_false_promotion_near_boundary,col=cols[i],lwd=2);points(z$per_dimension_coding_error,z$mean_false_promotion_near_boundary,pch=16,col=cols[i])};legend("right",gsub("_"," ",wf),col=cols,lwd=2,bty="n",cex=.55);tag("E")
  base_par(c(4.0,5.7,1.1,.8));z<-aggregate(decision_change_rate~workflow,trans,mean);dotchart(z$decision_change_rate,labels=gsub("_"," ",z$workflow),pch=21,bg=P["blue"],xlab="Mean bidirectional decision-change rate",cex=.62);abline(v=0,lty=2,col=P["rule"]);tag("F")
}
save_pub(draw7,"Figure_7_CORE_RESULTS",200,13)

cat(normalizePath(out,winslash="/"),"\n")
