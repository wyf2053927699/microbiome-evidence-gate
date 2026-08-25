options(stringsAsFactors = FALSE)
root <- Sys.getenv("HCC_PROJECT_ROOT", unset="")
if (!nzchar(root)) root <- normalizePath(".", winslash="/", mustWork=TRUE)
outdir <- file.path(root, "07_figures", "second_scenario_crc_v2")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

font <- "Times New Roman"
if (.Platform$OS.type == "windows") windowsFonts(`Times New Roman`=windowsFont("Times New Roman"))
width_mm <- 183; height_mm <- 170; pointsize <- 12
pal <- c(ink="#17242B", muted="#5D6B72", rule="#CFD8DC", pale="#F4F7F8",
         blue="#2F6F9F", blue_l="#E7F0F7", teal="#2E8378", teal_l="#E2F1EE",
         amber="#B4771C", amber_l="#F8EEDB", red="#AA4747", red_l="#F7E5E4",
         green="#357554", green_l="#E4F0E8", violet="#6B5795", violet_l="#ECE8F4")

v2 <- file.path(root, "06_results", "second_scenario_crc_v2")
sp <- read.csv(file.path(v2, "species_leakage_free_results.csv"), check.names=FALSE)
sp <- sp[tolower(as.character(sp$heldout_reproduced)) == "true", ]
bulk <- read.csv(file.path(root, "06_results", "second_scenario_crc", "host_validation", "bulk_program_replication.csv"), check.names=FALSE)
sc <- read.csv(file.path(v2, "host_validation", "single_cell_whole_atlas_exact_signflip.csv"), check.names=FALSE)
err <- read.csv(file.path(root, "06_results", "framework_benchmark_v2", "coding_error_summary.csv"), check.names=FALSE)
trans <- read.csv(file.path(root, "06_results", "framework_benchmark_v2", "bidirectional_transition_summary.csv"), check.names=FALSE)

panel_title <- function(x, y, tag, title) {
  text(x, y, tag, adj=c(0,1), font=2, cex=1.12, col=pal["ink"], xpd=NA)
  text(x+0.042, y-0.002, title, adj=c(0,1), font=2, cex=0.84, col=pal["ink"], xpd=NA)
}
card <- function(x0,y0,x1,y1,fill,border,head,body) {
  rect(x0,y0,x1,y1,col=fill,border=border,lwd=1)
  text((x0+x1)/2,y1-0.018,head,font=2,cex=.67,col=border)
  text((x0+x1)/2,y0+0.021,body,cex=.58,col=pal["ink"])
}

draw <- function() {
  par(family=font, mar=c(.2,.2,.2,.2), oma=c(0,0,0,0), xaxs="i", yaxs="i", fg=pal["ink"])
  plot.new(); plot.window(c(0,1),c(0,1))

  panel_title(.018,.985,"A","Outcome-blind selection and leakage-free split")
  card(.055,.875,.245,.947,pal["pale"],pal["muted"],"4 candidates","CRC | IBD | MASLD | T2D")
  arrows(.255,.911,.292,.911,length=.035,lwd=1,col=pal["muted"])
  card(.30,.875,.485,.947,pal["green_l"],pal["green"],"CRC selected","All entry gates available")
  arrows(.495,.911,.532,.911,length=.035,lwd=1,col=pal["muted"])
  card(.54,.875,.725,.947,pal["blue_l"],pal["blue"],"277 participants","127 healthy | 150 CRC")
  arrows(.735,.911,.772,.911,length=.035,lwd=1,col=pal["muted"])
  card(.78,.875,.965,.947,pal["violet_l"],pal["violet"],"Held-out split","165 discovery | 112 validation")
  text(.51,.852,"Prevalence filtering and feature selection used discovery data only",cex=.61,font=2,col=pal["red"])

  panel_title(.018,.810,"B","Four species pass multiplicity-controlled held-out validation")
  sp <- sp[order(sp$rank_biserial_validation),]
  labs <- gsub(" ", "\n", sp$feature, fixed=TRUE)
  y <- seq(.716,.592,length.out=nrow(sp)); x0 <- .205; x1 <- .655
  scale_x <- function(v) x0 + (x1-x0)*(v/0.60)
  for (v in seq(0,.6,.2)) {
    xx <- scale_x(v); segments(xx,.557,xx,.565,col=pal["ink"]); text(xx,.542,sprintf("%.1f",v),cex=.58)
  }
  segments(x0,.565,x1,.565,lwd=.8,col=pal["ink"])
  for (i in seq_len(nrow(sp))) {
    text(x0-.018,y[i],sp$feature[i],adj=c(1,.5),font=3,cex=.62)
    segments(scale_x(sp$rank_biserial_discovery_ci_low[i]),y[i]+.014,scale_x(sp$rank_biserial_discovery_ci_high[i]),y[i]+.014,col=pal["blue"],lwd=1.6)
    points(scale_x(sp$rank_biserial_discovery[i]),y[i]+.014,pch=21,bg=pal["blue"],col="white",cex=1.18,lwd=.8)
    segments(scale_x(sp$rank_biserial_validation_ci_low[i]),y[i]-.014,scale_x(sp$rank_biserial_validation_ci_high[i]),y[i]-.014,col=pal["teal"],lwd=1.6)
    points(scale_x(sp$rank_biserial_validation[i]),y[i]-.014,pch=22,bg=pal["teal"],col="white",cex=1.18,lwd=.8)
  }
  text((x0+x1)/2,.520,"Rank-biserial effect (CRC > healthy); 95% bootstrap CI",cex=.59,col=pal["muted"])
  points(.715,.742,pch=21,bg=pal["blue"],col="white",cex=1.15); text(.738,.742,"Discovery",adj=c(0,.5),cex=.61)
  points(.715,.706,pch=22,bg=pal["teal"],col="white",cex=1.15); text(.738,.706,"Held-out validation",adj=c(0,.5),cex=.61)
  card(.705,.605,.955,.674,pal["green_l"],pal["green"],"Robustness","CLR + covariates + 9/9 settings")
  text(.83,.574,"6 candidates tested in validation; Solobacterium q=0.054",cex=.55,col=pal["red"])

  panel_title(.018,.480,"C","Component ledger keeps missing bridges visible")
  fam <- c("Species","Microbial\northologues","Faecal\nmetabolites","Bulk host\nprogrammes","Whole-atlas\nsingle cell")
  tested <- c(219,5141,299,3,3); passed <- c(4,0,0,1,0); xx <- seq(.10,.70,length.out=5)
  for(i in seq_along(fam)) {
    fill <- if(passed[i]>0) pal["green_l"] else pal["pale"]; border <- if(passed[i]>0) pal["green"] else pal["rule"]
    rect(xx[i]-.052,.359,xx[i]+.052,.430,col=fill,border=border,lwd=1)
    text(xx[i],.405,passed[i],font=2,cex=.92,col=if(passed[i]>0) pal["green"] else pal["muted"])
    text(xx[i],.375,paste0("of ",tested[i]),cex=.56,col=pal["muted"])
    text(xx[i],.334,fam[i],cex=.57,font=2)
  }
  rect(.79,.340,.965,.430,col=pal["red_l"],border=pal["red"],lwd=1.1)
  text(.878,.404,"STOP",font=2,cex=.80,col=pal["red"])
  text(.878,.365,"No reproduced\nmetabolite bridge",cex=.58)
  text(.51,.302,"Single-cell = PARTIAL: donor-aware whole-atlas evidence, no validated cell-type annotation",cex=.58,font=2,col=pal["amber"])

  panel_title(.018,.260,"D","Host evidence supports bulk replication, not cell localisation")
  programs <- c("TLR4_NFKB_inflammatory_response","epithelial_EMT_barrier_response","autophagy_stress_response")
  labels <- c("TLR4/NF-kB","EMT/barrier","Autophagy/stress")
  mat <- matrix(NA,3,3); qmat <- matrix(NA,3,3)
  for(i in 1:3) {
    br <- bulk[bulk[[1]]==programs[i],]; sr <- sc[sc[[1]]==programs[i],]
    mat[i,] <- c(br$GSE74602_effect,br$GSE39582_effect,sr$mean_tumor_minus_normal)
    qmat[i,] <- c(br$GSE74602_q,br$GSE39582_q,sr$exact_signflip_q_three_programs)
  }
  xstart <- .155; cw <- .12; rh <- .045
  text(xstart+cw*(0:2)+cw/2,.217,c("GSE74602\nbulk","GSE39582\nbulk","GSE200997\nwhole atlas"),cex=.55,font=2)
  for(i in 1:3) {
    yy <- .175-(i-1)*rh; text(xstart-.012,yy,labels[i],adj=c(1,.5),cex=.59)
    for(j in 1:3) {
      sig <- qmat[i,j] < .05; fill <- if(mat[i,j]>=0) pal["blue_l"] else pal["amber_l"]
      bord <- if(sig) if(mat[i,j]>=0) pal["blue"] else pal["amber"] else pal["rule"]
      rect(xstart+(j-1)*cw,yy-.017,xstart+j*cw,yy+.017,col=fill,border=bord,lwd=if(sig) 1.5 else .8)
      text(xstart+(j-.5)*cw,yy,sprintf("%+.2f",mat[i,j]),cex=.57,font=if(sig) 2 else 1,col=if(sig) bord else pal["muted"])
    }
  }
  text(.535,.190,"Bulk autophagy/stress",font=2,cex=.61,col=pal["green"],adj=c(0,.5))
  text(.535,.158,"replicates in both cohorts",cex=.56,adj=c(0,.5))
  text(.535,.118,"Single cell: 0/3 at BH q<0.05",font=2,cex=.60,col=pal["red"],adj=c(0,.5))
  text(.535,.087,"7 paired donors; exact sign-flip",cex=.55,adj=c(0,.5))

  panel_title(.705,.260,"E","Decision-boundary trade-off")
  workflows <- c("any_positive","majority_dimensions","additive_score_ge_0_5","noncompensatory_all_complete")
  wlabs <- c("Any positive","Majority","Additive","Non-compensatory")
  cols <- c(pal["red"],pal["violet"],pal["blue"],pal["green"])
  ex <- c(.01,.05,.10,.20); px0 <- .745; px1 <- .960; py0 <- .075; py1 <- .185
  sx <- function(v) px0+(px1-px0)*(v-.01)/.19; sy <- function(v) py0+(py1-py0)*v
  segments(px0,py0,px1,py0); segments(px0,py0,px0,py1)
  for(v in ex){xx<-sx(v);segments(xx,py0,xx,py0-.006);text(xx,py0-.015,paste0(round(v*100),"%"),cex=.48)}
  for(v in c(0,.5,1)){yy<-sy(v);segments(px0-.006,yy,px0,yy);text(px0-.010,yy,sprintf("%.1f",v),adj=c(1,.5),cex=.48)}
  for(k in seq_along(workflows)) {
    z <- err[err$workflow==workflows[k],]; z <- z[order(z[[1]]),]
    lines(sx(z[[1]]),sy(z$mean_false_promotion_near_boundary),col=cols[k],lwd=1.8)
    points(sx(z[[1]]),sy(z$mean_false_promotion_near_boundary),pch=16,col=cols[k],cex=.65)
  }
  text((px0+px1)/2,py0-.038,"Per-dimension coding error",cex=.50)
  text(px0-.050,(py0+py1)/2,"False promotion",srt=90,cex=.50)
  text(.850,.224,"Near-boundary incomplete cases",cex=.50,font=2)
  legend(.710,.213,legend=wlabs,col=cols,lwd=1.8,pch=16,cex=.43,bty="n",ncol=2,x.intersp=.4,y.intersp=.7)
}

base <- file.path(outdir, "Figure7_CRC_leakage_free_framework_benchmark")
grDevices::svg(paste0(base,".svg"),width=width_mm/25.4,height=height_mm/25.4,pointsize=pointsize,family=font); draw(); dev.off()
cairo_pdf(paste0(base,".pdf"),width=width_mm/25.4,height=height_mm/25.4,pointsize=pointsize,family=font); draw(); dev.off()
grDevices::tiff(paste0(base,".tiff"),width=width_mm/25.4,height=height_mm/25.4,units="in",res=600,compression="lzw",pointsize=pointsize,family=font); draw(); dev.off()
grDevices::png(paste0(base,".png"),width=width_mm/25.4,height=height_mm/25.4,units="in",res=300,pointsize=pointsize,family=font); draw(); dev.off()
cat(base,"\n")
