options(stringsAsFactors = FALSE)
root <- Sys.getenv("HCC_PROJECT_ROOT", unset="")
if(!nzchar(root)) root <- normalizePath(".", winslash="/", mustWork=TRUE)
outdir <- file.path(root,"07_figures","second_scenario_crc")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

base_family <- "Times New Roman"
base_size <- 10
width_mm <- 183
height_mm <- 190
font <- base_family # user-approved manuscript-wide publication font
if(.Platform$OS.type=="windows") windowsFonts(`Times New Roman`=windowsFont("Times New Roman"))
pal <- c(ink="#17242B",muted="#66747B",rule="#D8E0E3",pale="#F5F8F9",
         blue="#2E6F9E",blue_l="#E8F1F7",teal="#31877C",teal_l="#E5F3F0",
         amber="#B7791F",amber_l="#F8F0DE",violet="#6F5A9B",violet_l="#EEEAF6",
         red="#B44D4D",red_l="#F8E8E7",green="#397A57",green_l="#E7F2EB")

bd <- file.path(root,"06_results","second_scenario_crc")
sp <- read.csv(file.path(bd,"species_primary_results.csv"),check.names=FALSE)
sp <- sp[tolower(as.character(sp$reproduced))=="true",]
host <- read.csv(file.path(bd,"host_validation","bulk_program_replication.csv"),check.names=FALSE)
sc <- read.csv(file.path(bd,"host_validation","single_cell_program_exact_signflip.csv"),check.names=FALSE)
gate <- read.csv(file.path(bd,"framework_portability","crc_candidate_gate_matrix.csv"),check.names=FALSE)
rules <- read.csv(file.path(bd,"framework_portability","crc_rule_comparison_summary.csv"),check.names=FALSE)
sens <- read.csv(file.path(bd,"species_threshold_sensitivity.csv"),check.names=FALSE)
names(sp)[1] <- "family"
names(host)[1] <- "program"
names(sc)[1] <- "program"
names(gate)[1] <- "candidate_family"
names(rules)[1] <- "workflow"
names(sens)[1] <- "prevalence_threshold"

dims <- c("microbial_reproducibility","human_metabolite_observed","microbe_metabolite_link","metabolite_target_link",
          "host_reproduced_two_datasets","single_cell_localisation","major_contradiction_resolved","robust_to_analytic_choices")
dimlab <- c("Microbial\nreplication","Human\nmetabolite","Microbe–\nmetabolite","Metabolite–\nhost","Host ×2\nbulk","Single-cell\nlocalisation","Contradiction\nresolved","Robustness")
case_lab <- c("Inflammatory","Autophagy","EMT/barrier")

panel <- function(x,y,lab,title){text(x,y,lab,adj=c(0,1),font=2,cex=1.08,col=pal["ink"],xpd=NA);text(x+.045,y-.003,title,adj=c(0,1),font=2,cex=.88,col=pal["ink"],xpd=NA)}
box <- function(x0,y0,x1,y1,fill,border,title,body){rect(x0,y0,x1,y1,col=fill,border=border,lwd=1);text(x0+.012,y1-.014,title,adj=c(0,1),font=2,cex=.67,col=border);text(x0+.012,y0+.018,body,adj=c(0,0),cex=.58,col=pal["ink"])}

draw <- function(){
  par(family=font,mar=c(.2,.2,.2,.2),oma=c(0,0,0,0),xaxs="i",yaxs="i",fg=pal["ink"])
  plot.new(); plot.window(c(0,1),c(0,1))

  # A — outcome-blind selection and locked cohort
  panel(.02,.985,"A","Outcome-blind screening selected a fully executable CRC scenario")
  box(.055,.855,.245,.935,pal["pale"],pal["muted"],"4 candidates","CRC | IBD | MASLD | T2D")
  arrows(.252,.895,.292,.895,length=.035,lwd=1,col=pal["muted"])
  box(.30,.855,.49,.935,pal["green_l"],pal["green"],"CRC: GO (98/100)","All non-compensatory entry gates passed")
  arrows(.497,.895,.537,.895,length=.035,lwd=1,col=pal["muted"])
  box(.545,.855,.73,.935,pal["blue_l"],pal["blue"],"277 participants","127 healthy | 150 CRC")
  arrows(.737,.895,.777,.895,length=.035,lwd=1,col=pal["muted"])
  box(.785,.855,.965,.935,pal["violet_l"],pal["violet"],"Locked split","165 discovery · 112 validation")
  text(.51,.833,"Selection used access, mapping, layer completeness and reproducibility - not effect direction or significance",cex=.60,font=2,col=pal["red"])

  # B — hero replicated species effect plot
  panel(.02,.792,"B","Five CRC-enriched species reproduced in the locked validation set")
  spp <- sp$feature
  ord <- order(sp$rank_biserial_validation)
  spp <- spp[ord]; dd <- sp$rank_biserial_discovery[ord]; vv <- sp$rank_biserial_validation[ord]
  y <- seq(.715,.565,length.out=length(spp)); x0 <- .19; x1 <- .56
  segments(x0,.535,x1,.535,col=pal["ink"],lwd=.8)
  for(v in c(0,.2,.4)){xx<-x0+(x1-x0)*v/.5;segments(xx,.535,xx,.527);text(xx,.517,sprintf("%.1f",v),adj=c(.5,1),cex=.56)}
  segments(x0,y[1]-.025,x0,y[length(y)]+.025,col=pal["rule"],lwd=1)
  for(i in seq_along(spp)){
    text(x0-.014,y[i],spp[i],adj=c(1,.5),cex=.59,font=3)
    xd <- x0+(x1-x0)*dd[i]/.5; xv <- x0+(x1-x0)*vv[i]/.5
    segments(xd,y[i],xv,y[i],col=pal["rule"],lwd=2)
    points(xd,y[i],pch=21,bg=pal["blue"],col="white",cex=1.05,lwd=.8)
    points(xv,y[i],pch=22,bg=pal["teal"],col="white",cex=1.05,lwd=.8)
  }
  text((x0+x1)/2,.495,"Rank-biserial effect (CRC > healthy)",cex=.58,col=pal["muted"])
  points(.60,.714,pch=21,bg=pal["blue"],col="white",cex=1.05);text(.62,.714,"Discovery: FDR < 0.05",adj=c(0,.5),cex=.58)
  points(.60,.679,pch=22,bg=pal["teal"],col="white",cex=1.05);text(.62,.679,"Validation: P < 0.05",adj=c(0,.5),cex=.58)
  text(.60,.625,"All five plotted features passed\nthe frozen same-direction rule",adj=c(0,.5),cex=.63,font=2,col=pal["green"])
  retained_n <- sum(tolower(as.character(sens$all_primary_hits_retained))=="true")
  text(.60,.575,sprintf("Threshold robustness: %d / %d settings\nretained all five primary hits",retained_n,nrow(sens)),adj=c(0,.5),cex=.60,col=pal["muted"])

  # C — complete evidence ledger
  panel(.02,.455,"C","Positive components coexist with explicit missing bridges")
  fam <- c("Species","Microbial KO","Metabolites","Bulk host\nprograms","Single-cell\nprograms")
  tested <- c(217,5124,299,3,3); passed <- c(5,0,0,1,1)
  xx <- seq(.12,.77,length.out=5)
  for(i in 1:5){
    rect(xx[i]-.055,.337,xx[i]+.055,.405,col=if(passed[i]>0) pal["green_l"] else pal["pale"],border=if(passed[i]>0) pal["green"] else pal["rule"],lwd=1)
    text(xx[i],.384,passed[i],font=2,cex=.92,col=if(passed[i]>0) pal["green"] else pal["muted"])
    text(xx[i],.352,paste0("of ",tested[i]),cex=.54,col=pal["muted"])
    text(xx[i],.322,fam[i],cex=.56,font=2)
  }
  rect(.825,.325,.965,.405,col=pal["red_l"],border=pal["red"],lwd=1)
  text(.895,.384,"MISSING",font=2,cex=.70,col=pal["red"])
  text(.895,.348,"validated metabolite\nbridge",cex=.56)
  text(.50,.285,"Component positives are retained; absent links are not imputed and cannot be rescued by literature count",cex=.58,font=2,col=pal["red"])

  # D — host program triangulation
  panel(.02,.245,"D","Host datasets distinguish replication, localisation and conflict")
  progs <- c("TLR4_NFKB_inflammatory_response","epithelial_EMT_barrier_response","autophagy_stress_response")
  plab <- c("TLR4/NF-kB","EMT/barrier","Autophagy/stress")
  mat <- matrix(NA,3,3); sig <- matrix(FALSE,3,3)
  for(i in 1:3){
    rr<-host[host$program==progs[i],]; ss<-sc[sc$program==progs[i],]
    mat[i,]<-c(rr$GSE74602_effect,rr$GSE39582_effect,ss$mean_tumor_minus_normal)
    sig[i,]<-c(rr$GSE74602_q<.05,rr$GSE39582_q<.05,ss$exact_signflip_q<.05)
  }
  hx0<-.19; hy0<-.065; cw<-.105; ch<-.047
  # coordinates are shifted below through positive plot space
  hybase <- .075
  text(hx0+cw*(0:2)+cw/2,.192,c("GSE74602\npaired bulk","GSE39582\nbulk","GSE200997\npaired donors"),cex=.54,font=2)
  for(i in 1:3){
    yy <- .162-(i-1)*ch
    text(hx0-.012,yy,plab[i],adj=c(1,.5),cex=.57)
    for(j in 1:3){
      fill <- if(mat[i,j]>0) pal["blue_l"] else pal["amber_l"]
      border <- if(sig[i,j]) if(mat[i,j]>0) pal["blue"] else pal["amber"] else pal["rule"]
      rect(hx0+(j-1)*cw,yy-ch*.36,hx0+j*cw,yy+ch*.36,col=fill,border=border,lwd=if(sig[i,j]) 1.4 else .7)
      text(hx0+(j-.5)*cw,yy,sprintf("%+.2f",mat[i,j]),cex=.54,font=if(sig[i,j]) 2 else 1,col=if(sig[i,j]) border else pal["muted"])
    }
  }
  text(.515,.158,"Autophagy/stress",font=2,cex=.58,col=pal["green"],adj=c(0,.5));text(.515,.132,"replicated in two bulk cohorts",cex=.51,adj=c(0,.5))
  text(.515,.098,"EMT/barrier",font=2,cex=.58,col=pal["red"],adj=c(0,.5));text(.515,.072,"single-cell positive;",cex=.51,adj=c(0,.5));text(.515,.052,"bulk directions conflict",cex=.51,adj=c(0,.5))

  # E — actual CRC gate and rule comparison
  text(.69,.245,"E",adj=c(0,1),font=2,cex=1.08,col=pal["ink"],xpd=NA)
  text(.82,.222,"Unchanged gates block\npartial evidence",adj=c(.5,1),font=2,cex=.55,col="#17232B",xpd=NA)
  gx0<-.735; gx1<-.965; gy_bottom<-.055; gy_top<-.165
  z <- as.matrix(gate[,dims]); nr<-nrow(z); nc<-ncol(z); cwg<-(gx1-gx0)/nc; rh<-(gy_top-gy_bottom)/nr
  for(i in 1:nr) for(j in 1:nc){
    v<-z[i,j]; fill<-if(v==1) pal["green"] else if(v>=.5) pal["amber"] else pal["blue_l"]
    y_hi <- gy_top-(i-1)*rh; y_lo <- gy_top-i*rh
    rect(gx0+(j-1)*cwg,y_lo,gx0+j*cwg,y_hi,col=fill,border="white")
    text(gx0+(j-.5)*cwg,(y_lo+y_hi)/2,format(v,nsmall=1),cex=.43,col=if(v>=.5) "white" else pal["muted"])
  }
  for(j in 1:nc) text(gx0+(j-.5)*cwg,gy_bottom-.012,dimlab[j],srt=55,adj=c(1,.5),cex=.39)
  for(i in 1:nr) text(gx0-.009,gy_top-(i-.5)*rh,case_lab[i],adj=c(1,.5),cex=.50)
  promo <- setNames(rules$observed_test_cases_promoted,rules$workflow)
  text(.735,.198,sprintf("Promoted: any %d/3 | majority %d/3 | additive %d/3",promo["any_positive"],promo["majority_dimensions"],promo["additive_score_ge_0_5"]),adj=c(0,.5),cex=.44,col=pal["red"])
  text(.735,.181,"Non-compensatory: 0/3; stable across 24 ablations",adj=c(0,.5),cex=.46,font=2,col=pal["green"])
}

base <- file.path(outdir,"Figure7_CRC_external_portability")
svglite::svglite(paste0(base,".svg"),width=width_mm/25.4,height=height_mm/25.4,pointsize=base_size,system_fonts=list(serif="Times New Roman"));draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=width_mm/25.4,height=height_mm/25.4,pointsize=base_size,family="Times New Roman");draw();dev.off()
raster_tiff <- file.path(tempdir(),"Figure7_CRC_external_portability.tiff")
raster_png <- file.path(tempdir(),"Figure7_CRC_external_portability.png")
ragg::agg_tiff(raster_tiff,width=width_mm/25.4,height=height_mm/25.4,units="in",res=600,compression="lzw",pointsize=base_size);draw();dev.off()
ragg::agg_png(raster_png,width=width_mm/25.4,height=height_mm/25.4,units="in",res=300,pointsize=base_size);draw();dev.off()
stopifnot(file.copy(raster_tiff,paste0(base,".tiff"),overwrite=TRUE))
stopifnot(file.copy(raster_png,paste0(base,".png"),overwrite=TRUE))
cat(base,"\n")
