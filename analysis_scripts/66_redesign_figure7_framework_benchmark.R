options(stringsAsFactors = FALSE)
root <- Sys.getenv("HCC_PROJECT_ROOT", unset="")
if(!nzchar(root)) root <- normalizePath(".",winslash="/",mustWork=TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR",unset=file.path(root,"07_figures","q4_information_architecture_redesign"))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

font <- "Arial"
if(.Platform$OS.type=="windows") windowsFonts(Arial=windowsFont("Arial"))
pal <- c(ink="#17242B",muted="#66747B",rule="#D8E0E3",pale="#F5F8F9",
         blue="#2E6F9E",blue_l="#E8F1F7",teal="#31877C",teal_l="#E5F3F0",
         amber="#B7791F",amber_l="#F8F0DE",violet="#6F5A9B",violet_l="#EEEAF6",
         red="#B44D4D",red_l="#F8E8E7",green="#397A57",green_l="#E7F2EB")

bd <- file.path(root,"06_results","framework_benchmark")
allp <- read.csv(file.path(bd,"benchmark_all_binary_patterns.csv"),check.names=FALSE)
obs <- read.csv(file.path(bd,"benchmark_observed_workflow_comparison.csv"),check.names=FALSE)
psum <- read.csv(file.path(bd,"benchmark_partial_code_summary.csv"),check.names=FALSE)
ssum <- read.csv(file.path(bd,"benchmark_structural_summary.csv"),check.names=FALSE)

methods <- c("any_positive","majority_dimensions","additive_score_ge_0_5","noncompensatory_all_complete")
method_lab <- c("Any positive","Majority","Additive >= 0.5","Non-compensatory")
short_candidate <- c("Bile acid / FXR","Aromatic / AHR","Indoles / AHR","Choline / lipid")
panel_label <- function(x,y,lab) text(x,y,lab,adj=c(0,1),font=2,cex=1.05,col=pal["ink"],xpd=NA)

draw <- function(){
  par(family=font,mar=c(.15,.15,.15,.15),oma=c(0,0,0,0),xaxs="i",yaxs="i",fg=pal["ink"])
  plot.new(); plot.window(c(0,1),c(0,1))

  # A: operational rule architecture
  panel_label(.018,.987,"A"); text(.065,.981,"Framework behaviour is benchmarked at three distinct levels",adj=c(0,1),font=2,cex=.87)
  cards <- data.frame(x0=c(.055,.365,.675),x1=c(.325,.635,.945),title=c("Component validation","Framework behaviour","External portability"),
                      body=c("Input contract\nunit / mapping / multiplicity","2^8 patterns\nstructural controls / perturbation","Independent disease context\nseparate Go/No-Go"),
                      fill=c(pal["blue_l"],pal["teal_l"],pal["violet_l"]),border=c(pal["blue"],pal["teal"],pal["violet"]))
  for(i in 1:3){
    rect(cards$x0[i],.847,cards$x1[i],.935,col=cards$fill[i],border=cards$border[i],lwd=.9)
    text(cards$x0[i]+.015,.916,cards$title[i],adj=c(0,1),font=2,cex=.69,col=cards$border[i])
    text(cards$x0[i]+.015,.878,cards$body[i],adj=c(0,.5),cex=.57,col=pal["ink"])
    if(i<3) arrows(cards$x1[i]+.008,.891,cards$x0[i+1]-.008,.891,length=.035,lwd=.8,col=pal["muted"])
  }
  text(.5,.824,"These tests answer different questions and are not interchangeable",cex=.60,font=2,col=pal["red"])

  # B: hero heatmap, all binary patterns ordered by completeness
  panel_label(.018,.785,"B"); text(.065,.779,"All 256 binary evidence patterns",adj=c(0,1),font=2,cex=.87)
  ord <- order(allp$n_complete,allp$pattern); z <- as.matrix(allp[ord,methods])
  # Reserve a dedicated title/header band.  The previous y1=.735 allowed the
  # two-line method headers to collide with the panel title at final size.
  x0 <- .075; x1 <- .58; y0 <- .385; y1 <- .700
  nr <- nrow(z); nc <- ncol(z); cw <- (x1-x0)/nc; rh <- (y1-y0)/nr
  for(j in 1:nc) for(i in 1:nr){
    fill <- if(z[i,j]==1) pal["red"] else pal["blue_l"]
    rect(x0+(j-1)*cw,y1-i*rh,x0+j*cw,y1-(i-1)*rh,col=fill,border=NA)
  }
  rect(x0,y0,x1,y1,border=pal["rule"],lwd=.8)
  for(j in 0:nc) segments(x0+j*cw,y0,x0+j*cw,y1,col="white",lwd=.8)
  cuts <- sapply(0:8,function(k) sum(allp$n_complete<=k));
  for(k in 0:8){
    yy <- y1-cuts[k+1]*rh
    segments(x0,yy,x1,yy,col="#FFFFFFAA",lwd=.45)
    mid <- if(k==0) cuts[1]/2 else (cuts[k]+cuts[k+1])/2
    text(x0-.012,y1-mid*rh,k,adj=c(1,.5),cex=.54,col=pal["muted"])
  }
  text(x0-.048,(y0+y1)/2,"Complete dimensions",srt=90,cex=.60,font=2)
  heat_lab <- c("Any\npositive","Majority\nrule","Additive\n>= 0.5","Non-\ncompensatory")
  for(j in 1:nc) text(x0+(j-.5)*cw,y1+.014,heat_lab[j],adj=c(.5,0),cex=.52,font=2,col=pal["ink"])
  rect(.095,.354,.115,.367,col=pal["blue_l"],border=pal["rule"]); text(.121,.360,"not promoted",adj=c(0,.5),cex=.55)
  rect(.235,.354,.255,.367,col=pal["red"],border=NA); text(.261,.360,"promoted",adj=c(0,.5),cex=.55)
  text(.418,.360,"Non-compensatory promotion occurs only for 11111111",adj=c(0,.5),cex=.57,font=2,col=pal["green"])

  # C: structural negative controls
  panel_label(.615,.785,"C"); text(.662,.779,"Structural false promotion",adj=c(0,1),font=2,cex=.87)
  vals <- ssum$structural_false_promotion_rate; names(vals) <- method_lab
  bx0 <- .68; bx1 <- .94; by <- c(.705,.635,.565,.495)
  for(i in 1:4){
    text(bx0-.012,by[i],method_lab[i],adj=c(1,.5),cex=.57,font=ifelse(i==4,2,1))
    segments(bx0,by[i],bx1,by[i],col=pal["rule"],lwd=6)
    segments(bx0,by[i],bx0+(bx1-bx0)*vals[i],by[i],col=if(i==4) pal["green"] else pal["red"],lwd=6)
    points(bx0+(bx1-bx0)*vals[i],by[i],pch=21,bg=if(i==4) pal["green"] else pal["red"],col="white",cex=1.2)
    text(bx1+.012,by[i],sprintf("%.0f%%",100*vals[i]),adj=c(0,.5),cex=.60,font=2,col=if(i==4) pal["green"] else pal["red"])
  }
  axis_y <- .455; segments(bx0,axis_y,bx1,axis_y,col=pal["ink"],lwd=.6)
  for(v in c(0,.5,1)){xx<-bx0+(bx1-bx0)*v;segments(xx,axis_y,xx,axis_y-.008);text(xx,axis_y-.018,paste0(v*100,"%"),adj=c(.5,1),cex=.52)}
  text((bx0+bx1)/2,.416,"8 single-dimension failures; complete case retained by all rules",cex=.55,col=pal["muted"])

  # D: partial-code perturbation
  panel_label(.018,.320,"D"); text(.065,.314,"Observed HCC candidates under 1,000 partial-code perturbations",adj=c(0,1),font=2,cex=.87)
  dx0 <- .125; dx1 <- .57; dy <- c(.255,.215,.175,.135)
  pord <- match(obs$candidate_family, psum$candidate_family)
  for(i in 1:4){
    rr <- psum[pord[i],]
    text(dx0-.012,dy[i],short_candidate[i],adj=c(1,.5),cex=.57)
    segments(dx0,dy[i],dx1,dy[i],col=pal["rule"],lwd=1)
    xx0 <- dx0+(dx1-dx0)*rr$additive_score_min/.6; xx1 <- dx0+(dx1-dx0)*rr$additive_score_max/.6
    segments(xx0,dy[i],xx1,dy[i],col=pal["amber"],lwd=5)
    abx <- dx0+(dx1-dx0)*obs$additive_score[i]/.6
    points(abx,dy[i],pch=21,bg=pal["amber_l"],col=pal["amber"],cex=1.05,lwd=1)
    text(dx1+.012,dy[i],paste0(rr$additive_promotions," / 1000 vs 0"),adj=c(0,.5),cex=.56,font=2,col=if(rr$additive_promotions>0) pal["red"] else pal["green"])
  }
  for(v in c(0,.25,.5)){xx<-dx0+(dx1-dx0)*v/.6;segments(xx,.105,xx,.112);text(xx,.097,sprintf("%.2f",v),adj=c(.5,1),cex=.52)}
  text((dx0+dx1)/2,.073,"Additive score range; right labels: additive promotions vs non-compensatory promotions",cex=.54,col=pal["muted"])

  # E: observed claim boundary
  panel_label(.615,.320,"E"); text(.662,.314,"Observed claim migration",adj=c(0,1),font=2,cex=.87)
  ex <- c(.69,.79,.89); labs <- c("Positive\nsignal","Integrated\ncandidate","Mechanistic\naxis")
  for(j in 1:3){
    rect(ex[j]-.042,.185,ex[j]+.042,.252,col=c(pal["blue_l"],pal["amber_l"],pal["red_l"])[j],border=c(pal["blue"],pal["amber"],pal["red"])[j],lwd=.8)
    text(ex[j],.219,labs[j],cex=.56,font=2,col=c(pal["blue"],pal["amber"],pal["red"])[j])
    if(j<3) arrows(ex[j]+.047,.219,ex[j+1]-.047,.219,length=.035,lwd=.8,col=pal["muted"],lty=2)
  }
  segments(.665,.145,.935,.145,col=pal["rule"],lwd=5)
  points(.79,.145,pch=21,bg=pal["amber"],col="white",cex=1.4)
  text(.79,.121,"4/4 observed families remain NO-GO",cex=.62,font=2,col=pal["red"])
  text(.80,.082,"Permitted: bounded associations and prior hypotheses\nNot permitted: participant-linked mediation or causal mechanism",cex=.56,col=pal["muted"])
  text(.5,.030,"Benchmark results describe framework behaviour, not clinical accuracy or biological truth",cex=.62,font=2,col=pal["red"])
}

base <- file.path(outdir,"Figure7_framework_benchmark_Q4")
width_mm<-180;height_mm<-165;w<-width_mm/25.4;h<-height_mm/25.4
png(paste0(base,".png"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=9,family=font,bg="white");draw();dev.off()
tiff(paste0(base,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,compression="lzw",pointsize=9,family=font,bg="white");draw();dev.off()
cairo_pdf(paste0(base,".pdf"),width=w,height=h,pointsize=9,family=font);draw();dev.off()
if(requireNamespace("svglite",quietly=TRUE)) {
  svglite::svglite(paste0(base,".svg"),width=w,height=h,pointsize=9,system_fonts=list(sans=font));draw();dev.off()
} else {
  svg(paste0(base,".svg"),width=w,height=h,pointsize=9,family=font);draw();dev.off()
}
