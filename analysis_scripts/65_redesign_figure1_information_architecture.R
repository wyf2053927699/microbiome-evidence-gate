options(stringsAsFactors = FALSE)

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset=file.path(root, "07_figures", "q4_information_architecture_redesign"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font <- "Arial"
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
pal <- c(ink="#17242B", muted="#66747B", rule="#D8E0E3", pale="#F5F8F9",
         blue="#2E6F9E", blue_l="#E8F1F7", teal="#31877C", teal_l="#E5F3F0",
         amber="#B7791F", amber_l="#F8F0DE", violet="#6F5A9B", violet_l="#EEEAF6",
         red="#B44D4D", red_l="#F8E8E7", green="#397A57", green_l="#E7F2EB")

registry <- read.csv(file.path(root,"07_figures","Figure1_source_data_cohort_registry.csv"), check.names=FALSE)
gates <- read.csv(file.path(root,"06_results","framework_validation","operational_gate_specification_v1.csv"), check.names=FALSE)

layer_col <- c("Microbiome"=unname(pal["blue"]), "Metabolome"=unname(pal["amber"]),
               "Bulk host"=unname(pal["teal"]), "Single cell"=unname(pal["violet"]))
role_levels <- c("discovery","locked validation","external sensitivity","longitudinal candidate discovery",
                 "cross-specimen tissue audit","bulk discovery","paired validation","direct external validation","cell-type localisation")
role_group <- c("Discovery","Validation","Sensitivity","Discovery","Audit","Discovery","Validation","Validation","Localisation")

panel_label <- function(x,y,lab) text(x,y,lab,adj=c(0,1),font=2,cex=1.05,col=pal["ink"],xpd=NA)
round_box <- function(x0,y0,x1,y1,fill,border=NA,r=.012) {
  symbols((x0+x1)/2,(y0+y1)/2,rectangles=matrix(c(x1-x0,y1-y0),1),inches=FALSE,add=TRUE,bg=fill,fg=if(is.na(border)) fill else border)
}
wrap2 <- function(x, width=27) paste(strwrap(x,width=width),collapse="\n")

draw <- function() {
  par(family=font, mar=c(.2,.2,.2,.2), oma=c(0,0,0,0), fg=pal["ink"], xaxs="i", yaxs="i")
  plot.new(); plot.window(c(0,1),c(0,1))

  # A — compact evidence streams
  panel_label(.018,.988,"A")
  text(.067,.982,"Independent evidence streams",adj=c(0,1),font=2,cex=.88,col=pal["ink"])
  streams <- data.frame(layer=c("Microbiome","Metabolome","Bulk host","Single cell"), n=c(3,2,3,1),
                        unit=c("259 participants","243 participants / 652 samples","468 participants / 50 pairs","10 donors"))
  xs <- c(.06,.29,.52,.75); w <- .195
  for(i in 1:4){
    lc <- layer_col[streams$layer[i]]
    round_box(xs[i],.825,xs[i]+w,.935,paste0(lc,"18"),lc)
    rect(xs[i],.825,xs[i]+.012,.935,col=lc,border=NA)
    text(xs[i]+.026,.913,streams$layer[i],adj=c(0,1),font=2,cex=.79,col=pal["ink"])
    text(xs[i]+.026,.873,paste0(streams$n[i]," dataset",ifelse(streams$n[i]>1,"s","")),adj=c(0,.5),cex=.72,col=lc,font=2)
    text(xs[i]+.026,.844,streams$unit[i],adj=c(0,.5),cex=.62,col=pal["muted"])
  }
  text(.5,.792,"No cross-layer arrow denotes matched participants or mediation",cex=.66,col=pal["red"],font=2)
  segments(.06,.773,.94,.773,col=pal["rule"],lwd=.8)

  # B — hero role matrix
  panel_label(.018,.744,"B")
  text(.067,.738,"Dataset role registry",adj=c(0,1),font=2,cex=.88,col=pal["ink"])
  x_dataset <- .06; x_layer <- .32; x_unit <- .47; x_roles <- c(.70,.765,.83,.895)
  y0 <- .685; rh <- .045
  headers <- c("Dataset","Layer","Independent biological unit")
  for(j in 1:3) text(c(x_dataset,x_layer,x_unit)[j],y0+.026,headers[j],adj=c(0,.5),font=2,cex=.64,col=pal["muted"])
  for(j in 1:4) text(x_roles[j],y0+.026,c("Discover","Validate","Stress-test","Localise")[j],srt=35,adj=c(.5,.5),font=2,cex=.58,col=pal["muted"])
  for(i in seq_len(nrow(registry))){
    y <- y0-(i-.5)*rh
    if(i%%2==0) rect(.052,y-rh/2,.945,y+rh/2,col=pal["pale"],border=NA)
    lc <- layer_col[registry$layer[i]]
    text(x_dataset,y,registry$dataset[i],adj=c(0,.5),cex=.66,font=2,col=pal["ink"])
    rect(x_layer,y-.011,x_layer+.018,y+.011,col=lc,border=NA)
    text(x_layer+.026,y,registry$layer[i],adj=c(0,.5),cex=.62,col=pal["ink"])
    unit <- sub("; ","\n",registry$independent_unit[i],fixed=TRUE)
    text(x_unit,y,unit,adj=c(0,.5),cex=.57,col=pal["muted"])
    rg <- role_group[i]; idx <- match(rg,c("Discovery","Validation","Sensitivity","Localisation"))
    if(is.na(idx)) idx <- ifelse(rg=="Audit",3,1)
    for(j in 1:4){
      symbols(x_roles[j],y,circles=.0075,inches=FALSE,add=TRUE,bg=if(j==idx) lc else "white",fg=if(j==idx) lc else pal["rule"])
    }
  }
  rect(.052,y0-9*rh,.945,y0,col=NA,border=pal["rule"],lwd=.8)

  # C — executable gate sequence and claim classes
  panel_label(.018,.258,"C")
  text(.067,.252,"Ordered operational gates",adj=c(0,1),font=2,cex=.88,col=pal["ink"])
  gx <- seq(.07,.77,length.out=6); gy <- .165
  gate_short <- c("Metadata","Biological\nunit","Feature\nharmonisation","Multiplicity","Independent\nreplication","Claim\nstrength")
  for(i in 1:6){
    fill <- if(i<6) pal["blue_l"] else pal["green_l"]; br <- if(i<6) pal["blue"] else pal["green"]
    symbols(gx[i],gy,rectangles=matrix(c(.105,.072),1),inches=FALSE,add=TRUE,bg=fill,fg=br)
    text(gx[i],gy,gate_short[i],cex=.61,font=2,col=br)
    text(gx[i],gy-.052,paste0("G",i),cex=.55,col=pal["muted"])
    if(i<6) arrows(gx[i]+.058,gy,gx[i+1]-.058,gy,length=.035,angle=22,lwd=.8,col=pal["muted"])
  }
  arrows(.805,.165,.855,.165,length=.04,angle=22,lwd=1.0,col=pal["muted"])
  rect(.866,.112,.955,.218,col=pal["red_l"],border=pal["red"],lwd=.9)
  text(.910,.188,"Claim class",font=2,cex=.61,col=pal["red"])
  text(.910,.151,"association\nprior hypothesis\nmechanism",cex=.55,col=pal["ink"])
  text(.51,.066,"PASS is required at every essential gate; FAIL or MISSING cannot be compensated by another layer",cex=.62,font=2,col=pal["red"])
}

base <- file.path(outdir,"Figure1_evidence_architecture_Q4")
width_mm <- 180; height_mm <- 150; w <- width_mm/25.4; h <- height_mm/25.4
png(paste0(base,".png"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=9,family=font,bg="white"); draw(); dev.off()
tiff(paste0(base,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,compression="lzw",pointsize=9,family=font,bg="white"); draw(); dev.off()
cairo_pdf(paste0(base,".pdf"),width=w,height=h,pointsize=9,family=font); draw(); dev.off()
if(requireNamespace("svglite",quietly=TRUE)) {
  svglite::svglite(paste0(base,".svg"),width=w,height=h,pointsize=9,system_fonts=list(sans=font)); draw(); dev.off()
} else {
  svg(paste0(base,".svg"),width=w,height=h,pointsize=9,family=font); draw(); dev.off()
}
