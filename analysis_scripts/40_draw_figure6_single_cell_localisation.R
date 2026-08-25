options(stringsAsFactors = FALSE)
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))

root <- Sys.getenv("HCC_PROJECT_ROOT", unset = "")
if (!nzchar(root)) root <- normalizePath(".", winslash = "/", mustWork = TRUE)
outdir <- Sys.getenv("HCC_OUTPUT_DIR", unset = "")
if (!nzchar(outdir)) outdir <- file.path(root, "07_figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(file.path(root, "06_results", "single_cell", "gse149614_celltype_counts_by_donor_site.csv"), check.names = FALSE)
results <- read.csv(file.path(root, "06_results", "single_cell", "paired_programs", "prespecified_program_paired_results.csv"), check.names = FALSE)
scores <- read.csv(file.path(root, "06_results", "single_cell", "paired_programs", "prespecified_program_scores.csv"), check.names = FALSE)
donor_loo <- read.csv(file.path(root, "06_results", "single_cell", "paired_programs", "prespecified_program_leave_one_donor_out.csv"), check.names = FALSE)
exact_sign <- read.csv(file.path(root, "06_results", "single_cell", "paired_programs", "prespecified_program_exact_signflip.csv"), check.names = FALSE)

celltypes <- c("B", "Endothelial", "Fibroblast", "Hepatocyte", "Myeloid", "T/NK")
programs <- c("hepatocyte_bile_acid_synthesis_transport", "hepatocyte_bile_acid_receptor_response",
              "ahr_aromatic_response", "choline_phospholipid_metabolism",
              "type_i_interferon_response", "cell_cycle")
prog_labs <- c("BA synthesis/\ntransport", "BA receptor\nresponse", "AHR/aromatic",
               "Choline/\nphospholipid", "Type I IFN", "Cell cycle")
names(prog_labs) <- programs
cell_cols <- c("B"="#999999", "Endothelial"="#76A5AF", "Fibroblast"="#C9A66B",
               "Hepatocyte"="#4477AA", "Myeloid"="#CC6677", "T/NK"="#228833")
up_col <- "#3A7D6B"; down_col <- "#C95A5A"; null_col <- "#B8B8B8"

# Derive donor-level tumour-minus-normal differences only for valid paired celltype-program tests.
paired_rows <- list(); k <- 1
for (i in seq_len(nrow(results))) {
  ct <- results$celltype[i]; pr <- results$program[i]
  z <- scores[scores$celltype == ct & scores$program == pr, c("patient", "site", "module_score")]
  wide <- reshape(z, idvar = "patient", timevar = "site", direction = "wide")
  if (all(c("module_score.Tumor", "module_score.Normal") %in% names(wide))) {
    wide <- wide[is.finite(wide$module_score.Tumor) & is.finite(wide$module_score.Normal), ]
    if (nrow(wide)) {
      paired_rows[[k]] <- data.frame(patient = wide$patient, celltype = ct, program = pr,
                                     paired_difference = wide$module_score.Tumor - wide$module_score.Normal)
      k <- k + 1
    }
  }
}
paired <- do.call(rbind, paired_rows); rownames(paired) <- NULL
write.csv(paired, file.path(outdir, "Figure6_source_data_donor_paired_differences.csv"), row.names = FALSE)
write.csv(donor_loo, file.path(outdir, "Figure6_source_data_leave_one_donor_out.csv"), row.names = FALSE)
write.csv(exact_sign, file.path(outdir, "Figure6_source_data_exact_signflip.csv"), row.names = FALSE)

panel_label <- function(lab) mtext(lab, side = 3, line = 0.7, at = par("usr")[1], adj = 0, font = 2, cex = 1.05)

draw_figure <- function() {
  layout(matrix(c(1,1,2,3,4,5), nrow=3, byrow=TRUE), heights=c(0.52,1.05,1.02), widths=c(1.03,1))
  par(family="Arial", fg="#222222", col.axis="#222222", col.lab="#222222")

  # a: design registry
  par(mar=c(1.0,1.0,2.0,1.0))
  plot(NA,xlim=c(0.5,4.5),ylim=c(0.05,1),axes=FALSE,xlab="",ylab="")
  vals <- c("71,915\ncells", "10 patients\n8 paired", "6 broad\ncell types", "Donor-aware\npseudobulk")
  subtitles <- c("post-QC atlas", "tumour-normal design", "prespecified classes", "cell-type-restricted scoring")
  for(i in 1:4){
    rect(i-0.42,0.16,i+0.42,0.86,col="#F4F7F8",border="#6F8794",lwd=1.2)
    text(i,0.59,vals[i],font=2,cex=0.80,col="#365A6C")
    text(i,0.29,subtitles[i],cex=0.58,col="#555555")
  }
  title("Single-cell analysis preserves the donor as the independent unit",adj=0.07,cex.main=0.93,font.main=2,line=0.55)
  panel_label("A")

  # b: cell composition by site
  par(mar=c(4.8,4.5,2.0,1.0))
  agg <- aggregate(n_cells~site+celltype,counts,sum)
  sites <- c("Normal","Tumor")
  mat <- matrix(0,nrow=length(celltypes),ncol=2,dimnames=list(celltypes,sites))
  for(i in seq_len(nrow(agg))) if(agg$site[i] %in% sites && agg$celltype[i] %in% celltypes) mat[agg$celltype[i],agg$site[i]]<-agg$n_cells[i]
  prop <- sweep(mat,2,colSums(mat),"/")
  barplot(prop,col=cell_cols[celltypes],border=NA,las=1,ylim=c(0,1),names.arg=sites,cex.names=0.72,ylab="Cell fraction")
  legend("top",legend=celltypes,fill=cell_cols[celltypes],bty="n",cex=0.57,ncol=3)
  title("Atlas composition by tissue site",adj=0.10,cex.main=0.90,font.main=2,line=0.55)
  mtext("Composition is descriptive; programme tests use donor pseudobulk",side=1,line=3.6,cex=0.56,col="#666666")
  panel_label("B")

  # c: all prespecified programme tests
  par(mar=c(6.5,6.4,2.0,1.0))
  m <- matrix(NA_real_,nrow=length(celltypes),ncol=length(programs),dimnames=list(celltypes,programs))
  q <- m
  for(i in seq_len(nrow(results))){
    if(results$celltype[i] %in% celltypes && results$program[i] %in% programs){
      m[results$celltype[i],results$program[i]]<-results$median_difference[i]
      q[results$celltype[i],results$program[i]]<-results$t_q_all_modules[i]
    }
  }
  lim <- max(abs(m),na.rm=TRUE)
  pal <- colorRampPalette(c("#C95A5A","#F4F4F4","#3A7D6B"))(101)
  plot(NA,xlim=c(0.5,6.5),ylim=c(0.5,6.5),axes=FALSE,xlab="",ylab="")
  axis(1,at=1:6,labels=prog_labs[programs],tick=FALSE,cex.axis=0.60)
  axis(2,at=6:1,labels=celltypes,las=2,tick=FALSE,cex.axis=0.67)
  for(i in 1:6) for(j in 1:6){
    yy<-7-i
    if(is.finite(m[i,j])){
      idx<-round((m[i,j]+lim)/(2*lim)*100)+1; idx<-max(1,min(101,idx))
      rect(j-0.47,yy-0.46,j+0.47,yy+0.46,col=pal[idx],border="white")
      text(j,yy,sprintf("%.2g",m[i,j]),cex=0.55,font=ifelse(q[i,j]<0.05,2,1),col=ifelse(abs(m[i,j])>0.65*lim,"white","#333333"))
      if(q[i,j]<0.05) points(j+0.34,yy+0.30,pch=8,cex=0.65,col="#222222")
    } else rect(j-0.47,yy-0.46,j+0.47,yy+0.46,col="#EEEEEE",border="white")
  }
  box(col="#BBBBBB")
  title("Tumour-minus-normal programme shifts",adj=0.10,cex.main=0.90,font.main=2,line=0.55)
  mtext("Numbers are median paired differences; star = global BH q<0.05",side=1,line=5.2,cex=0.56,col="#666666")
  panel_label("C")

  # d: donor-level differences for three paired-t BH signals, qualified by nonparametric sensitivities
  par(mar=c(6.1,4.6,2.0,1.0))
  sig <- results[results$t_q_all_modules<0.05,]
  sig$key <- paste(sig$celltype,sig$program,sep="|")
  sig_keys <- c("Hepatocyte|type_i_interferon_response","Myeloid|cell_cycle","T/NK|cell_cycle")
  sig <- sig[match(sig_keys,sig$key),]
  exact_sign$key <- paste(exact_sign$celltype, exact_sign$program, sep="|")
  exact_sig <- exact_sign[match(sig_keys, exact_sign$key),]
  dlabs <- c("Hepatocyte\nType I IFN","Myeloid\ncell cycle","T/NK\ncell cycle")
  pv <- paired[paste(paired$celltype,paired$program,sep="|") %in% sig_keys,]
  ylim<-range(pv$paired_difference,finite=TRUE)
  plot(NA,xlim=c(0.5,3.5),ylim=ylim,axes=FALSE,xlab="",ylab="Tumour - normal pseudobulk score")
  axis(1,at=1:3,labels=dlabs,tick=FALSE,cex.axis=0.66); axis(2,las=1,cex.axis=0.64); abline(h=0,lty=2,col="#999999")
  set.seed(20260818)
  for(i in 1:3){
    z<-pv[paste(pv$celltype,pv$program,sep="|")==sig_keys[i],]
    cc<-if(median(z$paired_difference)>0) up_col else down_col
    points(jitter(rep(i,nrow(z)),amount=0.10),z$paired_difference,pch=16,cex=0.60,col=adjustcolor(cc,0.55))
    segments(i-0.22,median(z$paired_difference),i+0.22,median(z$paired_difference),lwd=2.6,col=cc)
    text(i,ylim[2]-0.05*diff(ylim),paste0("paired-t q=",format(sig$t_q_all_modules[i],digits=2,scientific=TRUE)),cex=0.54)
    text(i,ylim[2]-0.14*diff(ylim),paste0("exact q=",format(exact_sig$exact_signflip_q_global_bh[i],digits=2)),cex=0.52,col="#9E2F2F")
  }
  title("Paired-t signals are not globally confirmed by exact sign-flip",adj=0.11,cex.main=0.82,font.main=2,line=0.55)
  mtext("Eight paired donors per test; Wilcoxon and exact sign-flip BH q=0.086",side=1,line=4.6,cex=0.56,col="#9E2F2F")
  panel_label("D")

  # e: axis-specific null ledger
  par(mar=c(4.4,8.7,2.0,1.0))
  keep <- (results$celltype=="Hepatocyte" & results$program %in% c("hepatocyte_bile_acid_synthesis_transport","hepatocyte_bile_acid_receptor_response","ahr_aromatic_response","choline_phospholipid_metabolism")) |
          (results$celltype=="Myeloid" & results$program %in% c("ahr_aromatic_response","choline_phospholipid_metabolism")) |
          (results$celltype=="T/NK" & results$program %in% c("ahr_aromatic_response","choline_phospholipid_metabolism"))
  z<-results[keep,]
  labs<-paste(z$celltype,prog_labs[z$program],sep=" | ")
  ord<-order(z$celltype,z$program); z<-z[ord,]; labs<-labs[ord]
  xr<-range(c(z$median_difference,0),finite=TRUE)+c(-0.15,0.15)
  plot(NA,xlim=xr,ylim=c(0.2,nrow(z)+1.2),axes=FALSE,xlab="Median tumour-normal difference",ylab="")
  axis(1,cex.axis=0.62); axis(2,at=rev(seq_len(nrow(z))),labels=labs,las=2,tick=FALSE,cex.axis=0.56); abline(v=0,lty=2,col="#999999")
  for(i in seq_len(nrow(z))){
    yy<-nrow(z)-i+1; points(z$median_difference[i],yy,pch=21,bg=null_col,col="white",cex=1.05)
    text(xr[2],yy,paste0("q=",format(z$t_q_all_modules[i],digits=2)),adj=1,cex=0.51,col="#555555")
  }
  rect(xr[1],0.25,xr[2],0.78,col="#FCE9E9",border="#C95A5A")
  text(mean(xr),0.52,"No proposed axis-specific programme passes global FDR",font=2,cex=0.66,col="#9E2F2F")
  title("Axis-specific localisation remains null",adj=0.14,cex.main=0.88,font.main=2,line=0.55)
  panel_label("E")

}

base<-file.path(outdir,"Figure6_single_cell_localisation")
width_mm<-180; height_mm<-150; width_in<-width_mm/25.4; height_in<-height_mm/25.4
if(requireNamespace("svglite",quietly=TRUE)) svglite::svglite(paste0(base,".svg"),width=width_in,height=height_in,system_fonts=list(Arial="Arial"),pointsize=8) else svg(paste0(base,".svg"),width=width_in,height=height_in,family="Arial",pointsize=8)
par(oma=c(0.2,0.2,0.4,0.2)); draw_figure(); dev.off()
cairo_pdf(paste0(base,".pdf"),width=width_in,height=height_in,family="Arial",pointsize=8); par(oma=c(0.2,0.2,0.4,0.2)); draw_figure(); dev.off()
tiff(paste0(base,".tiff"),width=width_mm,height=height_mm,units="mm",res=600,compression="lzw",pointsize=8,family="Arial"); par(oma=c(0.2,0.2,0.4,0.2)); draw_figure(); dev.off()
png(paste0(base,".png"),width=width_mm,height=height_mm,units="mm",res=600,pointsize=8,family="Arial"); par(oma=c(0.2,0.2,0.4,0.2)); draw_figure(); dev.off()
cat("Figure 6 exports written to",outdir,"\n"); cat("Donor-paired source rows:",nrow(paired),"; exact sign-flip rows:",nrow(exact_sign),"\n")
