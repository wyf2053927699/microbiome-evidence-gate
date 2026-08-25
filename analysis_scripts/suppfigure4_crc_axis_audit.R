options(stringsAsFactors=FALSE)
library(ggplot2);library(patchwork);library(scales)
theme_set(theme_classic(base_size=9.4,base_family="Times New Roman")+theme(axis.text=element_text(size=8.2,colour="#30343B"),axis.title=element_text(size=9.2),plot.title=element_text(size=10.1,face="bold",margin=margin(b=4)),strip.text=element_text(size=8.8,face="bold"),legend.text=element_text(size=8),legend.title=element_text(size=8.4,face="bold"),axis.line=element_line(linewidth=.4),axis.ticks=element_line(linewidth=.36),plot.margin=margin(7,8,7,8)))
navy<-"#1F4E79";teal<-"#2A9D8F";gold<-"#D99B2B";coral<-"#C95C54";grey<-"#8A94A1";ink<-"#20242A"
root <- Sys.getenv("MICROBIOME_GATE_ROOT", unset = getwd())
base<-file.path(root, "derived_data", "crc_axis_v1_1")
rb<-read.csv(file.path(base,"robustness","robustness_consensus.csv"),check.names=FALSE)
ji<-read.csv(file.path(base,"independence_audit","joint_species_independence_results.csv"),check.names=FALSE)
tc<-read.csv(file.path(base,"independence_audit","locked_species_rank_correlation.csv"),check.names=FALSE,row.names=1)
sc<-read.csv(file.path(base,"host_program_GSE241076","sample_logCPM_correlation.csv"),check.names=FALSE,row.names=1)
qa<-jsonlite::fromJSON(file.path(base,"host_program_GSE241076","quality_audit_summary.json"))
clean_met<-function(x){x<-sub("^C[0-9]+_","",x);gsub("_"," ",x,fixed=TRUE)}
clean_sp<-function(x)sub("Peptostreptococcus ","P. ",sub("Fusobacterium ","F. ",sub("Gemella ","G. ",sub("Parvimonas ","Parvimonas ",x))))

# A candidate attrition
at<-data.frame(stage=factor(c("All tests","Discovery","Held-out","Robustness","Joint attribution"),levels=rev(c("All tests","Discovery","Held-out","Robustness","Joint attribution"))),n=c(1008,223,96,72,6),class=c("screen","screen","screen","audit","priority"))
pA<-ggplot(at,aes(n,stage,colour=class))+geom_segment(aes(x=0,xend=n,yend=stage),linewidth=1.15,colour="#D7DBE0")+geom_point(size=3.5)+geom_text(aes(label=comma(n),x=n+35),hjust=0,size=3.2,family="Times New Roman",fontface="bold")+scale_colour_manual(values=c(screen=grey,audit=gold,priority=coral),guide="none")+coord_cartesian(xlim=c(0,1150),clip="off")+labs(x="Candidate pairs",y=NULL,title="Candidate attrition is explicit")

# B alternative profiler
mb<-rb[rb$motu_exact %in% c(TRUE,"True","TRUE"),]
pB<-ggplot(mb,aes(partial_spearman_validation,motu_r,colour=species))+geom_hline(yintercept=0,colour="#D3D7DC",linewidth=.35)+geom_vline(xintercept=0,colour="#D3D7DC",linewidth=.35)+geom_abline(slope=1,intercept=0,linetype=2,colour=grey,linewidth=.45)+geom_point(size=2.15,alpha=.82)+scale_colour_manual(values=c("Peptostreptococcus stomatis"=navy,"Gemella morbillorum"=teal,"Fusobacterium nucleatum"=gold),labels=clean_sp)+labs(x="MetaPhlAn validation effect",y="mOTU sensitivity effect",colour=NULL,title="Alternative taxonomic profiler")+theme(legend.position="bottom")

# C locked taxon correlation matrix
td<-as.data.frame(as.table(as.matrix(tc)));names(td)<-c("row","col","r");td$row<-factor(clean_sp(td$row),levels=rev(clean_sp(rownames(tc))));td$col<-factor(clean_sp(td$col),levels=clean_sp(colnames(tc)))
pC<-ggplot(td,aes(col,row,fill=r))+geom_tile(colour="white",linewidth=.8)+geom_text(aes(label=sprintf("%.2f",r)),size=3.0,family="Times New Roman")+scale_fill_gradient2(low="#4C78A8",mid="white",high="#C95C54",limits=c(-1,1),name="Rank r")+labs(x=NULL,y=NULL,title="Locked-taxon co-occurrence")+theme(axis.text.x=element_text(angle=30,hjust=1),axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")

# D top-six multi-analysis effect matrix
top<-ji[ji$crc_preferential_exploratory %in% c(TRUE,"True","TRUE"),]
em<-do.call(rbind,lapply(seq_len(nrow(top)),function(i)data.frame(candidate=paste(clean_sp(top$species[i]),clean_met(top$metabolite[i]),sep=" | "),analysis=c("Discovery","Held-out","mOTU","CRC only","Healthy only","Presence/absence"),effect=as.numeric(top[i,c("primary_discovery_r","primary_validation_r","motu_r","crc_univariable_r","beta_joint_healthy","presence_r")]))))
em$candidate<-factor(em$candidate,levels=rev(unique(em$candidate)));em$analysis<-factor(em$analysis,levels=c("Discovery","Held-out","mOTU","CRC only","Healthy only","Presence/absence"))
pD<-ggplot(em,aes(analysis,candidate,fill=effect))+geom_tile(colour="white",linewidth=.65)+geom_text(aes(label=sprintf("%.2f",effect)),size=2.75,family="Times New Roman")+scale_fill_gradient2(low="#4C78A8",mid="white",high="#C95C54",limits=c(-.65,.65),name="Effect")+labs(x=NULL,y=NULL,title="Direction and attribution audit")+theme(axis.text.x=element_text(angle=28,hjust=1),axis.text.y=element_text(size=7.6),axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")

# E host sample correlations
sd<-as.data.frame(as.table(as.matrix(sc)));names(sd)<-c("row","col","r");sd$row<-factor(sd$row,levels=rev(rownames(sc)));sd$col<-factor(sd$col,levels=colnames(sc))
pE<-ggplot(sd,aes(col,row,fill=r))+geom_tile(colour="white",linewidth=.75)+geom_text(aes(label=sprintf("%.3f",r)),size=2.65,family="Times New Roman")+scale_fill_gradient(low="#F2F4F6",high=navy,limits=c(.95,1),name="r")+labs(x=NULL,y=NULL,title="GSE241076 sample correlation")+theme(axis.text.x=element_text(angle=30,hjust=1),axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")

# F counts-versus-FPKM concordance
cf<-qa$counts_FPKM_mechanism_concordance
pF<-ggplot(cf,aes(counts_logCPM_log2FC,FPKM_log2FC,label=gene))+geom_abline(slope=1,intercept=0,linetype=2,colour=grey)+geom_point(size=3.2,colour=coral)+geom_text(nudge_y=.35,size=3.2,family="Times New Roman",fontface="bold")+coord_equal(xlim=c(4.5,9.5),ylim=c(4.5,9.5))+labs(x="Counts-derived log2FC",y="FPKM-derived log2FC",title="Expression-format concordance")

# G library-size audit
ls<-data.frame(sample=names(qa$library_sizes),reads=as.numeric(qa$library_sizes),condition=rep(c("Control","GCA"),each=3));ls$sample<-factor(ls$sample,levels=ls$sample)
pG<-ggplot(ls,aes(sample,reads/1e6,colour=condition))+geom_segment(aes(xend=sample,y=0,yend=reads/1e6),linewidth=1.1,colour="#D7DBE0")+geom_point(size=3.1)+scale_colour_manual(values=c(Control=grey,GCA=coral))+coord_cartesian(ylim=c(0,22))+labs(x=NULL,y="Library size (million reads)",colour=NULL,title="Comparable sequencing depth")+theme(legend.position="bottom")

fig<-(pA|pB)/(pC|pD)/(pE|pF)/pG+plot_layout(heights=c(.9,1.18,1.08,.78))+plot_annotation(tag_levels="A",theme=theme(plot.tag=element_text(family="Times New Roman",face="bold",size=13,colour=ink),plot.tag.position=c(.006,.992)))
outdir<-Sys.getenv("MICROBIOME_FIGURE_OUT", unset = file.path(root, "figures"));dir.create(outdir,recursive=TRUE,showWarnings=FALSE);outbase<-file.path(outdir,"Supplementary_Figure_4")
svglite::svglite(paste0(outbase,".svg"),width=183/25.4,height=218/25.4);print(fig);dev.off()
grDevices::cairo_pdf(paste0(outbase,".pdf"),width=183/25.4,height=218/25.4,family="Times New Roman");print(fig);dev.off()
ragg::agg_tiff(paste0(outbase,".tiff"),width=183/25.4,height=218/25.4,units="in",res=600,compression="lzw");print(fig);dev.off()
ragg::agg_png(paste0(outbase,".png"),width=183/25.4,height=218/25.4,units="in",res=300);print(fig);dev.off()
