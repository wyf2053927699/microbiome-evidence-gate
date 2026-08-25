options(stringsAsFactors = FALSE)
library(ggplot2)
library(patchwork)
library(scales)

theme_set(theme_classic(base_size = 9.4, base_family = "Times New Roman") +
  theme(axis.line = element_line(linewidth = 0.42, colour = "#30343B"),
        axis.ticks = element_line(linewidth = 0.38, colour = "#30343B"),
        axis.text = element_text(size = 8.4, colour = "#30343B"),
        axis.title = element_text(size = 9.4, colour = "#20242A"),
        legend.title = element_text(size = 8.7, face = "bold"),
        legend.text = element_text(size = 8.1),
        strip.text = element_text(size = 9.1, face = "bold"),
        plot.title = element_text(size = 10.3, face = "bold", margin = margin(b = 4)),
        panel.grid = element_blank(), plot.margin = margin(7, 8, 7, 8)))

navy <- "#1F4E79"; teal <- "#2A9D8F"; gold <- "#D99B2B"; coral <- "#C95C54"
grey <- "#8A94A1"; pale <- "#E8EEF3"; ink <- "#20242A"; green <- "#3A7D44"
root <- Sys.getenv("MICROBIOME_GATE_ROOT", unset = getwd())
base <- file.path(root, "derived_data", "crc_axis_v1_1")
lf <- read.csv(file.path(base, "species_metabolite_leakage_free_results.csv"), check.names = FALSE)
rb <- read.csv(file.path(base, "robustness", "robustness_consensus.csv"), check.names = FALSE)
ji <- read.csv(file.path(base, "independence_audit", "joint_species_independence_results.csv"), check.names = FALSE)
pca <- read.csv(file.path(base, "host_program_GSE241076", "sample_pca_scores.csv"), check.names = FALSE)
cnt <- read.delim(gzfile(file.path(root, "01_raw_data", "crc_axis_external", "GSE241076", "GSE241076_counts.txt.gz")), check.names = FALSE)

clean_met <- function(x) {
  x <- sub("^C[0-9]+_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x
}
clean_sp <- function(x) sub("Peptostreptococcus ", "P. ", sub("Fusobacterium ", "F. ", sub("Gemella ", "G. ", sub("Parvimonas ", "Parvimonas ", x))))

# A: design schematic
boxes <- data.frame(x = c(0.8, 2.5, 4.2, 5.9), y = 1,
                    label = c("347 paired\nparticipants", "Discovery\nn = 207", "Held-out validation\nn = 140", "Attribution +\nhost perturbation"),
                    fill = c(pale, "#DCEAF5", "#DDF1EC", "#F7E9CE"))
pA <- ggplot() +
  geom_segment(data = data.frame(x = c(1.45,3.15,4.85), xend = c(1.85,3.55,5.25), y = 1, yend = 1),
               aes(x=x,xend=xend,y=y,yend=yend), arrow=arrow(length=unit(2.2,"mm")), linewidth=.7, colour=grey) +
  geom_label(data=boxes,aes(x=x,y=y,label=label,fill=fill),size=3.15,family="Times New Roman",fontface="bold",linewidth=.35,label.padding=unit(.23,"lines"),lineheight=.95,show.legend=FALSE) +
  annotate("text",x=3.35,y=.30,label="4 locked oral taxa | 252 metabolites | 1,008 tests",size=3.1,family="Times New Roman",colour=ink) +
  scale_fill_identity() + coord_cartesian(xlim=c(0,6.7),ylim=c(.05,1.65),clip="off") + theme_void(base_family="Times New Roman") +
  ggtitle("Leakage-free paired-cohort design") + theme(plot.title=element_text(size=10.3,face="bold"))

# B: discovery versus validation effect map
rep <- lf[lf$heldout_reproduced %in% c(TRUE,"True","TRUE"),]
key <- paste(ji$species[ji$crc_preferential_exploratory %in% c(TRUE,"True","TRUE")], ji$metabolite[ji$crc_preferential_exploratory %in% c(TRUE,"True","TRUE")])
rep$key <- paste(rep$species, rep$metabolite)
rep$priority <- ifelse(rep$key %in% key, "Joint-taxon CRC-preferential", "Held-out reproduced")
pB <- ggplot(rep,aes(partial_spearman_discovery,partial_spearman_validation,colour=priority)) +
  geom_hline(yintercept=0,linewidth=.35,colour="#D2D6DB") + geom_vline(xintercept=0,linewidth=.35,colour="#D2D6DB") +
  geom_abline(slope=1,intercept=0,linetype=2,linewidth=.45,colour=grey) +
  geom_point(alpha=.84,size=2.15) +
  scale_colour_manual(values=c("Held-out reproduced"=grey,"Joint-taxon CRC-preferential"=coral)) +
  labs(x="Discovery partial Spearman",y="Validation partial Spearman",colour=NULL,title="Held-out effect concordance") +
  theme(legend.position="bottom",legend.key.width=unit(4,"mm"),legend.margin=margin(t=-2))

# C: six joint-taxon candidates
top <- ji[ji$crc_preferential_exploratory %in% c(TRUE,"True","TRUE"),]
top$label <- paste(clean_sp(top$species), clean_met(top$metabolite), sep=" | ")
top$lo <- top$beta_joint_crc - 1.96*top$se_joint_crc; top$hi <- top$beta_joint_crc + 1.96*top$se_joint_crc
top$label <- factor(top$label, levels=rev(top$label[order(top$beta_joint_crc)]))
pC <- ggplot(top,aes(beta_joint_crc,label,colour=species)) +
  geom_vline(xintercept=0,colour="#BFC5CC",linewidth=.45) +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.15,linewidth=.65) + geom_point(size=2.8) +
  scale_colour_manual(values=c("Peptostreptococcus stomatis"=navy,"Fusobacterium nucleatum"=gold),guide="none") +
  labs(x="Joint-taxon CRC coefficient (95% CI)",y=NULL,title="Six attribution-prioritised candidates") +
  theme(axis.text.y=element_text(size=8.1))

# D: glycocholate stratified joint estimates
g <- ji[ji$species=="Peptostreptococcus stomatis" & grepl("Glycocholate",ji$metabolite),]
gd <- data.frame(group=c("All paired participants","CRC only","Healthy only"),
                 beta=c(g$beta_joint_all,g$beta_joint_crc,g$beta_joint_healthy),
                 se=c(g$se_joint_all,g$se_joint_crc,g$se_joint_healthy),
                 q=c(g$q_joint_all,g$q_joint_crc,g$q_joint_healthy))
gd$lo <- gd$beta-1.96*gd$se; gd$hi <- gd$beta+1.96*gd$se
gd$group <- factor(gd$group,levels=rev(gd$group))
pD <- ggplot(gd,aes(beta,group)) + geom_vline(xintercept=0,colour="#BFC5CC",linewidth=.45) +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.13,linewidth=.7,colour=navy) + geom_point(size=3,colour=navy) +
  geom_text(aes(x=hi+.055,label=paste0("q=",formatC(q,format="f",digits=3))),hjust=0,size=3.0,family="Times New Roman") +
  coord_cartesian(xlim=c(min(gd$lo)-.03,max(gd$hi)+.28),clip="off") +
  labs(x="Joint-taxon coefficient (95% CI)",y=NULL,title="P. stomatis-associated glycocholate")

# E: external perturbation PCA
pca$condition <- factor(ifelse(tolower(pca$condition)=="gca","GCA","Control"),levels=c("Control","GCA"))
pE <- ggplot(pca,aes(PC1,PC2,colour=condition,shape=condition)) +
  geom_point(size=3.5,stroke=.7) +
  scale_colour_manual(values=c("Control"=grey,"GCA"=coral)) + scale_shape_manual(values=c(16,17)) +
  labs(x="PC1 (64.9% variance)",y="PC2 (20.8% variance)",colour=NULL,shape=NULL,title="GCA perturbation separates CRC cells") +
  theme(legend.position="bottom")

# F: mechanism-gene expression
rownames(cnt) <- cnt$gene_id; cnt$gene_id <- NULL
samples <- c("con_1","con_2","con_3","fa_1","fa_2","fa_3"); lib <- colSums(cnt[,samples])
genes <- c("SOX14","ZDHHC9","CD274")
ex <- do.call(rbind,lapply(genes,function(z) data.frame(gene=z,sample=samples,logCPM=log2(as.numeric(cnt[z,samples])/lib*1e6+.5),condition=rep(c("Control","GCA"),each=3))))
ex$gene <- factor(ex$gene,levels=genes)
pF <- ggplot(ex,aes(condition,logCPM,colour=condition)) +
  geom_point(position=position_jitter(width=.055,height=0),size=2.4,alpha=.9) +
  stat_summary(fun=mean,geom="crossbar",width=.46,linewidth=.62,colour=ink) +
  facet_wrap(~gene,scales="free_y",nrow=1) + scale_colour_manual(values=c("Control"=grey,"GCA"=coral),guide="none") +
  labs(x=NULL,y="log2(CPM + 0.5)",title="Locked host-response genes increase after GCA") +
  theme(axis.text.x=element_text(angle=0))

# G: evidence chain and claim boundary
chain <- data.frame(x=c(.9,3.4,5.9),y=1,label=c("P. stomatis\nabundance","Faecal\nglycocholate","SOX14-ZDHHC9-CD274\nhost response"),fill=c("#DCEAF5","#F7E9CE","#F5DDDA"))
pG <- ggplot() +
  geom_segment(aes(x=1.65,xend=2.65,y=1,yend=1),arrow=arrow(length=unit(2.4,"mm")),linewidth=.8,colour=navy) +
  geom_segment(aes(x=4.15,xend=5.15,y=1,yend=1),arrow=arrow(length=unit(2.4,"mm")),linewidth=.8,colour=coral) +
  geom_label(data=chain,aes(x=x,y=y,label=label,fill=fill),size=3.5,family="Times New Roman",fontface="bold",linewidth=.4,label.padding=unit(.3,"lines"),show.legend=FALSE) +
  annotate("text",x=2.15,y=1.36,label="paired-cohort association\n+ held-out replication",size=3.05,family="Times New Roman",colour=navy) +
  annotate("text",x=4.65,y=1.36,label="external GCA perturbation\n+ transcriptomic response",size=3.05,family="Times New Roman",colour=coral) +
  annotate("label",x=3.4,y=.32,label="TRIANGULATED CANDIDATE AXIS - not a same-participant causal mechanism",size=3.2,family="Times New Roman",fontface="bold",fill="#F2F3F5",linewidth=.3) +
  scale_fill_identity()+coord_cartesian(xlim=c(0,6.8),ylim=c(.08,1.65),clip="off")+theme_void(base_family="Times New Roman")+
  ggtitle("Evidence integration preserves the boundary of each link")+theme(plot.title=element_text(size=10.3,face="bold"))

fig <- (pA | pB) / (pC | pD) / (pE | pF) / pG +
  plot_layout(heights=c(.82,1.17,1.12,.82)) +
  plot_annotation(tag_levels="A",theme=theme(plot.tag=element_text(family="Times New Roman",face="bold",size=13,colour=ink),plot.tag.position=c(.006,.992)))

outdir <- Sys.getenv("MICROBIOME_FIGURE_OUT", unset = file.path(root, "figures")); dir.create(outdir, recursive=TRUE, showWarnings=FALSE); outbase <- file.path(outdir,"Figure_7")
svglite::svglite(paste0(outbase,".svg"),width=183/25.4,height=220/25.4); print(fig); dev.off()
grDevices::cairo_pdf(paste0(outbase,".pdf"),width=183/25.4,height=220/25.4,family="Times New Roman"); print(fig); dev.off()
ragg::agg_tiff(paste0(outbase,".tiff"),width=183/25.4,height=220/25.4,units="in",res=600,compression="lzw"); print(fig); dev.off()
ragg::agg_png(paste0(outbase,".png"),width=183/25.4,height=220/25.4,units="in",res=300); print(fig); dev.off()
