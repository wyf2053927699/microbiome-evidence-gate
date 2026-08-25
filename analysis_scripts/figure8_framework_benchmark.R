options(stringsAsFactors=FALSE)
library(ggplot2); library(patchwork); library(scales)
theme_set(theme_classic(base_size=9.5,base_family="Times New Roman")+
  theme(axis.text=element_text(size=8.4,colour="#30343B"),axis.title=element_text(size=9.4),
        plot.title=element_text(size=10.3,face="bold",margin=margin(b=4)),strip.text=element_text(size=9.1,face="bold"),
        legend.text=element_text(size=8.1),legend.title=element_text(size=8.5,face="bold"),
        axis.line=element_line(linewidth=.42,colour="#30343B"),axis.ticks=element_line(linewidth=.38),plot.margin=margin(7,8,7,8)))
navy<-"#1F4E79"; teal<-"#2A9D8F"; gold<-"#D99B2B"; coral<-"#C95C54"; grey<-"#8A94A1"; ink<-"#20242A"
root <- Sys.getenv("MICROBIOME_GATE_ROOT", unset = getwd())
b <- file.path(root, "derived_data", "framework_benchmark_v2")
states <- read.csv(file.path(b,"all_binary_patterns.csv"),check.names=FALSE)
trans <- read.csv(file.path(b,"bidirectional_transition_summary.csv"),check.names=FALSE); names(trans)[1] <- "workflow"
err <- read.csv(file.path(b,"coding_error_summary.csv"),check.names=FALSE); names(err)[1] <- "per_dimension_coding_error"
gate <- read.csv(file.path(root, "derived_data", "crc_axis_v1_1", "crc_axis_updated_gate_matrix.csv"),check.names=FALSE)
rules <- c("any_positive","majority_dimensions","additive_score_ge_0_5","noncompensatory_all_complete")
rule_lab <- c(any_positive="Any positive",majority_dimensions="Majority",additive_score_ge_0_5="Additive >=0.5",noncompensatory_all_complete="All complete")
rule_col <- c(any_positive=grey,majority_dimensions=gold,additive_score_ge_0_5=teal,noncompensatory_all_complete=navy)

# A gate matrix
dims <- c("microbial_reproducibility","human_metabolite_observed","microbe_metabolite_link","metabolite_host_perturbation","host_reproduced_two_datasets","single_cell_localisation","major_contradiction_resolved","robust_to_analytic_choices")
dim_lab <- c("Microbial\nreproducibility","Measured\nmetabolite","Microbe-\nmetabolite link","Metabolite-host\nperturbation","Host replication\n>=2 datasets","Single-cell\nlocalisation","Contradiction\nresolved","Analytic\nrobustness")
gm <- do.call(rbind,lapply(seq_len(nrow(gate)),function(i)data.frame(case=gate$candidate_family[i],dimension=dims,value=as.numeric(gate[i,dims]),claim=gate$claim_class[i])))
gm$case <- factor(gm$case,levels=rev(gate$candidate_family)); gm$dimension <- factor(gm$dimension,levels=dims,labels=dim_lab)
pA <- ggplot(gm,aes(dimension,case,fill=value)) + geom_tile(colour="white",linewidth=1.2) +
  geom_text(aes(label=ifelse(value==1,"PASS",ifelse(value==.5,"PARTIAL","STOP"))),size=3.0,family="Times New Roman",fontface="bold",colour=ink) +
  scale_fill_gradientn(colours=c("#F3D7D4","#F7E9CE","#D9ECE8"),values=c(0,.5,1),limits=c(0,1),guide="none") +
  labs(x=NULL,y=NULL,title="Cross-disease evidence gates preserve bounded claim classes") +
  theme(axis.text.x=element_text(size=7.7,angle=0,lineheight=.9),axis.text.y=element_text(size=8.1),axis.line=element_blank(),axis.ticks=element_blank(),plot.margin=margin(7,9,7,9))

# B exhaustive states
prom <- data.frame(workflow=rules,promotion_rate=sapply(states[rules],mean))
prom$workflow <- factor(prom$workflow,levels=rules,labels=rule_lab[rules])
pB <- ggplot(prom,aes(promotion_rate,workflow)) + geom_segment(aes(x=0,xend=promotion_rate,yend=workflow),linewidth=1.1,colour="#D4D8DD") +
  geom_point(size=3.2,aes(colour=workflow)) + geom_text(aes(label=percent(promotion_rate,accuracy=.1),x=promotion_rate+.035),hjust=0,size=3.1,family="Times New Roman") +
  scale_colour_manual(values=unname(rule_col[rules]),guide="none") + coord_cartesian(xlim=c(0,1.12),clip="off") +
  labs(x="Promotion frequency across 256 states",y=NULL,title="Exhaustive binary state space")

# C one-bit transitions
trans$workflow <- factor(trans$workflow,levels=rules,labels=rule_lab[rules]); trans$direction <- factor(trans$direction,levels=c("0->1","1->0"),labels=c("Evidence added","Evidence removed"))
pC <- ggplot(trans,aes(decision_change_rate,workflow,colour=direction)) +
  geom_point(size=2.8,position=position_dodge(width=.42)) + geom_segment(aes(x=0,xend=decision_change_rate,yend=workflow),position=position_dodge(width=.42),alpha=.45,linewidth=.8) +
  scale_colour_manual(values=c("Evidence added"=teal,"Evidence removed"=coral)) + scale_x_continuous(labels=percent_format(accuracy=1),limits=c(0,.32)) +
  labs(x="Decision-change rate",y=NULL,colour=NULL,title="2,048 bidirectional one-bit transitions") + theme(legend.position="bottom")

# D coding error trade-off
err$workflow <- factor(err$workflow,levels=rules,labels=rule_lab[rules])
e1 <- err[,c("per_dimension_coding_error","workflow","mean_false_promotion_near_boundary")]; names(e1)[3]<-"rate"; e1$outcome<-"False promotion near boundary"
e2 <- err[,c("per_dimension_coding_error","workflow","false_stop_complete_case")]; names(e2)[3]<-"rate"; e2$outcome<-"False stop of complete case"
ee<-rbind(e1,e2)
pD <- ggplot(ee,aes(per_dimension_coding_error,rate,colour=workflow,linetype=workflow)) + geom_line(linewidth=.85) + geom_point(size=1.9) +
  facet_wrap(~outcome,nrow=1) + scale_colour_manual(values=unname(rule_col[rules])) + scale_linetype_manual(values=c(2,3,4,1)) +
  scale_x_continuous(labels=percent_format(accuracy=1)) + scale_y_continuous(labels=percent_format(accuracy=1),limits=c(0,1)) +
  labs(x="Per-dimension coding-error probability",y="Simulated error rate",colour=NULL,linetype=NULL,title="Boundary-focused coding-error stress test") +
  theme(legend.position="bottom",legend.key.width=unit(7,"mm"))

# E operational interpretation
box <- data.frame(x=c(1,3.2,5.4),y=1,label=c("HCC\nSTOP / hypothesis only","CRC\ntriangulated candidate","Benchmark\nspecification conformance"),fill=c("#F3D7D4","#D9ECE8","#DCEAF5"))
pE <- ggplot() + geom_label(data=box,aes(x=x,y=y,label=label,fill=fill),size=3.6,family="Times New Roman",fontface="bold",linewidth=.4,label.padding=unit(.32,"lines"),show.legend=FALSE) +
  annotate("text",x=3.2,y=.28,label="The framework calibrates claim class; it does not certify causal truth or clinical accuracy.",size=3.4,family="Times New Roman",fontface="bold",colour=ink) +
  scale_fill_identity()+coord_cartesian(xlim=c(0,6.4),ylim=c(.05,1.55),clip="off")+theme_void(base_family="Times New Roman")+ggtitle("Operational interpretation")+theme(plot.title=element_text(size=10.3,face="bold"))

fig <- pA / (pB|pC) / pD / pE + plot_layout(heights=c(1.25,1,1.08,.63)) +
  plot_annotation(tag_levels="A",theme=theme(plot.tag=element_text(family="Times New Roman",face="bold",size=13,colour=ink),plot.tag.position=c(.006,.992)))
outdir <- Sys.getenv("MICROBIOME_FIGURE_OUT", unset = file.path(root, "figures")); dir.create(outdir,recursive=TRUE,showWarnings=FALSE); outbase<-file.path(outdir,"Figure_8")
svglite::svglite(paste0(outbase,".svg"),width=183/25.4,height=185/25.4);print(fig);dev.off()
grDevices::cairo_pdf(paste0(outbase,".pdf"),width=183/25.4,height=185/25.4,family="Times New Roman");print(fig);dev.off()
ragg::agg_tiff(paste0(outbase,".tiff"),width=183/25.4,height=185/25.4,units="in",res=600,compression="lzw");print(fig);dev.off()
ragg::agg_png(paste0(outbase,".png"),width=183/25.4,height=185/25.4,units="in",res=300);print(fig);dev.off()
