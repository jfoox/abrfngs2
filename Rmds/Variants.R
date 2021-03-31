# Dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(cowplot)
library(reshape2)
options(scipen=999)

theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_jon = theme_linedraw2 + 
  theme(legend.position="bottom", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=12), 
        axis.title.y=element_text(vjust=1, size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5), 
        strip.text = element_text(size=11), 
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0, 0.5, 0), "cm"))

platform.colors = c("BGISEQ500"="#fd8d3c", "HiSeq2500"="#6baed6", "HiSeq4000"="#3182bd", "HiSeqX10"="#08519c", "MGISEQ2000"="#e6550d", "NovaSeq2x150"="#f768a1", "NovaSeq2x250"="#ae017e", "PacBio" = "#74c476", "Nanopore"="#006d2c")
context.labels = c("Alu"="Alu", "L1"="L1", "L2"="L2",
                   "lowcomplexity"="Low Complexity",
                   "LTR"="LTR", "satellite"="Satellite", "simplerepeat"="Simple Repeat")


# Caller Comparison (precision vs sensitivity)

## Set up data matrix
df_sum = read.csv("tables/rtg_summary_v2.csv")
df_sum = separate(df_sum, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_sum$platform = factor(df_sum$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000", "PacBio", "Nanopore"))
df_sum$caller = factor(df_sum$caller, levels = c("dv", "gatk", "sen", "strelka2", "clair"))
df_sum

## Short-read values
gg_sumCaller = ggplot(subset(df_sum, caller != "clair"), aes(x=sensitivity, y=precision)) +
  geom_point(size=3, shape=21, aes(fill=caller)) +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.5) + geom_vline(xintercept=1, linetype="dashed", alpha=0.5) +
  theme_jon + theme(legend.position="top", legend.title=element_blank()) +
  scale_fill_discrete(name = "", labels = c("dv"="Deep Variant", "gatk"="GATK", "sen"="Sentieon", "strelka2"="strelka2")) +
  xlab("\nSensitivity") + ylab("Precision\n") +
  facet_wrap(~platform, nrow=1) +
  guides(fill = guide_legend(nrow = 1))

## Long-read values
gg_sumCaller2 = ggplot(subset(df_sum, caller=="clair"), aes(x=sensitivity, y=precision)) +
  geom_point(size=3, shape=21, aes(fill=platform)) +
  geom_hline(yintercept=1, linetype="dashed", alpha=0.5) + geom_vline(xintercept=1, linetype="dashed", alpha=0.5) +
  theme_jon + theme(legend.position="top", legend.title=element_blank()) +
  scale_fill_manual(values = platform.colors) +
  xlab("\n") + ylab("") +
  facet_wrap(~caller, nrow=1, labeller = labeller(caller = c("clair"="Long Reads"))) +
  guides(fill = guide_legend(nrow = 1))

## Combine
gg_caller = plot_grid(gg_sumCaller, gg_sumCaller2, nrow=1, rel_widths=c(10,2.75))

## Export
pdf(file="figures/Fig4a_AlgoComp.pdf", width=12, height=3.5)
gg_caller
dev.off()

# Inter- and intra-lab reproducibility (using F1)
gg_sumLabRep = ggplot(subset(df_sum, caller=="dv" | caller == "clair"), aes(x=genome, y=fmeasure)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform, shape=genome), show.legend=F, size=2, groupOnX = T) +
  geom_hline(yintercept=1) +
  theme_jon + theme(strip.text = element_text(size=9), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24)) +
  xlab("") + ylab("F1 Measure") +
  facet_wrap(~platform, scales="free_x", nrow=1)
gg_sumLabRep

# RTG Outputs /// ROC plots

## Set up matrix
rtg_snps = read.csv("tables/rtg_ROC_SNP.csv")
rtg_snps$VarType = "SNPs"
rtg_indels = read.csv("tables/rtg_ROC_INDEL.csv")
rtg_indels$VarType = "INDELs"
df_ROC = rbind(rtg_snps, rtg_indels)
df_ROC$VarType = factor(df_ROC$VarType, levels=c("SNPs", "INDELs"))
df_ROC = subset(df_ROC, df_ROC$score %in% c(0,20,40,60,80,99))
df_ROC = separate(df_ROC, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_ROC$platform = factor(df_ROC$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000"))

## Create ROC plot
gg_ROC_all = ggplot(subset(df_ROC), aes(x=sensitivity, y=precision, color=score)) +
  geom_point(size=2, shape=21) +
  geom_line(aes(group=sample)) +
  #geom_text(aes(label=score), hjust=-0.5, vjust=0) +
  theme_cem + 
  theme(legend.title=element_text(size=10),
      legend.text=element_text(size=10),
      axis.text.y = element_text(size=9),
      axis.text.x = element_text(size=9, angle=90, hjust=1, vjust=0.5),
      axis.title.y = element_text(size=10),
      axis.title.x = element_text(size=10),
      strip.text.x = element_text(size=10),
      strip.text.y = element_text(size=10),
      plot.margin = unit(c(1, 1, 0, 1), "cm"),
      legend.position = "bottom") +
  facet_grid(platform~VarType+caller) +
  xlab("\nSensitivity") + ylab("Precision\n")


# Variant Metrics in UCSC RepeatMasker Context

## Set up data matrix
df_ucsc = read.csv("tables/rtg_ucsc.csv")
df_ucsc = separate(df_ucsc, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_ucsc$platform = factor(df_ucsc$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000", "PacBio", "Nanopore"))

## Plot
gg_ucsc = ggplot(subset(df_ucsc, platform != "PromethION"), aes(x=sensitivity, y=precision)) +
  geom_point(size=3, aes(fill=platform, shape=genome), alpha=0.85) + 
  geom_hline(yintercept=1, linetype="dashed", alpha=0.5) + geom_vline(xintercept=1, linetype="dashed", alpha=0.5) +
  theme_jon + theme(legend.position="top") + #theme(legend.title=element_blank(), legend.position = c(0.95,0.07), legend.justification = c(0.95,0.07)) +
  scale_fill_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24)) +
  xlab("\nSensitivity") + ylab("Precision\n")  +
  facet_wrap(~context, nrow=1, labeller = labeller(context = context.labels)) +
  guides(fill=guide_legend(nrow=2, override.aes=list(shape=21))) 

## Export
pdf(file="figures/Fig4b_UCSC.pdf", width=12, height=4)
gg_ucsc
dev.off()

# Altuna Plots (presence/absence of TPs)

## Create matrix based on each platform's output
altuna_plots = list()
plotnum = 0
for(varType in c("SNPs", "nonSNPs")){
  for(context in c("Alu", "L1", "L2", "lowcomplexity", "LTR", "satellite", "simplerepeat")){
    for(genome in c("Son")){
      plotnum = plotnum + 1
      df_altuna = read.csv(paste0("rtg/outs/altuna/combined_",genome,"_",context,"_",varType,".matrix.tsv"), sep="\t")
      df_altuna[2:ncol(df_altuna)] = lapply(df_altuna[2:ncol(df_altuna)], function(x) ifelse(x == ".", 0, 1) )
      df_altuna = df_altuna %>% filter_at(vars(colnames(df_altuna[2:ncol(df_altuna)])), any_vars(. %in% c(1)))
      colnames(df_altuna) = gsub("PromethION", "Nanopore", colnames(df_altuna))
      for(platform in c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "NovaSeq2x150", "NovaSeq2x250", "PacBio", "Nanopore")){
        if(any(grepl(platform, colnames(df_altuna)))){
          idx = grepl(platform, colnames(df_altuna))
          if(sum(idx) > 1){
            df_altuna[, platform] = floor(rowMeans(df_altuna[, grepl(platform, colnames(df_altuna))]))  
          } else {
            df_altuna[, platform] = df_altuna[, idx]
          }
        }
      }
      df_altuna = df_altuna[, !grepl("_", colnames(df_altuna))]
      df_altuna2 = df_altuna[, 2:ncol(df_altuna)]
      df_altuna2 = df_altuna2[do.call(order, as.list(df_altuna2)), ]
      df_altuna2$id = 1:nrow(df_altuna2)
      df_altuna3 = melt(df_altuna2, "id")
      df_altuna3$value = factor(df_altuna3$value, levels=c(1, 0))
      
      if(plotnum %in% c(1,8)){
        curPlot = ggplot(df_altuna3, aes(x=rev(id), y=variable, fill=value, color=value)) + 
          geom_tile(show.legend=F) + scale_color_manual(values=c("#fed976", "#0570b0")) + 
          scale_fill_manual(values=c("#fed976", "#0570b0")) +  
          xlab("") + ylab("") + theme_jon + theme(plot.margin = unit(c(0.1, 0, 0.1, 0), "cm")) + ggtitle(context)
      } else {
        curPlot = ggplot(df_altuna3, aes(x=rev(id), y=variable, fill=value, color=value)) + 
          geom_tile(show.legend=F) + scale_color_manual(values=c("#fed976", "#0570b0")) + 
          scale_fill_manual(values=c("#fed976", "#0570b0")) +  
          xlab("") + ylab("") + theme_jon + theme(plot.margin = unit(c(0.1, 0, 0.1, 0), "cm"), axis.text.y = element_blank()) + ggtitle(context)
      }
      altuna_plots[[plotnum]] = curPlot
    }
  }
}

## Create rows of plots
altunarow_SNPs = plot_grid(altuna_plots[[1]], altuna_plots[[2]], altuna_plots[[3]], altuna_plots[[4]], altuna_plots[[5]], altuna_plots[[6]], altuna_plots[[7]], 
          nrow=1, rel_widths=c(3.2, 2, 2, 2, 2, 2, 2), align='h', axis='b')
altunarow_INDELs = plot_grid(altuna_plots[[8]]+ ggtitle(""), altuna_plots[[9]]+ ggtitle(""), altuna_plots[[10]]+ ggtitle(""), altuna_plots[[11]]+ ggtitle(""), altuna_plots[[12]]+ ggtitle(""), altuna_plots[[13]]+ ggtitle(""), altuna_plots[[14]]+ ggtitle(""), nrow=1, rel_widths=c(3.2, 2, 2, 2, 2, 2, 2), align='h', axis='b')
altunarows = plot_grid(altunarow_SNPs, altunarow_INDELs, nrow=2)

## Export as png (too big for pdf)
png(file="figures/Fig4cd_altuna.png", width=1200, height=500)
altunarows
dev.off()

# INDEL Size Distribution

## Set up matrix
df_lengths = read.csv("tables/rtg_allele_lengths.csv")
df_lengthsm = melt(df_lengths[, c(1,2,3,5)], id.vars = c("length","platform"))
df_lengthsm$length = ifelse(df_lengthsm$variable == "delete", -(df_lengthsm$length), df_lengthsm$length)
df_lengthsm[df_lengthsm$platform=="PromethION",]$platform = "Nanopore"
df_lengthsm

## Plot
gg_lengthsINDEL = ggplot(subset(df_lengthsm, value>0), aes(x=length, y=value, color=platform)) +
  geom_point(size=1, show.legend=T) +
  stat_smooth(alpha=0.3, se=FALSE, span=0.3, show.legend=F) +
  geom_vline(xintercept=c(-50, 50), linetype="dashed", alpha=0.5) +
  theme_jon + theme(legend.position="right", legend.title=element_blank()) +
  xlab("Size (bp)") + ylab("Count") +
  scale_y_continuous(trans="log10") +
  xlim(-50, 50) +
  scale_color_manual(values = platform.colors) +
  geom_hline(yintercept=1)

## Export
pdf(file="figures/Fig4e_INDEL.pdf", width=10, height=3)
gg_lengthsINDEL
dev.off()

