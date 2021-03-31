# Dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(reshape2)
library(ggsignif)
library(cowplot)
options(scipen=999)
theme_linedraw2 = theme_linedraw() + theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), strip.text.x=element_text(colour="black"), strip.text.y=element_text(colour="black"))
theme_jon = theme_linedraw2 + 
  theme(legend.position="bottom", 
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10), 
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5), 
        axis.text.y=element_text(size=10), 
        axis.title=element_text(size=10), 
        axis.title.y=element_text(vjust=1, size=9),
        axis.title.x = element_text(size=9),
        plot.title = element_text(size=12, hjust = 0.5), 
        strip.text = element_text(size=11), 
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank())

platform.colors = c("BGISEQ500"="#fd8d3c", "HiSeq2500"="#6baed6", "HiSeq4000"="#3182bd", "HiSeqX10"="#08519c", "MGISEQ2000"="#e6550d", "NovaSeq2x150"="#f768a1", "NovaSeq2x250"="#ae017e", "Genapsys"="#f3d04a", "S5"="#67a9cf", "Proton"="#1c9099", "PacBio" = "#8c510a", "PromethION"="#006d2c", "MinION"="#74c476", "Flongle"="#c7e9c0", "MiSeq"="#08519c", "S5"="#67a9cf", "PGM"="#1c9099")
# original blue-green for PacBio is #66c2a4; using a darker brown here to make clearer

context.labels = c("Alu"="Alu", "L1"="L1", "L2"="L2",
                   "lowcomplexity"="Low Complexity",
                   "LTR"="LTR", "satellite"="Satellite", "simplerepeat"="Simple Repeat")

## Set up data matrix
df = read.csv("tables/mismatch_UCSC.csv")
df = separate(df, col="sample", sep="_",  into=c("platform","lab","genome","rep"), remove=F)
df$platform = factor(df$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "Genapsys", "NovaSeq2x150", "NovaSeq2x250", "S5", "Proton", "PacBio", "Flongle", "MinION",  "PromethION"))

## Plot overall mismatch per context
gg_mismatchOverall = ggplot(subset(df), aes(x=platform, y=error)) + 
  geom_hline(yintercept = 0) +
  geom_bar(stat="identity", aes(fill=platform, group=sample), position="dodge", show.legend=F) +
  theme_jon + theme(axis.text.x = element_text(size=9), axis.title.y = element_text(size=12)) +
  scale_fill_manual(values=platform.colors) +
  ylab("Mismatch Rate") + xlab("") +
  coord_trans(y="sqrt") +
  scale_y_continuous(breaks=c(0.001, 0.01, 0.05, 0.1, 0.2)) +
  facet_wrap(~context, nrow=1, labeller = labeller(context = context.labels))

## Export
pdf(file="figures/Fig3a_MismatchByUCSC.pdf", width=12, height=3)
gg_mismatchOverall
dev.off()


# Mismatch by GC Content

## Set up data matrix
df_mhist = read.csv("tables/mhist.csv")
df_mhistm = melt(df_mhist, id.vars = c("sample","GCmax","base"), measure.vars = c("Match1", "Sub1", "Del1", "Ins1", "N1", "Other1"))
df_mhistm = df_mhistm %>% group_by(sample, GCmax, variable) %>%  summarize_at(vars(value), mean)
df_mhistm = separate(df_mhistm, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_mhistm$platform = factor(df_mhistm$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "Genapsys", "BGISEQ500", "MGISEQ2000", "PacBio", "Nanopore"))


## Plot
gg_mismatch_byGC = ggplot(subset(df_mhistm,  variable != "Other1" & variable != "N1"), 
       aes(x=as.factor(GCmax), y=value, fill=variable)) +
  geom_bar(stat="identity", position="fill", show.legend=F) +
  theme_jon + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=14)) +
  scale_fill_manual(name = "", labels = c("Match", "Substitution", "Deletion", "Insertion", "N", "Other"), 
                       values=c("Match1"="#b2df8a", "Sub1"="#e31a1c", "Del1"="#ff7f00", "Ins1"="#1f78b4", "N1"="#bdbdbd", "Other1"="#636363")) +
  facet_grid(.~platform) +
  coord_trans(y="sqrt") +
  xlab("GC%") + ylab("Mismatch Fraction") +
  scale_y_continuous(breaks=c(0.01, 0.05, 0.13, 0.25, 0.50, 1.0)) +
  scale_x_discrete(breaks = c(5, 30, 55, 80, 100), labels=c("0", "25", "50", "75", "100"))

## Export
pdf(file="figures/Fig3b_MismatchByGC.pdf", width=12, height=3)
gg_mismatch_byGC
dev.off()
```

# Mismatch by base

## Set up data matrix
df_mhistmBase = melt(df_mhist, id.vars = c("sample","GCmax","base"), measure.vars = c("Match1", "Sub1", "Del1", "Ins1", "N1", "Other1"))
df_mhistmBase = df_mhistmBase %>% group_by(sample, base, variable) %>%  summarize_at(vars(value), mean)
df_mhistmBase = separate(df_mhistmBase, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_mhistmBase$platform = factor(df_mhistmBase$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "Genapsys", "BGISEQ500", "MGISEQ2000", "PacBio", "Nanopore"))

## Plot
gg_mismatch_byBase = ggplot(subset(df_mhistmBase, variable != "Other1" & variable != "N1"), 
       aes(x=base, y=value, fill=variable)) +
  geom_bar(stat="identity", position="fill") +
  theme_jon + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=14)) +
  scale_fill_manual(name = "", labels = c("Match", "Substitution", "Deletion", "Insertion"), 
                       values=c("Match1"="#b2df8a", "Sub1"="#e31a1c", "Del1"="#ff7f00", "Ins1"="#1f78b4")) +
  facet_grid(.~platform, scales="free_x") +
  coord_trans(y="sqrt") +
  xlab("Base") + ylab("Mismatch Rate") +
  scale_y_continuous(breaks=c(0.01, 0.05, 0.13, 0.25, 0.50, 1.0))

## Export
pdf(file="figures/Fig3c_MismatchByBase.pdf", width=12, height=3)
gg_mismatch_byBase
dev.off()

# Homopolymer plot

## Set up data matrix
df_homo = read.csv("tables/mismatch_homo_v2.csv")
df_homo_m = melt(df_homo, id.vars = colnames(df_homo)[1:17], measure.vars = colnames(df_homo)[18:ncol(df_homo)], variable.name = "allinfo", value.name = "error" )
df_homo_m = separate(df_homo_m, col = "allinfo", into=c("platform", "lab", "genome", "rep"), sep = "_", remove=F)
df_homo_m

## Plot
gg_homo = ggplot(subset(df_homo_m, copyNum < 100), aes(x=copyNum, y=error, color=platform)) +
  stat_smooth(alpha=0.3, se=FALSE, span = 0.1, show.legend=F) +
  geom_hline(yintercept = 0) +
  theme_jon + theme(axis.title.y = element_text(size=12), axis.title.x = element_text(size=14)) +
  scale_color_manual(values=platform.colors) +
  scale_x_continuous(breaks=c(25,50,75,100)) +
  coord_trans(y="sqrt") + scale_y_continuous(breaks=c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.06, 0.09, 0.12)) +
  xlab("Homopolymer Length") + ylab("Mismatch Rate")

## Export
pdf(file="figures/Fig3d_Homopolymers.pdf", width=5, height=3)
gg_homo
dev.off()

# STRs plot

## Set up data matrix (visualize random 1M)
df_STRs = read.csv("tables/mismatch_strs_v2.csv")
df_STRs_m = melt(df_STRs, id.vars = colnames(df_STRs)[1:17], measure.vars = colnames(df_STRs)[18:ncol(df_STRs)], 
                 variable.name = "allinfo", value.name = "error" )
df_STRs_m = separate(df_STRs_m, col = "allinfo", into=c("platform", "lab", "genome", "rep"), sep = "_", remove=F)
df_STRs_m_1M = sample_n(subset(df_STRs_m, platform != "NextSeq"), 1000000)
df_STRs_m_1M$platform = factor(df_STRs_m_1M$platform, level=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "Genapsys", "BGISEQ500", "MGISEQ2000", "Flongle", "MinION", "PromethION", "PacBio"))
df_STRs_m_1M

## Plot
gg_STRs = ggplot(subset(df_STRs_m_1M, entropy < 10), aes(x=entropy, y=error, color=platform)) +
  stat_smooth(alpha=0.3, se=FALSE, span = 0.1) +
  geom_hline(yintercept = 0) +
  theme_jon + theme(legend.title=element_blank(), legend.position = "right", legend.text = element_text(size = 9), legend.key.size=unit(2, "point"), axis.title.y = element_text(size=12), axis.title.x = element_text(size=14)) +
  scale_color_manual(values=platform.colors) +
  coord_trans(y="sqrt") + scale_y_continuous(breaks=c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.06, 0.09, 0.12)) +
  xlab("Entropy") + ylab("")

## Export
pdf(file="figures/Fig3e_Entropy.pdf", width=6, height=3)
gg_STRs
dev.off()
