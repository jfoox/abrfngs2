# Dependencies
library(reshape2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(UpSetR)
library(grid)
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
        plot.title = element_text(size=10, hjust = 0.5), 
        strip.text = element_text(size=11), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.grid.major = element_line(colour = "grey98"), 
        panel.grid.minor = element_blank())
platform.colors = c("BGISEQ500"="#fd8d3c", "HiSeq2500"="#6baed6", "HiSeq4000"="#3182bd", "HiSeqX10"="#08519c", "MGISEQ2000"="#e6550d", "NovaSeq2x150"="#f768a1", "NovaSeq2x250"="#ae017e", "Genapsys"="#f3d04a", "S5"="#67a9cf", "Proton"="#1c9099", "PacBio" = "#66c2a4", "PromethION"="#006d2c", "MinION"="#74c476", "Flongle"="#c7e9c0", "MiSeq"="#08519c", "S5"="#67a9cf", "PGM"="#1fb597")

# Chromosomal mapping
df_chr = read.csv("tables/mapping_byChrType.csv")
df_chr = separate(df_chr, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_chr$decoyPct = df_chr$decoy / (df_chr$decoy + df_chr$canonical) * 100
df_chr$platform = factor(df_chr$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "Genapsys", "NovaSeq2x150", "NovaSeq2x250", "S5", "Proton", "PacBio", "Flongle", "MinION", "PromethION"))

## Plot
gg_chr = ggplot(df_chr, aes(x=platform, y=decoyPct)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(aes(fill=platform), show.legend=F) +
  geom_beeswarm(aes(fill=platform), dodge.width=0.75, shape=21, show.legend=F, size=3, cex=0.4) +
  scale_fill_manual(values=platform.colors) +
  ylab("% Total Reads") + xlab("") +
  theme_jon + theme(axis.title.y = element_text(size=14)) + ggtitle("Decoy Capture")

## Export
pdf(file="figures/Supp__ChrMapping.pdf", width=4, height=5)
gg_chr
dev.off()
