# Dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(cowplot)
library(ggpubr)

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
bacteria.colors = c("MiSeq"="#4daf4a", "Genapsys"="#ff7f00", "S5"="#377eb8", "PGM"="#984ea3")

normalized<-function(y) {
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)
  }


# Mapping Stats
read_labels = c("pctProperlyPaired"="Uniquely Mapped", "pctDuplicated"="Duplicates", "pctMultimap"="Multi-Mapped", "pctUnmap"="Unmapped")
readType_labels = c("short"="Short Reads", "long"="Long Reads", "shallow"="Shallow", "exome"="Exome", "bact"="Bacterial")
df_map = read.csv("tables/mapping_stats_all_v2.csv")
df_map = separate(df_map, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)

## label short, long, and exome reads
df_map$readType = ifelse(df_map$platform %in% c("PromethION", "PacBio"), "long", ifelse(df_map$platform %in% c("S5", "Proton"), "exome", ifelse(df_map$platform %in% c("Flongle", "MinION", "Genapsys"), "shallow", "short")))

## get total bases / mean depth
df_map$totalBases = as.numeric(df_map$readsTotal * df_map$avgLength)
df_map$meanDepth  = ifelse(df_map$readType == "exome", df_map$totalBases/60470914, df_map$totalBases/3100000000)

## get percentages
df_map$pctProperlyPaired = df_map$readsProperlyPaired / df_map$readsTotal * 100
df_map$pctDuplicated = df_map$readsDuplicated / df_map$readsTotal * 100
df_map$pctMultimap = df_map$readsMultimap / df_map$readsTotal * 100
df_map$pctUnmap = df_map$readsUnmap / df_map$readsTotal * 100

## re-order factors
df_map$readType = factor(df_map$readType, levels=c("short", "shallow", "long", "exome"))
df_map$platform = factor(df_map$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "Genapsys", "NovaSeq2x150", "NovaSeq2x250", "S5", "Proton", "PacBio", "Flongle", "MinION",  "PromethION"))

## adjust PB for aesthetics so points don't go over lines
df_map[df_map$platform=="PacBio", ]$pctProperlyPaired = df_map[df_map$platform=="PacBio", ]$pctProperlyPaired - 0.1

## Plot depth of coverage
gg_bwameth_depth_1 = ggplot(subset(df_map, readType=="short"), aes(x=platform, y=meanDepth)) +
  geom_bar(stat="identity", aes(fill=platform, color=sample), size=0.1, show.legend = F, position = "dodge") +
  theme_jon + theme(axis.text.x = element_text(size=10), axis.title.y = element_text(size=12)) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values = replicate(106, "#ffffff")) +
  facet_wrap(~readType, scales="free", nrow=1, labeller = labeller(readType = readType_labels)) +
  xlab("") + ylab("Mean Depth")
gg_bwameth_depth_2 = ggplot(subset(df_map, readType=="shallow"), aes(x=platform, y=meanDepth)) +
  geom_bar(stat="identity", aes(fill=platform, color=sample), size=0.1, show.legend = F, position = "dodge") +
  theme_jon + theme(axis.text.x = element_text(size=10)) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values = replicate(106, "#ffffff")) +
  facet_wrap(~readType, scales="free", nrow=1, labeller = labeller(readType = readType_labels)) +
  xlab("") + ylab("")
gg_bwameth_depth_3 = ggplot(subset(df_map, readType=="long"), aes(x=platform, y=meanDepth)) +
  geom_bar(stat="identity", aes(fill=platform, color=sample), size=0.1, show.legend = F, position = "dodge") +
  theme_jon + theme(axis.text.x = element_text(size=10)) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values = replicate(106, "#ffffff")) +
  facet_wrap(~readType, scales="free", nrow=1, labeller = labeller(readType = readType_labels)) +
  xlab("") + ylab("") 
gg_bwameth_depth_4 = ggplot(subset(df_map, readType=="exome"), aes(x=platform, y=meanDepth)) +
  geom_bar(stat="identity", aes(fill=platform, color=sample), size=0.1, show.legend = F, position = "dodge") +
  theme_jon + theme(axis.text.x = element_text(size=10)) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values = replicate(106, "#ffffff")) +
  facet_wrap(~readType, scales="free", nrow=1, labeller = labeller(readType = readType_labels)) +
  xlab("") + ylab("")
gg_bwameth_depth_5 = ggplot(subset(df_map_bact, readType=="single" | readType=="paired"), aes(x=platform, y=meanDepth)) +
  geom_bar(stat="identity", aes(fill=platform, color=sample), size=0.1, show.legend = F, position = position_dodge(width=1)) +
  theme_jon + theme(axis.text.x = element_text(size=10)) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values = replicate(106, "#ffffff")) +
  facet_wrap(~readType2, scales="free", nrow=1, labeller = labeller(readType2 = readType_labels)) +
  scale_y_log10() +
  xlab("") + ylab("") 
gg_bwameth_depth = plot_grid(gg_bwameth_depth_1, gg_bwameth_depth_2, gg_bwameth_depth_3, gg_bwameth_depth_4, gg_bwameth_depth_5, nrow=1, align = 'h', rel_widths=c(4, 1.5, 1.5, 2.5, 4))

## Export
pdf(file="figures/Fig1b_Depth.pdf", width=12, height=3)
gg_bwameth_depth
dev.off()


## Plot mapping stats
gg_bwameth_map = ggplot(df_map, aes(x=platform, y=pctProperlyPaired)) +
  geom_boxplot(aes(fill=platform, color=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform, shape=genome), alpha=0.5, show.legend=F, size=3) +
  geom_vline(xintercept = c(8.5, 10.5), linetype="dotted") +
  theme_jon +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24, "Eco"=21, "Ste"=21, "Pflu"=21, "Pool"=21)) +
  geom_hline(yintercept = 100) +
  xlab("") + ylab("% Read Pairs") +
  ggtitle("Properly Mapped")

gg_bwameth_multi = ggplot(df_map, aes(x=platform, y=pctMultimap)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform, shape=genome), alpha=0.5, show.legend=F, size=3) +
  geom_vline(xintercept = c(8.5, 10.5), linetype="dotted") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24, "Eco"=21, "Ste"=21, "Pflu"=21, "Pool"=21)) +
  theme_jon +
  xlab("") + ylab("") +
  ggtitle("Multi-Mapped")

gg_bwameth_dup = ggplot(df_map, aes(x=platform, y=pctDuplicated)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform, shape=genome), alpha=0.5, show.legend=F, size=3) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24, "Eco"=21, "Ste"=21, "Pflu"=21, "Pool"=21)) +
  theme_jon + 
  xlab("") + ylab("") +
  ggtitle("Duplicated")

gg_bwameth_unmapped = ggplot(df_map, aes(x=platform, y=pctUnmap)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform, shape=genome), alpha=0.5, show.legend=F, size=3, cex = 0.3) +
  geom_vline(xintercept = c(8.5, 10.5), linetype="dotted") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=platform.colors) +
  scale_color_manual(values=platform.colors) +
  scale_shape_manual(values = c("Mother" = 21, "Father" = 22, "Son" = 24, "Eco"=21, "Ste"=21, "Pflu"=21, "Pool"=21)) +
  theme_jon + 
  xlab("") + ylab("") +
  ggtitle("Unmapped")

human_row = plot_grid(gg_bwameth_map, gg_bwameth_multi, gg_bwameth_dup, gg_bwameth_unmapped, nrow=1, rel_widths = c(3.2, 3, 2, 3))
human_row

## Export
pdf(file="figures/Fig1c_MappingHuman.pdf", width=12, height=3)
human_row
dev.off()


# Mapping Rates for Bacterial Genomes

## bacterial data
df_map_bact = read.csv("tables/mapping_stats_bacteria.csv")
df_map_bact = separate(df_map_bact, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
## label paired and single end
df_map_bact$readType = ifelse(df_map_bact$platform == "MiSeq", "paired", "single")
df_map_bact$readType2 = "bact"
## get total bases / mean depth
df_map_bact$totalBases = as.numeric(df_map_bact$readsTotal) * as.numeric(df_map_bact$avgLength)
df_map_bact$meanDepth  = ifelse(df_map_bact$genome == "Eco", df_map_bact$totalBases/4641652, 
                                ifelse(df_map_bact$genome == "Ste", df_map_bact$totalBases/2499279, 
                                       ifelse(df_map_bact$genome == "Pflu", df_map_bact$totalBases/6511547, df_map_bact$totalBases/39879740)))
## get percentages
df_map_bact$pctProperlyPaired = df_map_bact$readsProperlyPaired / df_map_bact$readsTotal * 100
df_map_bact$pctDuplicated = df_map_bact$readsDuplicated / df_map_bact$readsTotal * 100
df_map_bact$pctMultimap = df_map_bact$readsMultimap / df_map_bact$readsTotal * 100
df_map_bact$pctUnmap = df_map_bact$readsUnmap / df_map_bact$readsTotal * 100
## re-order factors
df_map_bact$platform = factor(df_map_bact$platform, levels=c("MiSeq", "Genapsys", "S5", "PGM"))
df_map_bact$genome = factor(df_map_bact$genome, levels=c("Eco", "Ste", "Pflu", "Pool"))
df_map_bact$readType = factor(df_map_bact$readType, levels=c("paired", "single"))


## arrange plots
gg_bacteria_map = ggplot(subset(df_map_bact), aes(x=genome, y=pctProperlyPaired)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform), alpha=0.5, dodge.width=0.75, shape=21, show.legend=F, size=2) +
  theme_jon + 
  scale_fill_manual(values=bacteria.colors) +
  scale_color_manual(values=bacteria.colors) +
  geom_hline(yintercept = 100) +
  xlab("") + ylab("% Read Pairs") +
  ggtitle("Properly Mapped")

gg_bacteria_multi = ggplot(subset(df_map_bact), aes(x=genome, y=pctMultimap)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform), alpha=0.5, dodge.width=0.75, shape=21, show.legend=F, size=2) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=bacteria.colors) +
  scale_color_manual(values=bacteria.colors) +
  theme_jon + 
  xlab("") + ylab("") +
  ggtitle("Multi-Mapped")

gg_bacteria_dup = ggplot(subset(df_map_bact, readType=="paired"), aes(x=genome, y=pctDuplicated)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform), alpha=0.5, dodge.width=0.75, shape=21, show.legend=F, size=2) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=bacteria.colors) +
  scale_color_manual(values=bacteria.colors) +
  theme_jon + theme(legend.title = element_blank()) +
  xlab("") + ylab("") +
  ggtitle("Duplicated")

gg_bacteria_unmapped = ggplot(subset(df_map_bact), aes(x=genome, y=pctUnmap)) +
  geom_boxplot(aes(color=platform, fill=platform), alpha=0.5, show.legend=F) +
  geom_beeswarm(aes(fill=platform), alpha=0.5, dodge.width=0.75, shape=21, size=2, show.legend=F, cex = 0.3) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=bacteria.colors) +
  scale_color_manual(values=bacteria.colors) +
  theme_jon + theme(legend.title = element_blank()) +
  xlab("") + ylab("") +
  ggtitle("Unmapped")
bact_row = plot_grid(gg_bacteria_map, gg_bacteria_multi, gg_bacteria_dup, gg_bacteria_unmapped,  nrow=1, align='h', axis='b', rel_widths = c(3.2, 3, 2, 3))
bact_row

## Export
pdf(file="figures/Fig1d_MappingBact.pdf", width=12, height=2.5)
bact_row
dev.off()


# Insert Size Distrbution

## set up matrix
df_insert = read.csv("tables/insert_size_all.csv")
df_insert = df_insert %>% group_by(sampleOverall, size) %>% mutate(freq = mean(freq))
df_insert = df_insert[, c("sampleOverall", "size", "freq")]
colnames(df_insert) = c("sample", "size", "freq")
df_insert = separate(df_insert, col="sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
df_insert$platform = factor(df_insert$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "NovaSeq2x150", "NovaSeq2x250"))

label_insert = c("HiSeq2500" = "HiSeq\n2500", 
                 "HiSeq4000" = "HiSeq\n4000", 
                 "HiSeqX10" = "HiSeq\nX10", 
                 "BGISEQ500" = "BGISEQ\n500",
                 "MGISEQ2000" = "MGISEQ\n2000", 
                 "NovaSeq2x150" = "NovaSeq\n2x150", 
                 "NovaSeq2x250" = "NovaSeq\n2x250")

## Plot
gg_insert = ggplot(df_insert, aes(x=size, y=sample, height=freq, fill=platform)) +
  geom_ridgeline(scale=1, show.legend = F, alpha=0.7) +
  scale_fill_discrete(guide = guide_legend(reverse=TRUE)) +
  theme_bw() + 
  theme(axis.title.y=element_blank(), 
    axis.ticks.y=element_blank(), 
    panel.spacing.y=unit(0,"pt"), 
    axis.text.y=element_blank(), 
    text = element_text(size=13),
    legend.text=element_text(size=10)) + 
  facet_grid(platform+lab~., scales="free_y", switch="y",
             labeller = labeller(platform = label_insert)) + 
  xlab("Insert size")

## Export
pdf(file="figures/InsertSize.pdf", width=6, height=10)
gg_insert
dev.off()


# GC Content Heatmap (Human)

## Set up matrix
df_gc = read.csv("tables/GC_rawCoverage_all.csv")
df_gc = df_gc %>% group_by(sampleOverall, GC) %>% mutate(normcov = mean(normcov))
df_gc = df_gc[, c("sampleOverall", "GC", "normcov")]
colnames(df_gc) = c("sample", "GC", "normcov")

## go wide
df_gc_tmp = dcast(data = df_gc, formula = ...~sample, value.var="normcov", fun.aggregate=sum)

## normalize values
df_gc2 = apply(df_gc_tmp[, colnames(df_gc_tmp)[2:ncol(df_gc_tmp)]],2,normalized)

## melt again
df_gc2 = melt(df_gc2)
colnames(df_gc2) = c("GC", "Sample", "Reads_Used")
df_gc2$GC = df_gc2$GC - 1
df_gc2 = separate(df_gc2, col=Sample, into=c("Platform", "Lab", "Genome", "Rep"), remove=F)

## go wide again
df_gc = dcast(data = df_gc2, formula = ...~GC, value.var="Reads_Used",  fun.aggregate=mean)

## create numeric matrix
df_gc_m = as.matrix(df_gc[,6:ncol(df_gc)])
df_gc_m[is.na(df_gc_m)] = 0
rownames(df_gc_m) = df_gc[,1]
df_gc_m[1:5, 1:10]

## add WGS versus exome
df_gc$capture = ifelse(df_gc$Platform %in% c("Proton", "S5"), "exome", "WGS")
df_gc$capture = factor(df_gc$capture, levels=c("WGS", "exome"))

## heatmap annotations
ha_names = HeatmapAnnotation(Platform = df_gc$Platform,
                             col = list(Platform = c("HiSeq2500" = "#9ecae1", 
                                                     "HiSeq4000" = "#6baed6", 
                                                     "HiSeqX10" = "#3182bd", 
                                                     "NovaSeq2x150" = "#ae017e", 
                                                     "NovaSeq2x250" = "#ae017e", 
                                                     "BGISEQ500" = "#fd8d3c", 
                                                     "MGISEQ2000" = "#e6550d", 
                                                     "Genapsys" = "#f3d04a", 
                                                     #"Flongle" = "#a1d99b", 
                                                     #"GridION" = "#41ab5d", 
                                                     #"PromethION" = "#006d2c", 
                                                     #"PacBioCLR" = "#ef3b2c", 
                                                     #"PacBioCCS" = "#a50f15",
                                                     "Proton" = "#fc8d59", "S5" = "#fee08b")), which="row")
ha_genomes = HeatmapAnnotation(Genome = df_gc$Genome,
                             col = list(Genome = c("Mother"="#FBB4AE", 
                                                   "Father"="#B3CDE3", 
                                                   "Son"="#CCEBC5")), which="row")
## plot
hm = Heatmap(subset(df_gc_m),
             col = colorRamp2(c(0,0.5,1), c( "#2c7bb6","#ffffbf","#d7191c")),
             cluster_columns = FALSE,
             show_column_names = TRUE,
             column_title = "Human Normalized Coverage across GC (WGS and Exome)",
             split = df_gc$capture,
             cluster_rows = T,
             column_title_gp = gpar(fontsize = 14),
             column_names_side = "top",
             column_labels = c("0",".",".",".",".",".",".",".",".",".",".","10",".",".",".",".",".",".",".",".",".","20",".",".",".",".",".",".",".",".",".","30",".",".",".",".",".",".",".",".",".","40",".",".",".",".",".",".",".",".",".","50",".",".",".",".",".",".",".",".",".","60",".",".",".",".",".",".",".",".",".","70",".",".",".",".",".",".",".",".",".","80",".",".",".",".",".",".",".",".",".","90",".",".",".",".",".",".",".",".","100"),
            
                 column_names_rot = 45,
                 column_names_gp = gpar(fontsize=12),
                 show_row_names = TRUE,
                 #bottom_annotation = ha_dens,
                 show_heatmap_legend = T,
                 heatmap_legend_param = list(
                   title = "Normalized\nCoverage"
                 ))


## Export
pdf(file="figures/Supp_HeatmapHuman.pdf", width=7.5, height=10)
draw(hm + ha_names + ha_genomes, merge_legend = TRUE, heatmap_legend_side="right")
dev.off()


# GC Content Heatmap (Bacteria)

## set up matrix
df_gc_bact = read.csv("tables/bacteria_GC_rawCoverage_all.csv")
#df_gc_bact = read.csv("tables/GC_rawCoverage_all.csv")
df_gc_bact = df_gc_bact %>% group_by(sampleOverall, GC) %>% mutate(normcov = mean(normcov))
df_gc_bact = df_gc_bact[, c("sampleOverall", "GC", "normcov")]
colnames(df_gc_bact) = c("sample", "GC", "normcov")

df_bac_gc_tmp = dcast(data = df_gc_bact, formula = ...~Sample, value.var="normcov", fun.aggregate=sum)

df_bac_gc2 = apply(df_bac_gc_tmp[, colnames(df_bac_gc_tmp)[2:ncol(df_bac_gc_tmp)]],2,normalized)

df_bac_gc2 = melt(df_bac_gc2)
colnames(df_bac_gc2) = c("GC", "Sample", "Reads_Used")
df_bac_gc2$GC = df_bac_gc2$GC - 1
df_bac_gc2 = separate(df_bac_gc2, col=Sample, into=c("Platform", "Lab", "Genome", "Rep"), remove=F)

## convert tall to wide
df_gc_bact = dcast(data = df_bac_gc2, formula = ...~GC, value.var="Reads_Used",  fun.aggregate=mean)

## create metadata
#df_gc_bact = separate(df_gc_bact, col="Sample", sep="_", into=c("platform", "lab", "genome", "rep"), remove=F)
#df_gc_bact$platform = factor(df_gc_bact$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "NovaSeq2x150", "NovaSeq2x250", "Genapsys", "Proton", "S5"))

# create numeric matrix
df_gc_bact_m = as.matrix(df_gc_bact[,6:ncol(df_gc_bact)])
df_gc_bact_m[is.na(df_gc_bact_m)] = 0
rownames(df_gc_bact_m) = df_gc_bact[,1]

## heatmap annotations
ha_names = HeatmapAnnotation(Platform = df_gc_bact$Platform,
                             col = list(Platform = c("MiSeq"="#7fc97f", "S5"="#fdc086", "PGM"="#386cb0", "Genapsys"="#f0027f")), which="row")

ha_genomes = HeatmapAnnotation(Genome = df_gc_bact$Genome,
                             col = list(Genome = c("Mother"="#FBB4AE", 
                                                   "Father"="#B3CDE3", 
                                                   "Son"="#CCEBC5")), which="row")


## Plot
hm = Heatmap(subset(df_gc_bact_m),
             col = colorRamp2(c(0,0.5,1), c( "#2c7bb6","#ffffbf","#d7191c")),
             cluster_columns = FALSE,
             show_column_names = TRUE,
             column_title = "Bacterial Normalized Coverage across GC",
             cluster_rows = T,
             split = factor(df_gc_bact$Genome, levels=c("Ste", "Eco", "Pflu", "Pool")),
             column_title_gp = gpar(fontsize = 14),
             column_names_side = "top",
             column_labels = c("0",".",".",".",".",".",".",".",".",".",".","10",".",".",".",".",".",".",".",".",".","20",".",".",".",".",".",".",".",".",".","30",".",".",".",".",".",".",".",".",".","40",".",".",".",".",".",".",".",".",".","50",".",".",".",".",".",".",".",".",".","60",".",".",".",".",".",".",".",".",".","70",".",".",".",".",".",".",".",".",".","80",".",".",".",".",".", "87"),
            
                 column_names_rot = 45,
                 column_names_gp = gpar(fontsize=8),
                 show_row_names = TRUE,
                 #bottom_annotation = ha_dens,
                 show_heatmap_legend = T,
                 heatmap_legend_param = list(
                   title = "Normalized\nCoverage"
                 ))

## Export
pdf(file="figures/Supp__HeatmapBact.pdf", width=7.5, height=7)
draw(hm + ha_names, merge_legend = TRUE, heatmap_legend_side="right")
dev.off()

# Entropy plot
bact_entropy = as.data.frame(matrix(nrow=19, ncol=4))
colnames(bact_entropy) = c("Ste", "Pool", "Eco", "Pflu")

entropy <- function(p) rowSums(-(p * log(p)), na.rm = T)

bact_entropy$Ste = entropy(df_gc_bact[df_gc_bact$Genome=="Ste", 6:82])[1:19]
bact_entropy$Pool = entropy(df_gc_bact[df_gc_bact$Genome=="Pool", 6:82])[1:19]
bact_entropy$Eco = entropy(df_gc_bact[df_gc_bact$Genome=="Eco", 6:82])[1:19]
bact_entropy$Pflu = entropy(df_gc_bact[df_gc_bact$Genome=="Pflu", 6:82])[1:19]

bact_entropym = melt(bact_entropy)
colnames(bact_entropym) = c("Genome", "entropy")

## Export
pdf(file="figures/Supp__HeatmapBactEntropy.pdf", width=5, height=3)
ggplot(bact_entropym, aes(x=Genome, y=entropy)) + geom_boxplot(aes(fill=Genome)) + geom_point(shape=21) + theme_jon
dev.off()
