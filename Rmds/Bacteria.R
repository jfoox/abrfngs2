# Dependencies 
library(reshape2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(circlize)
library(ggbeeswarm)
library(viridis)
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
platform.colors = c("Genapsys"="#f3d04a", "S5"="#67a9cf", "Proton"="#1c9099", "PGM"="#1c9099", "MinION"="#74c476", "Flongle"="#c7e9c0")

# Metagenomic Pool classification

## Set up data matrix
df_bacmap = read.csv("tables/bacteria_pool_mosdepth_nosupp_v2.csv")
df_bacmap = separate(df_bacmap, col="sample", sep="_", into=c("platform","lab","pool","rep"), remove=F)
df_bacmap$platform = factor(df_bacmap$platform, levels=c("MiSeq","Genapsys","S5","PGM","Flongle","MinION"))
df_bacmap$taxon = factor(df_bacmap$taxon, levels=c("S_epidermidis","E_faecalis", "P_haloplanktis", "H_halophilus", "B_subtilis", "E_coli", "P_fluorescens",
                                                   "C_violaceum", "H_volcanii", "M_luteus"))

## Plot
gg_bacmap = ggplot(df_bacmap, aes(x=sample, y=pct, fill=taxon)) +
  geom_bar(stat="identity", position = position_stack(reverse=TRUE)) +
  facet_grid(platform~., scales="free", space = "free") +
  theme_jon +  theme(legend.position = "bottom", legend.title = element_blank(),
                     axis.text.x = element_text(angle=0, hjust=0.5),
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     strip.text.y = element_text(angle=0)) +
  scale_fill_brewer(palette="Spectral", labels = c(expression(italic("S. epidermidis")),
                                                   expression(italic("E. faecalis")),
                                                   expression(italic("P. haloplanktis")),
                                                   expression(italic("H. halophilus")),
                                                   expression(italic("B. subtilis")),
                                                   expression(italic("E. coli")),
                                                   expression(italic("P. fluorescens")),
                                                   expression(italic("C. violaceum")),
                                                   expression(italic("H. volcanii")),
                                                    expression(italic("M. luteus")))) +
  scale_y_continuous(position = "right") +
  xlab("") + ylab("") +
  coord_flip()

## Export
pdf(file="figures/Figure6a_barplots.pdf", width=7.5, height=5)
gg_bacmap
dev.off()

## Reformat
df_bacmap2 = dcast(df_bacmap, formula = taxon~platform, value.var = "pct", fun.aggregate = mean)
df_bacmap2_cor = cor(df_bacmap2[2:ncol(df_bacmap2)])

## Plot and export
pdf(file="figures/Figure6c_heatmap.pdf", width=7.5, height=3.5)
Heatmap(df_bacmap2_cor,
        col=viridis(100),
        name="Spearman",
         cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(sprintf("%.3f", df_bacmap2_cor[i, j]), x, y, gp = gpar(fontsize = 10, fontcolor="white"))
})
dev.off()


# Boxplot per genome

## Data
df_bacmap2 = df_bacmap[, c(2,6,8)] %>% group_by(platform, taxon) %>% summarise_all(funs(mean))

## Plot
gg_bacmapTaxon = ggplot(df_bacmap, aes(x=taxon, y=pct)) +
  geom_boxplot(aes(fill=taxon)) +
  #geom_beeswarm(aes(color=platform), dodge.width=0.75, show.legend=F, size=2, cex = 0.1) +
  geom_point(aes(fill=taxon), size=2, shape=21) +
  theme_jon + theme(legend.position="none", axis.text.x = element_text(angle=0, size=10, hjust=0.5, vjust=0.5), axis.ticks.x = element_blank()) +
  scale_x_discrete(labels=c("S_epidermidis"="S. epid\n28.0%\n(+)","E_faecalis"="E. faec\n37.2%\n(+)", "P_haloplanktis"="P. halo\n40.1%\n(-)", "H_halophilus"="H. halo\n41.6%\n(+)", "B_subtilis"="B. subt\n43.5%\n(+)", "E_coli"="E. coli\n50.8%\n(-)", "P_fluorescens"="P. fluo\n61.4%\n(-)","C_violaceum"="C. viol\n64.8%\n(-)", "H_volcanii"="H. volc\n65.5%\n(+)", "M_luteus"="M. lute\n72.0%\n(+)")) +
  geom_hline(yintercept = 10, linetype="dashed", alpha=0.5) +
  xlab("") + ylab("") +
  scale_fill_brewer(palette="Spectral") +
  scale_color_manual(values = platform.colors)

## Export
pdf(file="figures/Figure6c_boxplots.pdf", width=7.5, height=3)
gg_bacmapTaxon
dev.off()
```
