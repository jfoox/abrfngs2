# Dependencies 
options(scipen=999)
library(reshape2)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(viridis)
library(UpSetR)
library(grid)

# Variant Allele Heatmap

## Data
df_GT = read.csv("vafrobe/all_WGS.chr1.GT.toMatrix", sep="\t")
df_GT_m = as.matrix(df_GT[, 2:ncol(df_GT)])
rownames(df_GT_m) = df_GT[,1]

## Metadata
df_GT_meta = as.data.frame(colnames(df_GT)[2:ncol(df_GT)])
colnames(df_GT_meta) = "sample"
df_GT_meta = separate(df_GT_meta, col="sample", into=c("platform", "lab"," genome", "rep"), remove=F)

platform.colors = c("BGISEQ500"="#fd8d3c", "HiSeq2500"="#6baed6", "HiSeq4000"="#3182bd", "HiSeqX10"="#08519c", "MGISEQ2000"="#e6550d", "NovaSeq2x150"="#f768a1", "NovaSeq2x250"="#ae017e", "Genapsys"="#f3d04a", "S5"="#67a9cf", "Proton"="#1c9099", "PacBio" = "#66c2a4", "PromethION"="#006d2c", "MinION"="#74c476", "Flongle"="#c7e9c0", "MiSeq"="#08519c", "S5"="#67a9cf", "PGM"="#1fb597")
genome.colors = c ("Mother" = "#fb9a99", "Father"="#a6cee3", "Son"="#b2df8a")
lab.colors = c("LAB01" = "#fdd49e", "LAB02"="#fc8d59", "LAB03"="#990000")
rep.colors = c("REP01" = "#f7f7f7", "REP02" = "#d9d9d9", "REP03"="#969696")

## Plot
hm = Heatmap(df_GT_m[1:5000, ],
        name="GT", show_heatmap_legend = T,
        #col = colorRamp2(c(0,0.5,1), c( "#2c7bb6","#ffffbf","#d7191c")),
        col = structure(c( "#2c7bb6","#fc8d59","#d7191c"), names = c(0, 0.5, 1)),
        cluster_rows = T,
        show_row_names = F,
        top_annotation = HeatmapAnnotation(Genome = df_GT_meta$` genome`,
                                           Platform = df_GT_meta$platform,
                                           Lab = df_GT_meta$lab,
                                           Rep = df_GT_meta$rep,
                                           col = list(Platform = platform.colors,
                                                      Lab = lab.colors,
                                                      Genome = genome.colors,
                                                      Rep = rep.colors),
                                           show_legend = T),
        na_col = "#737373",
        cluster_columns = T,
        show_column_names = F,
        column_names_gp = gpar(fontsize=8))
        #col = variant_colors,
        #heatmap_legend_param = list(title = "Variant", at = c(-1, -0.5, 0, 0.5, 1), 
        #labels = c("No Info", "Reference (Low)", "Reference", "Variant (Low)", "Variant"), border = "black"))


## Export
pdf(file="figures/VAFrobe.pdf", width=7.5, height=7)
draw(hm)
dev.off()

# Correlation plot

## Calculate
df_GT_cor_HG002 = cor(df_GT_m[, grepl("Son", colnames(df_GT_m))], use = "complete.obs")
Heatmap(df_GT_cor_HG002,
        col=viridis(100),
        name="Spearman",
        show_row_names = F,
        show_column_names = F,
                top_annotation = HeatmapAnnotation(Genome = subset(df_GT_meta, ` genome` == "Son")$` genome`,
                                           Platform = subset(df_GT_meta, ` genome` == "Son")$platform,
                                           Lab = subset(df_GT_meta, ` genome` == "Son")$lab,
                                           Rep = subset(df_GT_meta, ` genome` == "Son")$rep,
                                           col = list(Platform = platform.colors,
                                                      Lab = lab.colors,
                                                      Genome = genome.colors,
                                                      Rep = rep.colors),
                                           show_legend = T))


# UpSet plots

## Set up data matrix
genomes   = c("Mother", "Father", "Son")
df_GT_upset = df_GT
rownames(df_GT_upset) = df_GT_upset[, 1]
df_GT_upset = df_GT_upset[, 2:ncol(df_GT_upset)]
df_GT_upset[is.na(df_GT_upset)] = 0
df_GT_upset[1:ncol(df_GT_upset)] = lapply(df_GT_upset[1:ncol(df_GT_upset)], function(x) ifelse(x > 0, 1, 0) )

## loop to make plots and export
genomes = c("Mother", "Father", "Son")
for(genome in genomes) {
  
  ### list platforms by genome
  if(genome=="Son") {
    platforms = c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq", "BGISEQ500", "MGISEQ2000", "PromethION", "PacBio")
  } else {
    platforms = c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq", "BGISEQ500", "MGISEQ2000")
  }
  
  ### create a subset just for a given genome
  df_GT_upset_sub = df_GT_upset[, grepl(genome, colnames(df_GT_upset))]
  df_GT_upset_sub[is.na(df_GT_upset_sub)] = 0
  
  ### summarize variants by platform
  for(platform in platforms) {
    df_GT_upset_sub[[platform]] = rowMeans(df_GT_upset_sub[grepl(platform, names(df_GT_upset_sub))])
    df_GT_upset_sub[[platform]] = df_GT_upset_sub[[platform]] + 0.01  # to account for "round to even" rule; normally, round(0.5) = 0
    df_GT_upset_sub[[platform]] = round(df_GT_upset_sub[[platform]])
  }
  df_GT_upset_sub_platforms = df_GT_upset_sub[, platforms]
  df_GT_upset_sub_platforms[is.na(df_GT_upset_sub_platforms)] = 0
  
  ### create UpSet plot
  pdf(file=paste0("figures/upset__GT__",genome,".pdf"), width=7.5, height=6, onefile=F)
  print(upset(df_GT_upset_sub_platforms,
        sets = rev(platforms),
        keep.order = T,
        nsets = 6, nintersects = 10,
        order.by = "freq",
        text.scale = c(0, 2, 0, 1, 1.7, 1.3),
        point.size = 3, line.size = 1.5))
  grid.text(paste0("Variants across ",genome, " Samples (Chr 1)"), x=0.65, y=0.95, gp=gpar(fontsize=14))
  dev.off() 
}
