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
strata = c("SNP", "INS_5", "INS_6to15", "INS_15", "DEL_5", "DEL_6to15", "DEL_15")
platforms = c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000")

## create plots in a loop
for(stratum in strata){
  ### read in matrix
  df = read.csv(paste0("mendelian_violation/outs/forUpset__",stratum,".csv"))
  rownames(df) = df[,1]
  df = df[,-c(1)]
  
  ### summarize by platform
  df_sub = df
  for(platform in platforms){
    df_sub[[platform]] = rowMeans(df_sub[grepl(platform, names(df_sub))])
    df_sub[[platform]] = df_sub[[platform]] + 0.01  # to account for "round to even" rule; normally, round(0.5) = 0
    df_sub[[platform]] = round(df_sub[[platform]])
  }
  
  df_sub_platforms = df_sub[, platforms]
  df_sub_platforms[is.na(df_sub_platforms)] = 0
  
  ### remove all zeroes
  df_sub_platforms = df_sub_platforms[rowSums(df_sub_platforms) > 0, ]
  
  ### plot and export
  pdf(file=paste0("figures/Mendel_",stratum,".pdf"), width=6, height=4, onefile = F)
  print(upset(df_sub_platforms, 
      nsets=7,
      #nsets=6, 
      nintersects=12, 
      sets=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000"),
      show.numbers="no", 
      main.bar.color="#ea5d4e", 
      sets.bar.color="#317eab",
      text.scale = c(0, 1.5, 0, 1.2, 1.5, 0),
      point.size = 3, line.size = 1.5,
      #empty.intersections=NULL, 
      #order.by = "freq",
      keep.order = T,
      mainbar.y.label ="No. of Intersections", sets.x.label ="Set size"))
  grid.text(stratum,x = 0.65, y=0.95, gp=gpar(fontsize=20))
  dev.off()
}
