# Dependencies
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(ComplexHeatmap)
library(circlize)
library(ggbeeswarm)
library(reshape2)
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
platform.colors = c("BGISEQ500"="#fd8d3c", "HiSeq2500"="#6baed6", "HiSeq4000"="#3182bd", "HiSeqX10"="#08519c", "MGISEQ2000"="#e6550d", "NovaSeq2x150"="#f768a1", "NovaSeq2x250"="#ae017e", "PacBio" = "#74c476", "PromethION"="#006d2c")
context.labels = c("Alu"="Alu", "L1"="L1", "L2"="L2",
                   "lowcomplexity"="Low Complexity",
                   "LTR"="LTR", "satellite"="Satellite", "simplerepeat"="Simple Repeat")
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

# Per-region Coverage comparisons

## Create data matrix
df = read.csv("tables/ucsc-coverage.csv")
df$locus = paste0(df$chr, "_", df$start, "_", df$end)
df = df[, c("context", "sample", "locus", "cov")]
df = separate(df, col = "sample", sep="_", into=c("platform","lab","genome","rep"), remove=F)

## Plot
gg_covWithSig = ggplot(df, aes(x=platform, y=ihs(cov))) +
  geom_violin(aes(fill=platform, alpha=0.85), show.legend=F) +
  theme_jon + theme(axis.title.y = element_text(size=12)) +
  ylab("Mean Cov (IHS)") + xlab("") +
  scale_y_continuous(breaks=c(0,ihs(10), ihs(25), ihs(50),ihs(100)), labels = c("0","10","25","50","100")) +
  coord_cartesian(ylim=c(0,6)) +
  scale_fill_manual(values=platform.colors) +
  stat_compare_means(method="wilcox.test", ref.group=".all.", label="p.signif", method.args = list(alternative = "greater"), hide.ns=T, size=3, label.y = 5.5) +
  facet_wrap(~context, nrow=2, labeller = labeller(context = context.labels)) 


## Export
pdf(file="figures/Fig2a_CovIHS_v2.pdf", width=7.5, height=5)
gg_covWithSig
dev.off()


## Create matrix of A versus B, pct of sites that A has better coverage than B
contexts = c("Alu", "L1", "L2", "lowcomplexity", "LTR", "satellite", "simplerepeat")
platforms = c("HiSeq2500", "HiSeq4000", "HiSeqX10", "BGISEQ500", "MGISEQ2000", "NovaSeq2x150", "NovaSeq2x250", "PacBio", "PromethION")
ava_matrix_all = data.frame()
for(context in contexts){
  ava_matrix = matrix(, nrow=9, ncol=9)
  for(plat1 in 1:8){
    for(plat2 in (plat1+1):9) {
      platform1 = platforms[plat1]
      platform2 = platforms[plat2]
      df2 = subset(df, (platform == platform1 | platform == platform2))
      df2 = df2[df2$context == context , ]
      df2 = df2 %>% group_by(platform, locus) %>% summarise(cov = mean(cov))
      df2 = dcast(df2, formula = locus ~ platform)
      df2$comp = ifelse(df2[, platform1] > df2[, platform2], platform1, ifelse(df2[, platform1] == df2[, platform2], "Even", platform2))
      df2[is.na(df2)] = 0
      df2[ df2[[platform1]]==0 & df2[[platform2]]==0, ]$comp = "BothNA"
      platform1count = nrow(subset(df2, comp == platform1))
      platform2count = nrow(subset(df2, comp == platform2))
      evencount = nrow(subset(df2, comp == "Even"))
      ava_matrix[plat1,plat2] = platform1count / (platform1count + platform2count + evencount)
      ava_matrix[plat2,plat1] = platform2count / (platform1count + platform2count + evencount)  
    }
  }
  rownames(ava_matrix) = platforms
  colnames(ava_matrix) = platforms
  
  ava_matrix_m= melt(ava_matrix)
  names(ava_matrix_m) = c("Sample1", "Sample2", "value")
  ava_matrix_m$context = context
  
  ava_matrix_all = rbind(ava_matrix_all, ava_matrix_m)
}
ava_matrix_all$barColor = ifelse(ava_matrix_all$value >= 0.5, "greater", "lesser")

## Plot
gg_avaCovAll = ggplot(ava_matrix_all, aes(x=Sample2, y=value)) +
  geom_bar(stat="identity", aes(fill=barColor, color=barColor), width = 1, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5) +
  theme_bw() +
  theme(legend.position="bottom",
      text = element_text(size=12),
      axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text.y = element_text(angle=0),
      plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("greater"="#3288bd", "lesser"="#d53e4f")) +
  scale_color_manual(values = replicate(106, "#ebebeb")) +
  xlab("") + ylab("") +
  coord_flip() +
  facet_grid(context~Sample1, labeller = labeller(context = context.labels)) 

## Export
pdf(file="figures/Supp__CovAVA.pdf", width=11, height=8)
gg_avaCovAll
dev.off()

## Plot all vs all
gg_avaCov = ggplot(ava_matrix_all, aes(x=context, y=value)) +
  geom_point(aes(fill=barColor), shape=21, show.legend = F, size=3, alpha=0.85) +
  geom_hline(yintercept = 0.5, linetype="dashed", alpha=0.5) +
  theme_jon +
  scale_fill_manual(values = c("greater"="#3288bd", "lesser"="#d53e4f")) +
  xlab("") + ylab("") +
  coord_flip() +
  facet_wrap(~Sample1, nrow=1)#, labeller = labeller(context = context.labels))

## Export
pdf(file="figures/Fig2b_CovAVA_v2.pdf", width=7.5, height=5)
gg_avaCov
dev.off()

# Coefficient of Coverage

## Set up matrix
df_cv = df %>% group_by(context, platform, locus) %>% summarise(coeffvar = sd(cov,na.rm=TRUE) / mean(cov,na.rm=TRUE))
df$platform = factor(df$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000", "PacBio", "PromethION"))
df_cv$platform = factor(df_cv$platform, levels=c("HiSeq2500", "HiSeq4000", "HiSeqX10", "NovaSeq2x150", "NovaSeq2x250", "BGISEQ500", "MGISEQ2000", "PacBio", "PromethION"))
df_cv

## Plot
gg_covCV = ggplot(subset(df_cv), aes(x=platform, y=coeffvar)) +
  geom_boxplot(aes(fill=platform), outlier.shape = NA, show.legend=F) +
  theme_jon +
  scale_fill_manual(values=platform.colors) +
  xlab("") + ylab("CV") +
  scale_y_continuous(breaks=seq(0,0.6,0.1)) +
  coord_trans(ylim=c(0,0.7)) +
  stat_compare_means(method="wilcox.test", ref.group=".all.", label="p.signif", method.args = list(alternative = "greater"), hide.ns=T, size=3, label.y = 5.5) +
  facet_wrap(.~context, nrow=2, labeller = labeller(context = context.labels)) 

## Export
pdf(file="figures/Fig2c_CovCV_v2.pdf", width=7.5, height=5)
gg_covCV
dev.off()
