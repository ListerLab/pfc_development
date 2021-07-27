### Distance to the TSS of CRE element 

library(SingleCellExperiment)
library(Polychrome)
library(rcartocolor)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(cba)
library(viridis)
library(cowplot)
library(ggpubr)

pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")


g1 <- ggplot(pearson_corr_sig , aes(x=dist+1)) + geom_density(fill="black") + 
    theme_classic() + scale_x_log10(breaks=c(1,50,1000,10000, 250000),
    labels=c(0,50,1000,10000, 250000)) + xlab("Distance to TSS (bp)") +
    theme(axis.title.y = element_blank())

df <- pearson_corr_sig %>% select("dist", "peak_name")
df$peak_name <- sapply(df$peak_name, function(x) 
    strsplit(x, ".", fixed=T)[[1]][1])
df <- df %>% mutate(type=case_when(dist==0 ~ "Geneic",
                                   (dist>0 & dist<=2000) ~ "Proximal",
                                   dist>1000 ~ "Distal"))

df1 <- df %>% group_by(type) %>% summarise(n=n())
df1$fraction = df1$n/ sum(df1$n)

# Make donut plot
g2 <- ggdonutchart(data=df1, "fraction", label="type",
                   color = "white", fill="type",
                   palette = brewer.pal(3, "Set3")) + guides(fill=FALSE)

gg <- ggdraw(g1) + 
    draw_plot(g2, 0.2, .25, .6, .6) 

ggsave(gg, file="paper_figures/Fig3h_TSSdist.svg", width=2.5*3.5, height=2.5*2)
