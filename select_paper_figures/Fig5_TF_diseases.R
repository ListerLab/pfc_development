## TF for diseases

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

tf <- read.table("snATACseq/processed_data/tf_diseases.tsv", sep="\t")
tf <- tf[tf$padjust<0.05,]
tf$feature_short <- sapply(tf$feature, function(x) strsplit(x, "_")[[1]][1])

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")
a <- acast(tf,  feature_short ~ cell_type , mean, fill=1,
           drop=T, value.var="padjust")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
tf$feature_short <- factor(tf$feature_short,
                           levels=levels_motifs)


g1 <- ggplot(tf, aes(x=feature_short, y=cell_type, fill=padjust)) + geom_tile() +
    theme_classic() + 
    scale_fill_gradientn(colours =
      rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) 

ggsave(g1, file="paper_figures/Fig5_TF_diseases.svg", height=5, width=7)