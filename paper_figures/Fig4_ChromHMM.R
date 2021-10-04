# ChromHMM plot

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

dat2 <- read.table("snATACseq/processed_data/ChromHMM_AllPeaks.txt", sep="\t")
dat2$Tissue <- rownames(dat2)
dat2 <- melt(dat2, id="Tissue")
dat2$Tissue <- factor(dat2$Tissue, levels=rev(c("Brain Angular Gyrus", 
    "Brain Anterior Caudate", "Brain Cingulate Gyrus", "Brain Germinal Matrix",
    "Brain Hippocampus Middle", "Brain Inferior Temporal Lobe", 
    "Brain Dorsolateral Prefrontal Cortex", "Brain Substantia Nigra",
    "Fetal Brain Male", "Fetal Brain Female", "Fetal Lung", "Fetal Thymus",
    "Esophagus", "Left Ventricle", "Liver", "Lung", "Ovary", "Pancreas", 
    "Placenta", "Thymus",  "Spleen")))

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
  "#08306B")

g1 <- ggplot(dat2, aes(y=Tissue, x=variable, fill=value)) + 
    geom_tile() + theme_classic() + 
    scale_fill_gradientn(colours =
         grDevices::colorRampPalette(colors=blue_colours)(20)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) 
ggsave(g1, file="paper_figures/Fig4e=_ChromHMM.svg", width=3.5, height=4)
