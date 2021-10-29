# Peak length

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

processed_data <- readRDS(file=
   "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_dge_stage_gr.Rds")

all_peaks <- as(processed_data, "GRangesList")
all_peaks <- unlist(all_peaks, recursive = TRUE)

dat3 <- data.frame(id=names(all_peaks), width=width(all_peaks),
                   cell_type=sapply(names(all_peaks), function(x) 
                       strsplit(x, ".", fixed=T)[[1]][1]))

g1 <- ggplot(dat3, aes(x=width)) + geom_density(fill="black") + 
    theme_classic() + ggtitle("All peaks") + scale_x_log10()
ggsave(g1, file="supp_figures/Fig4_PeakLength.svg", width=3.5, height=2)
