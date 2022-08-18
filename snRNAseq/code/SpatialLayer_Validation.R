library(spatialLIBD)
library(ggplot2)
library(dplyr)

# 1. Load spatial data
sce <- fetch_data(type = 'sce')

# 2. Make pictures for HTR2A

dat <- data.frame(Layer=sce$layer_guess_reordered_short,
    HTR2A=sce@assays@data$logcounts[which(rowData(sce)$gene_name=="HTR2A"),],
                  Sample=sce$sample_name)

dat <- dat %>% group_by(Sample, Layer) %>% summarize(HTR2A=mean(HTR2A)) 
dat <- dat %>% filter(!is.na(Layer))

g1 <- ggplot(dat, aes(x=Layer, y=HTR2A)) + geom_boxplot() + 
    geom_point() + theme_classic() 
ggsave(g1, file=paste0("supp_figures/SuppFig5_HTR2A_spatial.svg"),
       height=4, width=7)

