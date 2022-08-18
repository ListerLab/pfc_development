library(readxl)
library(ggplot2)

paths <- "snRNAseq/processed_data/"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
dev_no <- table(gene_trends$major_clust)
names(dev_no) <- gsub("-", "/", names(dev_no), fixed=TRUE)

no_trajectory <- read_xlsx("annotation/number_trajectories.xlsx")
ind <- match(no_trajectory$`Developmental Trajectory`, names(dev_no))
no_trajectory$`Number of devDEGs` <- dev_no[ind]

col <- c(Astro= '#ffc857', CGE_dev='#c6d5c0', ID2= '#558140', 
         `L2/3_CUX2`= '#6e0614', L4_RORB= '#8b3843', `L5/6_THEMIS`='#a86a72',
         `L5/6_TLE4`='#c59ba1', LAMP5_CA1='#8eab80', MGE_dev='#ecd1c8',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', PN_dev='#e2cdd0',
         PV='#c77459', PV_SCUBE3='#daa290',`Poor-Quality`='#ffffff',
         SST='#b44622',VIP='#1c5701', Vas='#a3a3a3')

gg1 <- ggplot(no_trajectory, aes(x=`Number of Nuclei`, y=`Number of devDEGs`)) +
    geom_point(aes(col=`Developmental Trajectory`), size=4) + theme_minimal() +
    scale_color_manual(values=col) + guides(col=FALSE)

ggsave(gg1, file="supp_figures/SuppFig3_devDEGsvsNuclei.svg", height=4, width=4)


