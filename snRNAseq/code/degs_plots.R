library(ggplot2)
library(dplyr)

all_res <- readRDS("snRNAseq/processed_data/All_trajectories_major.RDS")

res <- sapply(all_res, function(x) length(x))
dat <- data.frame(Trajectory=names(res), 
                  `Number of differentially expressed genes`=res,
                  check.names = FALSE)
dat$Trajectory[dat$Trajectory=="L2-3_CUX2"] <-"L2/3_CUX2"
dat$Trajectory[dat$Trajectory=="L5-6_TLE4"] <-"L5/6_TLE4"
dat$Trajectory[dat$Trajectory=="L5-6_THEMIS"] <-"L5/6_THEMIS"
dat$Trajectory <- factor(dat$Trajectory, 
  levels=c("L2/3_CUX2", "L4_RORB", "L5/6_TLE4", "L5/6_THEMIS",
    "PV", "SST", "PV_SCUBE3", "VIP", "ID2", "LAMP5_CA1",
    "Astro", "Oligo", "OPC", "Micro", "Vas"))

col <- c(Astro= '#ffc857', CGE_dev='#c6d5c0', ID2= '#558140', 
         `L2/3_CUX2`= '#6e0614', L4_RORB= '#8b3843', `L5/6_THEMIS`='#a86a72',
         `L5/6_TLE4`='#c59ba1', LAMP5_CA1='#8eab80', MGE_dev='#ecd1c8',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', PN_dev='#e2cdd0',
         PV='#c77459', PV_SCUBE3='#daa290',`Poor-Quality`='#ffffff',
         SST='#b44622',VIP='#1c5701', Vas='#a3a3a3')

col <- col[match(dat$Trajectory,names(col))]
gg1 <- ggplot(dat, aes(x=Trajectory, 
                       y=`Number of differentially expressed genes`, fill=Trajectory)) +
    geom_col() + theme_classic() + scale_fill_manual(values=col) + 
    guides(fill=FALSE) +
    geom_text(aes(label=`Number of differentially expressed genes`), 
              position=position_dodge(width=0.9), vjust=-0.25) + guides(fill=FALSE)
ggsave(gg1, file="main_figures/Number_DEG.svg", height=5, width=6)


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

ggsave(gg1, file="supp_figures/devDEGsvsNuclei.svg", height=4, width=4)

g1 <- ggplot(data = colData_sce, mapping = aes(y = batch, fill = major_clust)) + 
  geom_bar(position = "fill", col="black", size=0.2) + 
  scale_fill_manual(values = col) + 
  theme_classic() +
  guides(fill=FALSE)







