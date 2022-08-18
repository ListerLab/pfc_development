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
library(FNN)

source("code/predict_functions.R")


umap_original <- readRDS("processed_data/UMAP_GeneScore_all.rds")
commap <- readRDS(file="processed_data/Organoids_Integrated_All.RDS")

df <- umap_original[[2]]

col <- c(Astro= '#ffc857', `IN dev`='#c6a785', 
         `L2/3`= '#6e0614', L4= '#8b3843', `L5/6`='#a86a72',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN dev`='#e2cdd0',
         `MGE der`='#b44622',`CGE der`='#1c5701', Vas='#a3a3a3',
         `Low quality`='lightgrey')

ind <- match(rownames(commap), rownames(df))
commap$Age <- NA
commap$Age[!is.na(ind)] <- df$Age[ind[!is.na(ind)]]
commap$Anno <- NA
commap$Anno[!is.na(ind)] <- df$Anno[ind[!is.na(ind)]]
commap$Stage <- NA
commap$Stage[!is.na(ind)] <- df$Stage[ind[!is.na(ind)]]
commap$Stage <- factor(commap$Stage, levels=c("Fetal",
      "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")) 

commap_ref <- commap[commap$dataset=="ref", ]
commap_query <- commap[commap$dataset=="query", ]
commap_query$Stage <- "query"

g1 <- ggplot(commap_ref[sample(1:nrow(commap_ref),), ], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) +
    geom_point(data=commap_query, alpha=0.5, size=1, col="#990000") + 
    theme_classic() + 
    geom_density_2d(data=commap_query, color="#ff7b00", size=0.5, bins=15,
                    adjust=1/2)

ggsave(g1, file="supp_figures/SuppFig7_organoid_integration_density_lines_all.png",
       width=7, height=7)

g1 <- ggplot(commap_ref[sample(1:nrow(commap_ref),), ], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) + theme_classic()

ggsave(g1, file="supp_figures/SuppFig7_organoid_integration_outline_all.png",
       width=7, height=7)
