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

umap_original <- readRDS("snATACseq/processed_data/UMAP_GeneScore_rm.rds")
commap <- readRDS(file="snATACseq/processed_data/Organoids_Integrated_RL1738.RDS")

source("code/predict_functions.R")

df <- umap_original[[3]]

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

predict_all(commap, test=commap$dataset=="query", 
            train=commap$dataset=="ref", knn=50, age_method="mode",
            method_all="mode")
predict_celltypes(commap, test=commap$dataset=="query", 
                  train=commap$dataset=="ref", knn=50, age_method="mode",
                  method_celltypes="mode")

table(predict_cells(commap, test=commap$dataset=="query", 
  train=commap$dataset=="ref", knn=50, col="Anno"))/sum(commap$dataset=="query")

p1 <- plot_histogram(commap, test=commap$dataset=="query", 
                     train=commap$dataset=="ref", knn=50, method="mode")
ggsave(p1, file="paper_figures/Fig6_histogram_organoid_predicted_age.svg",
       height=5, width=8)

commap_ref <- commap[commap$dataset=="ref", ]
commap_query <- commap[commap$dataset=="query", ]
commap_query$Stage <- "query"

g1 <- ggplot(commap_ref[sample(1:nrow(commap_ref),), ], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) +
    #scale_colour_manual(values=c(Fetal="#512568", 
    #    Neonatal="#443682", Infancy="#3D6B93",
    #     Childhood="#20988C", 
    #     Adolescence="#98CA43", Adult="#F9E51B")) +
    theme_classic() + 
    geom_point(data=commap_query, alpha=0.5, size=1, col="#990000")  +
    geom_density_2d(data=commap_query, color="#ff7b00", size=0.5, bins=15,
                    adjust=1/2) 

ggsave(g1, file="paper_figures/Fig6_organoid_integration_density_lines.png",
       width=7, height=7)

g1 <- ggplot(commap_ref[sample(1:nrow(commap_ref),), ], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) + theme_classic()

ggsave(g1, file="paper_figures/Fig6_organoid_integration_outline.png",
       width=7, height=7)
