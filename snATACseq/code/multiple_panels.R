# UMAP plot and composition plot

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

save_objs <- readRDS("snATACseq/processed_data/UMAP_GeneScore_rm.rds")
df <- save_objs[[3]]
df$Sample <- save_objs[[2]]$Sample
df$Stage <- save_objs[[2]]$Stage
df$Cluster <- save_objs[[2]]$Clusters
df$Anno <- save_objs[[2]]$Anno1
df$Group <- save_objs[[2]]$predictedGroup
colnames(df)[1:2] <- c("UMAP1", "UMAP2")
df$Age <- save_objs[[2]]$arcsin_ages
df$Age <- as.numeric(as.character(df$Age))

set.seed(10)
df <- df[sample(1:nrow(df), nrow(df)),]
df$Sample <- as.character(df$Sample)

## Subplot a)

### Annotation

col <- c(Astro= '#ffc857', `IN dev`='#c6a785', 
         `L2/3`= '#6e0614', L4= '#8b3843', `L5/6`='#a86a72',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN dev`='#e2cdd0',
         `MGE der`='#b44622',`CGE der`='#1c5701', Vas='#a3a3a3',
         `Low quality`='lightgrey')

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Anno), size=0.285, alpha=0.5) +
    theme_classic() +  scale_color_manual(values=col) + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))
ggsave(g1, file="paper_figures/Fig4_Anno.png", height=7, width=7)

### Maturation 

ages_dat <- data.frame(c("ga22","ga38","60d","1yr","10yr","20yr"),
   c(-0.59206997, -0.06903606,  0.29393111,  1.35923668,  3.59432467,
     4.28690546))
ages_dat <- ages_dat[order(as.numeric(ages_dat[,2])),]

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Age), size=0.285, alpha=0.5) +
    theme_classic() + scale_colour_gradientn(colours =
        grDevices::colorRampPalette(colors=c("#512568", "#443682", "#3D6B93",
        "#20988C", "#98CA43", "#F9E51B"))(20),
        breaks=as.numeric(ages_dat[,2]), labels=ages_dat[,1]) 

ggsave(g1, file="paper_figures/Fig4_Maturation.png", height=5.3, width=4.9)

## Subplot b)

### Composition

df1 <- as.data.frame(save_objs[[2]])[,c("Stage","Clusters", "Sample")]
df1$Anno <- cell_types[match(df1$Cluster, names(cell_types))]

relevel_dat <- unique(save_objs[[2]][, c("Sample", "arcsin_ages")])
new_levels <- relevel_dat$Sample[order(relevel_dat$arcsin_ages, decreasing=FALSE)]

df1$Stage <- factor(df1$Stage, levels = c("Fetal", "Neonatal", "Infancy", 
                                          "Childhood", "Adolescence", "Adult"))
df1$Anno <- factor(df1$Anno, levels=rev(c('L2/3', 'L4', 
                                          'L5/6', 'PN dev', 'MGE der', 'CGE der', 'IN dev', 
                                          'Astro', 'Oligo', 'OPC', 'Micro', 'Vas')))
df1$Sample <- factor(df1$Sample, levels = as.character(new_levels))
g1 <- ggplot(data = df1, mapping = aes(y = Stage, fill = Anno)) + 
    geom_bar(position = "fill", col="black", size=0.2) + 
    scale_fill_manual(values = col) + 
    scale_y_discrete(limits = rev(levels(df1$Stage))) + theme_classic() +
    guides(fill=FALSE)

ggsave(g1, file="paper_figures/Fig4_Composition.svg", 
       height=2.81*2, width=8.98*2, units="cm")

