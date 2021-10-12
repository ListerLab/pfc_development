# Plot motif enrichment
library(readxl)
library(ggplot2)
library(chromVAR)
library(JASPAR2018)
library(seqLogo)
library(viridis)
library(cba)
library(reshape2)
library(dplyr)
library(survcomp)

trajectories <- c("L2_3", "L4", "L5_6","CGE_der", "MGE_der", "Astro", "Oligo")


blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

dat <- read.csv("supp_table/tf_general_trends_trajectories_sig.csv")


for(i in 1:length(trajectories)){
    
    print(i)
    not_exprs <- tf_expr[trajectories[i]]
    tmp_feature <- unique(dat$feature[dat$trajectory==trajectories[i]])
    tmp_feature <- tmp_feature[tmp_feature %in% not_exprs]
    
    if(length(tmp_feature)>0){
        
        ind <- which(dat$feature %in% tmp_feature &
                         dat$trajectory == trajectories[i])
        dat <- dat[-ind, ]
    }
    
    
}

dat$feature_short <- gsub("..", "+", dat$feature, fixed=T)
dat$feature_short <- gsub(".var.2", "", dat$feature_short, fixed=T)

levels_trend <-  paste0(rep(trajectories, each=4), "_",
         rep(c("up", "interup", "interdown", "down"), length(trajectories)))
dat$trajectory_trend <- factor(paste0(dat$trajectory, "_",
     dat$general_trend), levels=levels_trend)

a <- acast(dat,  feature ~ trajectory_trend , mean, fill=1,
           drop=T, value.var="adjusted_pvalue")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
dat$feature <- factor(dat$feature, levels=levels_motifs)

g1 <- ggplot(dat, aes(y=feature, x=trajectory_trend, fill=adjusted_pvalue)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
        rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=4), 
          axis.text.x=element_text(size=5, angle=90,  vjust = 0.5, hjust=1),
          panel.background=element_rect(fill="white", colour="white")) +
    scale_x_discrete(drop=FALSE) +
    xlab("Gene Trend")

ggsave(g1, file="supp_figures/SuppFig5_Motif_enrichment_trajectories_trends_all.svg", height=15, width=6.99)
