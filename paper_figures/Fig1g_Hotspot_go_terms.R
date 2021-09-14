# Hotspot dotplot

library(ggplot2)
library(reshape2)
library(ggsci)
library(rcartocolor)
library(cba)
library(palettetown)
library(dplyr)

dat <- read.csv("snRNAseq/processed_data/hotspot_GOs-vs-modules_dataframe.csv", 
                check.names = FALSE)
dat_inter <- read.csv("snRNAseq/processed_data/hotspot_GOs-vs-modules_dataframe_intersect.csv", 
                      check.names = FALSE)
rownames(dat) <- dat[,1]
colnames(dat)[1] <- "modules"
rownames(dat_inter) <- dat_inter[,1]
colnames(dat_inter)[1] <- "modules"

colnames(dat) <- gsub("extracellular matrix", "ECM",colnames(dat))
colnames(dat) <- gsub("transmembrane transporter activity", "TTA", colnames(dat))
colnames(dat) <- gsub("G-protein coupled receptor", "GPCR", colnames(dat))
colnames(dat) <- gsub("plasma membrane", "PM", colnames(dat))

order_row <- c(14, 4, 11, 6, 7, 12, 3, 5, 1, 10, 2, 9, 8, 13)
hc <- hclust(dist(t(dat[order_row,-1])))
order_col <- order.optimal(dist(t(dat[order_row,-1])), hc$merge)


dat <- melt(dat, id.vars="modules")
dat_inter <- melt(dat_inter, id.vars="modules")
colnames(dat) <- c("modules", "GO terms", "p-values")
dat$intersect <- dat_inter[,3]
dat[dat[,3]==0,3] <- NA
dat[dat[,4]==0,4] <- NA
dat$modules <- sapply(dat$modules, function(x) strsplit(x, "-")[[1]][2])

order_modules <- c(1:14)[order_row]
dat$modules <- factor(dat$modules, levels=rev(order_modules))

order_go_terms <- names(order_col$order)[order_col$order]
dat$`GO terms` <- factor(dat$`GO terms`, levels=order_go_terms)


dat <- dat[!is.na(dat[,3]),]
dat <- dat[dat$`p-values`>log10(0.05)*-1,]

gg_final <- ggplot(dat, aes(y=modules, x=`GO terms`, 
    fill=`p-values`, size=intersect)) + 
    geom_point(colour="black", pch=21) + theme_light() +
    scale_fill_material("grey") + 
    scale_x_discrete(guide = guide_axis(angle = -45)) +
    scale_y_discrete(expand=expansion(mult = 0, add = 0.75)) +
    scale_size(range = c(0, 5)) +
    theme(plot.margin=unit(c(5.5, 40, 5.5, 5.5), "points"),
          axis.text.x = element_text(size=7)) 

  
ggsave(gg_final, file="paper_figures/Fig1g_Hotspot_go_terms.svg", height=3.5, width=11)
