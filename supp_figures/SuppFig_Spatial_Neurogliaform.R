library(spatialLIBD)
library(ggplot2)
library(dplyr)
library(scater)
library(hexbin)

sce <- fetch_data(type = 'sce')

genes <-  c("LAMP5", "NOS1", "ADARB1", "LHX6", "CA1", "CALB2", "NDNF", "GAD1", "GAD2", "CCK")

tmp <- sce@assays@data$logcounts[match(genes, rowData(sce)$gene_name),]
rownames(tmp) <- genes
tmp <- as.matrix(tmp)
tmp <- t(tmp)
tmp <- as.data.frame(tmp)
tmp$layer <- sce$layer_guess_reordered_short
tmp <- tmp[!tmp$layer %in% c("WM"),]
tmp <- tmp[!is.na(tmp$layer),]

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(tmp, aes(x=GAD1, y=LAMP5))  + theme_classic() +
    stat_binhex(bins=40) +
    scale_fill_gradientn(colours =
       grDevices::colorRampPalette(colors=blue_colours)(20)[-c(1:3)])
ggsave(g1, file="supp_figures/SuppFig_LAMP5_GAD1.svg", height=8, width=7)


ind <- tmp[,"GAD1"]>0 & tmp[,"LAMP5"]>0
sum(ind)/dim(tmp)[1]

tmp <- tmp[ind,]

g2a <- ggplot(tmp, aes(x=NOS1, y=NDNF)) + theme_classic() +
    stat_binhex(bins=40) +
    scale_fill_gradientn(colours =
      grDevices::colorRampPalette(colors=blue_colours)(20)[-c(1:3)])
g2b <- ggplot(tmp, aes(x=NOS1, y=CCK)) + theme_classic() +
    stat_binhex(bins=40) +
    scale_fill_gradientn(colours =
       grDevices::colorRampPalette(colors=blue_colours)(20)[-c(1:3)])
g2c <- ggplot(tmp, aes(x=NDNF, y=CCK)) + theme_classic() +
    stat_binhex(bins=40) +
    scale_fill_gradientn(colours =
     grDevices::colorRampPalette(colors=blue_colours)(20)[-c(1:3)])
ggsave(g2a, file="supp_figures/SuppFig_NOS1_NDNF.svg", height=8, width=7)
ggsave(g2b, file="supp_figures/SuppFig_NOS1_CCK.svg", height=8, width=7)
ggsave(g2c, file="supp_figures/SuppFig_NDNF_CCK.svg", height=8, width=7)

tmp$layer <- factor(tmp$layer, levels=paste0("L", 1:6))
res <- rbind(table(tmp$layer[tmp$NOS1==0 & tmp$NDNF==0 & !tmp$CCK==0]), 
             table(tmp$layer[(!tmp$NOS1==0 & tmp$NDNF==0) |
                    (!tmp$CCK==0 & !tmp$NOS1==0)]),
             table(tmp$layer[(tmp$NOS1==0 & !tmp$NDNF==0) |
                    (tmp$NOS1==0 & tmp$CCK==0)]))
rownames(res) <- c("CCK+", "NOS1+", "NDNF+")
apply(res, 1, function(x) x/colSums(res))
apply(res, 2, function(x) x/rowSums(res))
a <- chisq.test(res)
res1 <- a$observed/a$expected
rownames(res1) <- c("CCK+", "NOS1+", "NDNF+")

res1 <- reshape2::melt(res1)
colnames(res1) <- c("Cell type", "Layer", "Stat")
g1 <- ggplot(res1, aes(y=`Cell type`, x=Layer, fill=Stat)) + geom_tile() + 
    theme_classic() +  scale_fill_gradientn(colours = terrain.colors(10)[-10])
ggsave(g1, file="supp_figures/SuppFig_Location_Neurogliaform.svg", 
       width=7, height=7)
    
chisq.test(res[1:2,])
chisq.test(res[c(1,3),])
chisq.test(res[c(2,3),])

