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

motif_anno <- read_xlsx("annotation/motif_annotations.xlsx", sheet=2)
motif_anno <- motif_anno[motif_anno$Database=="Jaspar2018",]
motif_anno$Motif_short <- sapply(motif_anno$Motif, function(x)
    strsplit(x, "_")[[1]][1])
motif_anno$Motif_short <- toupper(motif_anno$Motif_short)
tf_not <- readRDS("processed_data/tf_not_rna.RDS")
tf_not <- gsub("::", "+", tf_not)

motif_cluster <- read_xlsx("annotation/motif_annotations.xlsx", sheet=1)

all_tfs <- read.table("processed_data/tf_gene_trends.tsv", sep="\t")
all_tfs$feature_short <- sapply(all_tfs$feature, function(x)
    strsplit(x, "_")[[1]][1])
all_tfs$feature_short <- gsub("..", "+", all_tfs$feature_short, fixed=T)
all_tfs$feature_short <- gsub(".var.2", "", all_tfs$feature_short, fixed=T)

index <- lapply(all_tfs$feature_short, function(x) 
    which(motif_anno$Motif_short %in% x))
all_tfs <- all_tfs[!lengths(index)==0,]
index <- index[!lengths(index)==0]

all_tfs_new <- list()

for(i in 1:length(index)){
    
    tmp <- all_tfs[i,]
    tmp <- matrix(as.character(tmp), nrow=length(index[[i]]), 
                  ncol=dim(all_tfs)[2] ,byrow = TRUE)
    tmp <- as.data.frame(tmp)
    tmp$V15 <- motif_cluster$Name[match(
        motif_anno$Cluster_ID[index[[i]]], motif_cluster$Cluster_ID)]
    all_tfs_new[[i]] <- tmp
    
}

all_tfs_new <- do.call(rbind, all_tfs_new)

colnames(all_tfs_new) <- c(colnames(all_tfs), "motif_cluster")     
ind <- !(all_tfs_new$feature_short %in% tf_not)
all_tfs_new <- all_tfs_new[ind,]
all_tfs_new$padjust <- as.numeric(all_tfs_new$padjust)

levels_trend <- c(7, 0, 3, 6, 11, 5, 13, 12, 8, 9, 1, 2, 10, 4)
all_tfs_new$cell_type <- factor(all_tfs_new$cell_type, levels=levels_trend) 

tmp <- all_tfs_new %>% group_by(motif_cluster, cell_type) %>% 
    summarize(p.summary=combine.test(padjust))
tmp <- tmp[tmp$p.summary<0.05,]

a <- acast(tmp,  motif_cluster ~ cell_type , mean, fill=1,
           drop=T, value.var="p.summary")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
tmp$motif_cluster <- factor(tmp$motif_cluster, levels=levels_motifs)

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")
tmp <- tmp[tmp$p.summary<0.05,]

g1 <- ggplot(tmp, aes(y=motif_cluster, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
     rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=5), 
          panel.background=element_rect(fill="white", colour="white")) +
    xlab("Gene Trend")
ggsave(g1, file="supp_figures/SuppFig4_enrichment_gene_trends.svg", height=6, width=4)

g1 <- all_tfs %>% filter(padjust<0.05) %>%
    ggplot(aes(y=feature_short, x=cell_type, fill=padjust)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
            rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=4.5), 
          panel.background=element_rect(fill="white", colour="white")) +
    ylab("Motif") +
    xlab("Gene Trend")
ggsave(g1, file="supp_figures/SuppFig4_enrichment_gene_trends_full.svg", height=12, width=4)

subset_motifs <- tmp %>% group_by(motif_cluster, cell_type) %>% 
    group_by(motif_cluster) %>%
    summarise(n = n()) %>% filter(n<5) %>% ungroup() %>% 
    dplyr::pull(motif_cluster) 

tmp_sub <- tmp[tmp$motif_cluster %in% subset_motifs,]
all_tfs_sub <- all_tfs_new[all_tfs_new$motif_cluster %in% subset_motifs,]
all_tfs_sub <- all_tfs[all_tfs$padjust<0.05,]
saveRDS(all_tfs_sub, file="processed_data/tf_gene_trends_sub.RDS")

a <- acast(tmp_sub,  motif_cluster ~ cell_type , mean, fill=1,
           drop=T, value.var="p.summary")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
tmp_sub$motif_cluster <- factor(tmp_sub$motif_cluster, 
    levels=levels_motifs)

subset_motifs <- all_tfs %>% filter(padjust<0.05) %>%
    group_by(feature_short, cell_type) %>% group_by(feature_short) %>%
    summarise(n = n()) %>% filter(n<5) %>% ungroup() %>% 
    dplyr::pull(feature_short) 

g1 <- all_tfs[all_tfs$feature_short %in% subset_motifs,] %>% 
    filter(padjust<0.05) %>%
    ggplot(aes(y=feature_short, x=cell_type, fill=padjust)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
      rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=4.5), 
          panel.background=element_rect(fill="white", colour="white")) +
    ylab("Motif") +
    xlab("Gene Trend")
ggsave(g1, file="supp_figures/SuppFig4_enrichment_gene_trends_subset_full.svg", height=12, width=4)

