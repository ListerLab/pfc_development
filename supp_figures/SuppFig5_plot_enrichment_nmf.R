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

tf_not <- readRDS("processed_data/tf_not_rna.RDS")
tf_not <- gsub("::", "+", tf_not)

all_tfs <- read.table("processed_data/TF_nmf.tsv", sep="\t")
all_tfs$feature_short <- sapply(all_tfs$feature, function(x)
    strsplit(x, "_")[[1]][1])
all_tfs$feature_short <- gsub("..", "+", all_tfs$feature_short, fixed=T)
all_tfs$feature_short <- gsub(".var.2", "", all_tfs$feature_short, fixed=T)
ind <- !(all_tfs$feature_short %in% tf_not)
all_tfs <- all_tfs[ind,]

all_tfs <- all_tfs[all_tfs$padjust<0.05,]

a <- acast(all_tfs,  feature_short ~ cell_type , mean, fill=1,
           drop=T, value.var="padjust")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
hc <- hclust(dist(t(a[order_row$order,])))
order_col <- order.optimal(dist(t(a[order_row$order,])), hc$merge)
levels_modules <- names(order_col$order)[order_col$order]
all_tfs$feature_short <- factor(all_tfs$feature_short, levels=levels_motifs)
all_tfs$cell_type <- factor(all_tfs$cell_type, levels=levels_modules)

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(all_tfs, aes(y=feature_short, x=cell_type, fill=padjust)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
      rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=5), 
          panel.background=element_rect(fill="white", colour="white")) +
    xlab("Module") 
ggsave(g1, file="supp_figures/SuppFig5_enrichment_nmf_distinct.svg", height=6, width=4)

all_tfs <- read.table("processed_data/TF_nmf_max_mod.tsv", sep="\t",
    stringsAsFactors = FALSE)
all_tfs$feature_short <- sapply(all_tfs$feature, function(x)
    strsplit(x, "_")[[1]][1])
all_tfs$feature_short <- gsub("..", "+", all_tfs$feature_short, fixed=T)
all_tfs$feature_short <- gsub(".var.2", "", all_tfs$feature_short, fixed=T)

motif_anno <- read_xlsx("annotation/motif_annotations.xlsx", sheet=2)
motif_anno <- motif_anno[motif_anno$Database=="Jaspar2018",]
motif_anno$Motif_short <- sapply(motif_anno$Motif, function(x)
    strsplit(x, "_")[[1]][1])
motif_anno$Motif_short <- toupper(motif_anno$Motif_short)
tf_not <- readRDS("processed_data/tf_not_rna.RDS")
tf_not <- gsub("::", "+", tf_not)

motif_cluster <- read_xlsx("annotation/motif_annotations.xlsx", sheet=1)

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

tmp <- all_tfs_new %>% group_by(motif_cluster, cell_type) %>% 
    summarize(p.summary=combine.test(padjust))
tmp <- tmp[tmp$p.summary<0.05,]

a <- acast(tmp,  motif_cluster ~ cell_type , mean, fill=1,
           drop=T, value.var="p.summary")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
hc <- hclust(dist(t(a[order_row$order,])))
order_col <- order.optimal(dist(t(a[order_row$order,])), hc$merge)
levels_modules <- names(order_col$order)[order_col$order]
tmp$motif_cluster <- factor(tmp$motif_cluster, levels=levels_motifs)
tmp$cell_type <- factor(tmp$cell_type, levels=levels_modules)

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(tmp, aes(y=motif_cluster, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=5), 
          panel.background=element_rect(fill="white", colour="white")) +
    xlab("Module") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(g1, file="supp_figures/SuppFig5_enrichment_nmf_all.svg", height=12, width=8)
