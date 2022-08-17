library(scran)
library(scater)
library(scuttle)
library(batchelor)
library(cowplot)
library(ggpubr)
library(Seurat)
library(viridis)
library(ggplot.multistats)
library(reticulate)
library(Seurat)
library(uwot)
library(FNN)
library(stringr)
library(ggridges)
library(forcats)
library(ComplexHeatmap)
library(bluster)
library(readxl)


markers <- read.table("snRNAseq//processed_data/Organoids_Brain_markers_PN.csv",
      sep=",")
markers_up <- rownames(markers)[markers$FDR<0.05 & 
        markers$summary.logFC>0]
markers_down <- rownames(markers)[markers$FDR<0.05 & 
        markers$summary.logFC<0]

motif_anno <- read_xlsx("annotation/motif_annotations.xlsx", sheet=2)
motif_anno <- motif_anno[motif_anno$Database=="Jaspar2018",]
motif_anno$motif <- sapply(motif_anno$Motif, function(x)
  strsplit(x, "_")[[1]][1])
motifs <- toupper(motif_anno$motif)
motifs <- unlist(sapply(motifs, function(x) strsplit(x, "+", fixed=TRUE)[[1]]))
motifs <- unique(motifs)
motifs <- c(motifs, c("EWSR1","FLI1")) 

all_networks <- readRDS("snATACseq/processed_data/all_networks.rds")
all_networks <- all_networks[2:4]

markers_tf_down <- intersect(markers_down, motifs)
markers_tf_up <- intersect(markers_up, motifs)

markers_targets_down <- setdiff(markers_down, motifs)
markers_targets_up <- setdiff(markers_up, motifs)

all_networks_down <- lapply(all_networks, function(x) 
  x[x$CRE %in% markers_targets_down,])
all_networks_up <- lapply(all_networks, function(x) 
  x[x$CRE %in% markers_targets_up,])

count_cres_motif <- function(markers_tf, all_networks){
  
  mean(sapply(all_networks, function(x) sum(x$motif==markers_tf)))
   
}


tf_down <- sapply(markers_tf_down, function(x) 
  count_cres_motif(x, all_networks_down))
tf_down <- (tf_down-mean(tf_down))/sd(tf_down)
tf_down_sel<- sort(tf_down, decreasing=TRUE)[1:10]
tf_up <- sapply(markers_tf_up, function(x) 
  count_cres_motif(x, all_networks_up))
tf_up <- (tf_up-mean(tf_up))/sd(tf_up)
tf_up_sel<- sort(tf_up, decreasing=TRUE)[1:3]

dat <- data.frame(motif=names(tf_down_sel), zscore=tf_down_sel)
dat$motif <- factor(dat$motif, level=rev(names(tf_down_sel)))
dat$motif_family <- factor(c("MZF1  MZF1", "RARA  NR/3", "JUN  CREB/ATF/2", 
          "JUND  AP1/1", "RORA   NR/3", "SREBF2  Ebox/CACGTG/1",
          "JDP2  CREB/ATF/2", "ZNF263  GC-tract", "ZNF384  ZNF384/2",
          "ZNF354C  ZNF53"), levels=rev(c("MZF1  MZF1", "RARA  NR/3", "JUN  CREB/ATF/2", 
          "JUND  AP1/1", "RORA   NR/3", "SREBF2  Ebox/CACGTG/1",
          "JDP2  CREB/ATF/2", "ZNF263  GC-tract", "ZNF384  ZNF384/2",
          "ZNF354C  ZNF53")))
dat$family <- c("bZIP", "bZIP", "bZIP", "bZIP", "bZIP", "bHLH", "bZIP",
                "bHLH", "C2H2", "C2H2")

gene_trends <- read.csv("snRNAseq//processed_data/gene_cluster_ids.csv")
gene_trends[gene_trends$gene_name %in% dat$motif,]

g1 <- ggplot(dat, aes(y=motif_family, x=zscore,
          fill=family)) + geom_col() + theme_classic() + 
          scale_fill_manual(values=c(bHLH="#429488", bZIP="#75FB4C", 
          C2H2="#EA33F7"))

ggsave(plot=g1, 'main_figures/tf_brain.svg', width=4, height=4)
