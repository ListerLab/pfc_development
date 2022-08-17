################################################################################
#                                                                              #
#                                                                              #
# All plot for transcription factor analysis                                   #
#                                                                              #
#                                                                              #    
################################################################################

library(ggplot2) 
library(readxl)
library(reshape2)
library(dplyr)
library(cba)
library(cowplot)
library(survcomp)

source("snATACseq/R/functions_TF.R")

## 1. TF for NMF distinct peaks

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

tf_dat <- read.table("snATACseq/processed_data/TF_nmf.tsv")
tf_dat$feature_short <- sapply(tf_dat$feature, function(x)
    strsplit(x, "_")[[1]][1])
tf_dat$feature_short <- gsub("..", "+", tf_dat$feature_short, fixed=T)
tf_dat$feature_short <- gsub(".var.2", "", tf_dat$feature_short, fixed=T)

tf_dat <- attach_family_name(tf_dat, motif_cluster)

tf_dat <- remove_non_expressed_tf(tf_dat, tf_not)
tf_dat$padjust <- as.numeric(tf_dat$padjust)

levels_trend <- c(9,11,14,26,21,16,23,30,28,6, 
  12,5,10,15,17,7,18,8,4,13,25,24,22,31,20)
tf_dat$cell_type <- factor(tf_dat$cell_type, levels=levels_trend)

tf_dat_comb <- tf_dat %>% group_by(family_name, cell_type) %>% 
    summarize(p.summary=combine.test(padjust))
tf_dat_comb <- tf_dat_comb[tf_dat_comb$p.summary<0.05,]

levels_motif <- find_motif_order(tf_dat_comb)
tf_dat_comb$family_name <- factor(tf_dat_comb$family_name, levels=levels_motif)

g1 <- ggplot(tf_dat_comb, 
             aes(y=family_name, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=8, angle=0, hjust=1), 
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) +
    xlab("NMF Cluster") + background_grid()

ggsave(g1, file="supp_figures/tf_nmf_clusters.svg", height=8, width=4)

## 2. TF for CREs

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

### 2a. TF for CREs in major trends

tf_dat <- read.table("snATACseq/processed_data/gene_trends_tf.tsv", sep="\t")
tf_dat$feature_short <- sapply(tf_dat$feature, function(x)
    strsplit(x, "_")[[1]][1])
tf_dat$feature_short <- gsub("..", "+", tf_dat$feature_short, fixed=T)
tf_dat$feature_short <- gsub(".var.2", "", tf_dat$feature_short, fixed=T)

tf_dat <- attach_family_name(tf_dat, motif_cluster)

tf_dat <- remove_non_expressed_tf(tf_dat, tf_not)
tf_dat$adjusted_pvalue <- as.numeric(tf_dat$padjust)

cell_types <- paste0("G", 1:14)
tf_dat$cell_type  <- factor(tf_dat$cell_type , levels=cell_types)

tf_dat_comb <- tf_dat %>% group_by(family_name, cell_type) %>% 
    summarize(p.summary=combine.test(adjusted_pvalue))
tf_dat_comb <- tf_dat_comb[tf_dat_comb$p.summary<0.05,]

levels_motif <- find_motif_order(tf_dat_comb)
tf_dat_comb$family_name <- factor(tf_dat_comb$family_name, levels=levels_motif)

g1 <- ggplot(tf_dat_comb, 
             aes(y=family_name, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=8, angle=0, hjust=1), 
          axis.text.x=element_text(size=8, angle=90, hjust=1),
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) +
    xlab("") + background_grid()

ggsave(g1, file="main_figures/tf_gene_trends_all.svg", height=10, width=6)

rare_family <- table(tf_dat_comb$family_name)
rare_family <- names(rare_family)[rare_family<5]

g1 <- ggplot(tf_dat_comb[tf_dat_comb$family_name %in% rare_family,], 
             aes(y=family_name, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=8, angle=0, hjust=1), 
          axis.text.x=element_text(size=8, angle=90, hjust=1),
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) +
    xlab("") + background_grid()

ggsave(g1, file="main_figures/tf_gene_trends_sub.svg", height=6, width=6)

### 2b. TF for CREs in major trends per trajectory

tf_dat <- read.csv("snATACseq/processed_data/tf_general_trends_trajectories_sig.csv")
tf_dat$feature_short <- sapply(tf_dat$feature, function(x)
    strsplit(x, "_")[[1]][1])
tf_dat$feature_short <- gsub("..", "+", tf_dat$feature_short, fixed=T)
tf_dat$feature_short <- gsub(".var.2", "", tf_dat$feature_short, fixed=T)

tf_dat <- attach_family_name(tf_dat, motif_cluster)

tf_dat <- remove_non_expressed_tf_cell_type(tf_dat, tf_not_cell_type)
tf_dat$adjusted_pvalue <- as.numeric(tf_dat$adjusted_pvalue)

cell_types <- c("L2_3", "L4", "L5_6", "CGE_der", "MGE_der", "Astro", "Oligo")
general_trends <- c("up", "interup", "interdown", "down")

levels_trend <- paste0(rep(cell_types, each=length(general_trends)), "..",
        rep(general_trends, length(cell_types)))
tf_dat$cell_type <- paste0(tf_dat$trajectory, "..", tf_dat$general_trend)
tf_dat$cell_type  <- factor(tf_dat$cell_type , levels=levels_trend)

tf_dat_comb <- tf_dat %>% group_by(family_name, cell_type) %>% 
    summarize(p.summary=combine.test(adjusted_pvalue))
tf_dat_comb <- tf_dat_comb[tf_dat_comb$p.summary<0.05,]

levels_motif <- find_motif_order(tf_dat_comb)
tf_dat_comb$family_name <- factor(tf_dat_comb$family_name, levels=levels_motif)

g1 <- ggplot(tf_dat_comb, 
             aes(y=family_name, x=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=8, angle=0, hjust=1), 
          axis.text.x=element_text(size=8, angle=90, hjust=1),
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) +
    xlab("") + background_grid()

ggsave(g1, file="supp_figures/tf_trajectories.svg", height=10, width=6)

## 3. TF for ion transport

tf_dat <- read.csv("snATACseq/processed_data/tf_ion_trends.csv")
tf_dat$feature_short <- sapply(tf_dat$feature, function(x)
    strsplit(x, "_")[[1]][1])
tf_dat$feature_short <- gsub("..", "+", tf_dat$feature_short, fixed=T)
tf_dat$feature_short <- gsub(".var.2", "", tf_dat$feature_short, fixed=T)

tf_dat <- attach_family_name(tf_dat, motif_cluster)

tf_dat <- remove_non_expressed_tf(tf_dat, tf_not_cell_type)
tf_dat$adjusted_pvalue <- as.numeric(tf_dat$adjusted_pvalue)

cell_types <- c("L2_3", "L4", "L5_6", "CGE_der", "MGE_der", "Astro", "Oligo")
general_trends <- c("up", "interup", "interdown", "down")

levels_trend <- paste0(rep(cell_types, each=length(general_trends)), "..",
        rep(general_trends, length(cell_types)))
tf_dat$cell_type <- paste0(tf_dat$trajectory, "..", tf_dat$general_trend)
tf_dat$cell_type  <- factor(tf_dat$cell_type , levels=levels_trend)

tf_dat_comb <- tf_dat %>% group_by(family_name, cell_type) %>% 
    summarize(p.summary=combine.test(adjusted_pvalue))
tf_dat_comb <- tf_dat_comb[tf_dat_comb$p.summary<0.05,]

levels_motif <- find_motif_order(tf_dat_comb)
tf_dat_comb$family_name <- factor(tf_dat_comb$family_name, levels=levels_motif)

g1 <- ggplot(tf_dat_comb, 
             aes(x=family_name, y=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=8, angle=0, hjust=1), 
          axis.text.x=element_text(size=8, angle=90, hjust=1),
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) + xlab("") +
    background_grid()

ggsave(g1, file="main_figures/tf_ion_gos.svg", width=8, height=4.5)

## 4. TF for disease enrichment

tf_dat <- data.table::fread("snATACseq/processed_data/tf_diseases_new.tsv", 
    sep="\t", header=T, data.table = FALSE)
tf_dat$feature_short <- sapply(tf_dat$feature, function(x)
    strsplit(x, "_")[[1]][1])
tf_dat$feature_short <- gsub("..", "+", tf_dat$feature_short, fixed=T)
tf_dat$feature_short <- gsub(".var.2", "", tf_dat$feature_short, fixed=T)

tf_dat <- attach_family_name(tf_dat, motif_cluster)

tf_dat <- remove_non_expressed_tf(tf_dat, tf_not_cell_type)
tf_dat$adjusted_pvalue <- as.numeric(tf_dat$adjusted_pvalue)

cell_types <- c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "PV", "SST")
general_trends <- c("up", "interup", "down")
pathway <- c("Alzheimer's Disease", "Dementia", "Amnesia", "Impaired cognition",
    "Mild cognitive disorder", "Schizophrenia", "Bipolar Disorder",
    "Major Depressive Disorder", "Autism Spectrum Disorders", 
    "Attention deficit hyperactivity disorder", "Narcolepsy", 
    "Huntington Disease", "Alcoholic Intoxication, Chronic")

levels_trend <- paste0(rep(general_trends, each=length(cell_types)), "..",
        rep(cell_types, length(general_trends)))
levels_trend <- paste0(rep(pathway, each=length(levels_trend)), "..",
      rep(levels_trend, length(pathway)))
tf_dat$cell_type <- paste0(tf_dat$pathway,  "..", tf_dat$general_trend,
                           "..", tf_dat$trajectory)
tf_dat$cell_type  <- factor(tf_dat$cell_type , levels=levels_trend)

tf_dat_comb <- tf_dat %>% group_by(family_name, cell_type) %>% 
    summarize(p.summary=combine.test(adjusted_pvalue))
tf_dat_comb <- tf_dat_comb[tf_dat_comb$p.summary<0.05,]

levels_motif <- find_motif_order(tf_dat_comb)
tf_dat_comb$family_name <- factor(tf_dat_comb$family_name, levels=levels_motif)

g1 <- ggplot(tf_dat_comb, 
             aes(x=family_name, y=cell_type, fill=p.summary)) + 
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
       rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.y=element_text(size=7, angle=0, hjust=1), 
          axis.text.x=element_text(size=8, angle=90, hjust=1),
          panel.background=element_rect(fill="white", colour="white"),
          axis.title.y=element_blank()) + xlab("") +
    background_grid()

ggsave(g1, file="main_figures/tf_disease_new.svg", height=7, width=6)



