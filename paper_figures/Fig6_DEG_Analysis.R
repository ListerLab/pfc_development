library(scater)
library(scran)
library(SingleCellExperiment)
library(batchelor)
library(FNN)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
library(ggplot2)

source("snATACseq/R/functions_nmf.R")

# Volcano Plot

markers_PN <- read.table(file="snATACseq/processed_data/Organoids_PN_markers.csv",
             sep=",", header = T)

markers_PN$p.adjust <- -log10(markers_PN$FDR)
g1 <- ggplot(markers_PN, aes(x=summary.logFC, y=p.adjust)) + 
    geom_point(col="black", alpha=0.5) +
    geom_point(data=markers_PN[markers_PN$summary.logFC > 0 &
      markers_PN$FDR < 0.05,], col="#990099", alpha=0.5) +
    geom_point(data=markers_PN[markers_PN$summary.logFC < 0 &
        markers_PN$FDR < 0.05,], col="#009900", alpha=0.5) +
    theme_classic() #+
    #geom_text_repel(data=markers_PN[markers_PN$summary.logFC > 0.75 &
    #    markers_PN$FDR < 0.05,], col="#100010", aes(label=Symbol),
    #    size=2, max.overlaps = 20) +
    #geom_text_repel(data=markers_PN[markers_PN$summary.logFC < -0.75 &
    #   markers_PN$FDR < 0.05,], col="#001000", aes(label=Symbol),
    #   size=2, max.overlaps = 20)

ggsave(g1, file="paper_figures/Fig6_Volcano_PN.png", width=7,
       height=7)

# GO term analysis and plots

go_pn_up <- go_term_analysis(markers_PN$ID[markers_PN$FDR<0.05 & 
        markers_PN$summary.logFC>0], markers_PN$ID)
go_pn_down <- go_term_analysis(markers_PN$ID[markers_PN$FDR<0.05 & 
        markers_PN$summary.logFC<0], markers_PN$ID)
go_pn <- go_term_analysis(markers_PN$ID[markers_PN$FDR<0.05], markers_PN$ID)

write.table(go_pn_up@result[,-8], file="snATACseq/processed_data/GO_PN_up.tsv",
            quote=F, sep="\t")
write.table(go_pn@result[,-8], file="snATACseq/processed_data/GO_PN.tsv",
            quote=F, sep="\t")
write.table(go_pn_down@result[,-8], file="snATACseq/processed_data/GO_PN_down.tsv",
            quote=F, sep="\t")

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- barplot(go_pn, showCategory=20) + scale_fill_gradientn(limits=c(0, 0.05), 
   colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20)))
ggsave(g1, file="paper_figures/Fig6_Go_PN.svg", width=7, height=7)

g1 <- barplot(go_pn_up, showCategory=20) + scale_fill_gradientn(limits=c(0, 0.05), 
   colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20)))
ggsave(g1, file="paper_figures/Fig6_Go_PN_up.svg", width=7, height=7)

g1 <- barplot(go_pn_down, showCategory=20) + scale_fill_gradientn(limits=c(0, 0.05), 
    colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20)))
ggsave(g1, file="paper_figures/Fig6_Go_PN_down.svg", width=7, height=7)

# UMAP plots

path1 <- "path_to_umaps"
all_samples <- c("Luc9228.RDS", "org2290.RDS", "org2432.RDS")

atlas <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")
colData(atlas) <- DataFrame(Sample=atlas$batch, 
                            Barcode=colnames(atlas), predcelltype=atlas$major_clust, 
                            predarcsin=atlas$arcsin_ages, predstage=atlas$stage_ids)

tmp <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")
umap_org <- as.data.frame(reducedDim(tmp, "UMAP"))
umap_org$major_clust <- tmp$major_clust
umap_org$stage <- tmp$stage_ids

rm(tmp)

all_sce <- lapply(all_samples, function(x) readRDS(paste0(path1, x)))

mature_PN <- lapply(all_sce, function(x) colnames(x)[which(
    colData(x)[, "predcelltype"] 
    %in% c("L2/3_CUX2", "L4_RORB", "L5/6_THEMIS_TLE4", "PN_dev")
    & colData(x)[, "predstage"] != "Fetal")])

path2 <- "path_to_umaps_organoids"
umaps <- lapply(all_samples, function(x) 
    readRDS(paste0(path2, gsub(".RDS", "_umap.RDS", x))))

universe <- intersect(rowData(atlas)$gene_ids, rowData(all_sce[[1]])$ID)
atlas <- atlas[match(universe, rowData(atlas)$gene_ids),]
rowData(atlas) <- rowData(all_sce[[1]])[
    match(universe, rowData(all_sce[[1]])$ID),]
rownames(atlas) <- rowData(atlas)$ID

neighbours_PN <- list()

umaps_PN <- list()

umaps_PN_im <- list()

immature_PN <- lapply(all_sce, function(x) colnames(x)[which(
    colData(x)[, "predcelltype"] 
    %in% c("L2/3_CUX2", "L4_RORB", "L5/6_THEMIS_TLE4", "PN_dev")
    & colData(x)[, "predstage"] == "Fetal")])

for(i in 1:length(umaps)){
    
    umaps_PN[[i]] <- umaps[[i]][match(mature_PN[[i]], umaps[[i]]$barcode),]
    umaps_PN_im[[i]] <- umaps[[i]][match(immature_PN[[i]], umaps[[i]]$barcode),]
    PN_id <- umap_org$major_clust %in%  
        c("L2/3_CUX2", "L4_RORB", "L5/6_THEMIS_TLE4", "PN_dev")
    
}

umaps_PN <- do.call(rbind, umaps_PN)
umaps_PN_im <- do.call(rbind, umaps_PN_im)

colnames(umap_org)[1:2]<- c("UMAP1", "UMAP2")

g1 <- ggplot(umap_org[,1:2], aes(x=UMAP1, y=UMAP2)) +
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) +
    theme_classic() +
    geom_point(data=umaps_PN_im[,1:2], alpha=0.5, size=1, col="#009900") +
    geom_point(data=umaps_PN[,1:2], alpha=0.5, size=1, col="#990099")
ggsave(g1, file="paper_figures/Fig6_UMAP_PN.png", width=7, height=7)     



