# PCA plots

library(bsseq)
library(GenomicRanges)
library(stringr)
library(magrittr)
library(data.table)
library(GenomicFeatures)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(rtracklayer)
library(limma)
library(edgeR)
library(parallel)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

peak_mat_list <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
mdat <- readxl::read_xlsx("snATACseq/annotation/scATACseq_neuronal_maturation.xlsx")
mdat <- mdat[mdat$Used=="Yes",]


col <- c(Astro= '#ffc857', 
         `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
         Micro='#484848', Oligo='#255f85', 
         `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')

plot_pca <- function(cell_types, title, fpkm=TRUE){
    
    cell_dat <- lapply(peak_mat_list[cell_types], function(x)
        x@elementMetadata)
    
    # Extract data
    cell_dat <- lapply(cell_dat, function(x)
        x[,grepl(paste0("_RL"), colnames(x))])
    cell_dat <- lapply(cell_dat, function(x)  x%>% as.data.frame() %>%
        limma::normalizeBetweenArrays(method = "quantile"))
    for(i in 1:length(cell_dat)){
        
        colnames(cell_dat[[i]]) <- gsub(paste0(names(cell_dat)[i], "_"), "", 
            colnames(cell_dat[[i]]))
    }
    
    pr <- lapply(cell_dat, function(x) prcomp(t(x), scale.=TRUE))
    pc1 <- lapply(pr, function(x) 
        (summary(x)$importance[2, 1] * 100) %>% round(digits = 2))
    tmp <- unlist(lapply(1:length(cell_types), function(xx) rep(cell_types[xx],
        length=ncol(cell_dat[[xx]]))))
    
    pca_df <- data.frame(PC1=unlist(lapply(pr, function(xx) xx$x[ ,1])),
        Cell_type=tmp,
        Sample=unlist(lapply(cell_dat, function(x) colnames(x))))

    ind <- match(pca_df$Sample, mdat$Sample)
    stage <- mdat$Stage[ind] 
    arc_age <- mdat1[ind,2]
    
    pca_df$Stage <- stage
    pca_df$Age <- arc_age
    pca_df$Stage <- factor(pca_df$Stage, levels=c("Fetal", "Neonatal",
        "Infancy", "Childhood", "Adolescence", "Adult"))
    
    pca_df_mean <- pca_df %>% group_by(Stage, Cell_type) %>% 
        summarise(PC1=mean(PC1))
    
    pc1_per <- sapply(1:length(pc1), function(x)
        paste0(names(pc1)[x], ":", pc1[x], "%"))
    pc1_per <- paste0(pc1_per, collapse=", ")
    
    gg_pca <- ggplot(data = pca_df, mapping = aes(x = Stage, y = PC1,
        col=Cell_type)) + geom_point() +
        geom_line(dat=pca_df_mean, aes(x = Stage, y = PC1,
           group=Cell_type), alpha = 1, size = 2) +
        ggtitle(title) + xlab("") +
        ylab(str_c("PC1 (", pc1_per, ")")) + theme_classic() +
        scale_color_manual(values=col) + theme(axis.title.y = element_text(
            size=10))
    
    return(gg_pca)
    
}


g1 <- plot_pca(c("Astro", "Oligo", "Micro"), "Glia")
ggsave(g1, file="paper_figures/Fig3f_Glia_PC1.svg", height=3, width=7)
g2 <- plot_pca(c("L2_3", "L4", "L5_6"), "Principal Neurons")
ggsave(g2, file="paper_figures/Fig3f_PN_PC1.svg", height=3, width=7)
g3 <- plot_pca(c("CGE_der", "MGE_der"), "Inhibitory Neurons")
ggsave(g3, file="paper_figures/Fig3f_IN_PC1.svg", height=3, width=7)
