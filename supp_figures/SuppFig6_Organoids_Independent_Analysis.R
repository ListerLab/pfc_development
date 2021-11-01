library(scran)
library(scater)
library(scuttle)
library(batchelor)
library(cowplot)

organoids <- c("snRNAseq/processed_data/orgpred/Luc9228.RDS",
               "snRNAseq/processed_data/orgpred/org2290.RDS",
               "snRNAseq/processed_data/orgpred/org2432.RDS")
samples <- c("LucH9228", "RL2290", "RL2432")

# Downsampling and individual clustering approach
sce <- lapply(organoids, function(x) readRDS(x))
dec <- list()
for(i in 1:length(sce)){
    
    bycol = 1000/colSums(counts(sce[[i]]))            #ensures downsample target is 1000
    set.seed(100)
    downsampled = downsampleMatrix(x = counts(sce[[i]]), 
        prop = bycol, bycol = TRUE)
    counts(sce[[i]]) <- downsampled
    sce[[i]]$Sample <- samples[i]
    sce[[i]] <- logNormCounts(sce[[i]])
    dec[[i]] <- modelGeneVar(sce[[i]])
    
    
}

org_sce <- readRDS(
    "snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

hvg<- lapply(dec, function(x) x$bio>0)
for(i in 1:length(sce)){
    sce[[i]] <- runPCA(sce[[i]], subset_row=hvg[[i]],
                      BSPARAM=BiocSingular::RandomParam())
    sce[[i]] <- runUMAP(sce[[i]], dimred="PCA")

}

for(i in 1:length(sce)){
    df <- reducedDim(sce[[i]], "UMAP")
    colnames(df) <- c("UMAP1", "UMAP2")
    df <- as.data.frame(df)
    df$cell_type <- sce[[i]]$predcelltype
    df$stage <- sce[[i]]$predstage
    df$batch <- sce[[i]]$batch
    set.seed(30)
    df <- df[sample(1:dim(df)[1]),]
    
    col <- c(Astro= '#ffc857', CGE_dev='#c6d5c0', ID2= '#558140', 
             `L2/3_CUX2`= '#6e0614', L4_RORB= '#8b3843', `L5/6_THEMIS_TLE4`='#a86a72',
             `L5/6_TLE4`='#c59ba1', LAMP5_CA1='#8eab80', MGE_dev='#ecd1c8',
             Micro='#484848',OPC='#92afc2', Oligo='#255f85', PN_dev='#e2cdd0',
             PV='#c77459', PV_SCUBE3='#daa290',`Poor-Quality`='lightgrey',
             SST='#b44622',VIP='#1c5701', Vas='#a3a3a3')
    col_stage <- c(Fetal="#512568", Neonatal="#443682", Infancy="#3D6B93",
                   Childhood="#20988C", Adolescence="#98CA43", Adult="#F9E51B")
    
    g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=cell_type)) + 
        geom_point(size=0.5, alpha=0.5) + theme_classic() +
        scale_color_manual(values=col) + 
        guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))
    
    g2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=stage)) + 
        geom_point(size=0.5, alpha=0.5) + theme_classic() +
        scale_color_manual(values=col_stage) +
        guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))
    
    p1 <- plot_grid(g1, g2)
    ggsave(paste0("supp_figures/SuppFig7_Organoids_Single_",
        samples[i], ".png"), p1,  height=5, width=10, device="png")
    
}
    
    
    
# Combined analysis     
sce <- lapply(organoids, function(x) readRDS(x))
dec <- list()
for(i in 1:length(sce)){
    
    sce[[i]]$Sample <- samples[i]
    sce[[i]] <- logNormCounts(sce[[i]])
    dec[[i]] <- modelGeneVar(sce[[i]])
    
}

rescaled <- multiBatchNorm(sce[[1]], sce[[2]], sce[[3]])
combined.dec <- combineVar(dec[[1]], dec[[2]], dec[[3]])
chosen.hvgs <- combined.dec$bio > 0
rescaled <- do.call(cbind, rescaled)
rescaled$Batch <- rescaled$Sample

set.seed(0010101010)
uncorrected <- runPCA(rescaled, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())
uncorrected <- runUMAP(uncorrected, dimred="PCA")
corrected <- rescaleBatches(rescaled, batch=rescaled$Batch)
corrected<- runPCA(corrected, subset_row=chosen.hvgs, 
                   exprs_values="corrected",
                   BSPARAM=BiocSingular::RandomParam())
corrected <- runUMAP(corrected, dimred="PCA", spread=0.2, min_dist=0.075)
corrected$predcelltype <- uncorrected$predcelltype
corrected$predstage <- uncorrected$predstage

df <- reducedDim(corrected, "UMAP")
colnames(df) <- c("UMAP1", "UMAP2")
df <- as.data.frame(df)
df$cell_type <- corrected$predcelltype
df$stage <- corrected$predstage
df$batch <- corrected$batch
set.seed(30)
df <- df[sample(1:dim(df)[1]),]
df$stage <- factor(df$stage, levels=c("Fetal", "Neonatal", "Infancy",
    "Childhood", "Adolescence", "Adult'"))

col <- c(Astro= '#ffc857', CGE_dev='#c6d5c0', ID2= '#558140', 
         `L2/3_CUX2`= '#6e0614', L4_RORB= '#8b3843', `L5/6_THEMIS_TLE4`='#a86a72',
         `L5/6_TLE4`='#c59ba1', LAMP5_CA1='#8eab80', MGE_dev='#ecd1c8',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', PN_dev='#e2cdd0',
         PV='#c77459', PV_SCUBE3='#daa290',`Poor-Quality`='lightgrey',
         SST='#b44622',VIP='#1c5701', Vas='#a3a3a3')
col_stage <- c(Fetal="#512568", Neonatal="#443682", Infancy="#3D6B93",
               Childhood="#20988C", Adolescence="#98CA43", Adult="#F9E51B")

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=cell_type)) + 
    geom_point(size=0.5, alpha=0.5) + theme_classic() +
    scale_color_manual(values=col) + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

g2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=stage)) + 
    geom_point(size=0.5, alpha=0.5) + theme_classic() +
    scale_color_manual(values=col_stage) + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

g3 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=batch)) + 
    geom_point(size=0.5, alpha=0.5) + theme_classic()  + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

p1 <- plot_grid(g1,g2,g3)
ggsave(p1, file="supp_figures/SuppFig7_Organoid_Combined.png", height=7, width=7)

