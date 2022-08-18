#organoid + glioblastoma umaps

library(BiocManager)
library(scater)
library(scran)
library(scuttle)
library(reticulate)
library(Seurat)
library(uwot)
library(FNN)
library(stringr)
library(ggridges)
library(forcats)

reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS')
initembedding = reducedDim(reference, "UMAP")

latent = readRDS('snRNAseq/processed_data/rscvi_downsampled_ref.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
    init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
    min_dist = 0.5, verbose = TRUE, n_epochs = 0)
refmap = as.data.frame(reference_umap$embedding)
refmap$dataset = rep('reference', 154748)

sce = list(Luc9228 = 'Luc9228', 
           RL2290 = 'RL2290', 
           RL2432 = 'RL2432', 
           GBM = 'GBM')

allqumaps = list()

for(i in names(sce)){
    
  qlatent = readRDS(paste0('snRNAseq/processed_data/', i, '_qlatent.RDS'))
  qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)
  qumap = as.data.frame(qumap)
  qumap$dataset= rep(paste0(i), nrow(qumap))
  qumap$pred = knn.reg(train = refmap[,1:2], test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
  qumap$cluster = knn(train = refmap[,1:2], test = qumap[,1:2], cl = reference$cluster, k = 10)
  qumap$stage = knn(train = refmap[,1:2], test = qumap[,1:2], cl = reference$stage_ids, k =10)
  g1 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
    geom_point(col='black', size=1.5) +
    geom_point(col='white', size=0.8) +
    theme_classic() +
    geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
    geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                    adjust=1/2)
  
  ggsave(plot=g1, paste0('paper_figures/', i, '_umap.png'))
  
}


