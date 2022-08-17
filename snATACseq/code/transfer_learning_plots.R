################################################################################
#                                                                              #
#                                                                              #
# Organoid plots                                                               #
#                                                                              #
#                                                                              #    
################################################################################

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
library(ggrepel)

# 1. Transfer learning plots

### 1a. Samples used to Verify Age Predictions

# add Velmeshev

reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS')
initembedding = reducedDim(reference, "UMAP")

#velmeshev plotting

query = readRDS('snRNAseq/processed_data/transfer_learning/datasets/velmeshev_dataset.RDS')
# this data can be downloaded from https://www.science.org/doi/abs/10.1126/science.aav8130

qumap <- readRDS("snRNAseq/processed_data/transfer_learning/datasets/velmeshev_umap.RDS")
refmap <- initembedding

qumap = as.data.frame(qumap)
refmap$dataset = character(154748)
refmap$dataset[1:154748] = 'reference'
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'velmeshev'

myknn = function(ref_embed, ref_class, test_embed, test_class, k){
  
  output = list()
  kest = knn(train = ref_embed, test = test_embed, 
        cl = as.factor(ref_class), prob=TRUE, k =k)
  result = character(nrow(test_embed))
  
  for(i in 1:nrow(test_embed)){
    if(kest[i] == as.factor(test_class)[i]) {     
      result[i] = "correct"
    }
    else {
      result[i] = "incorrect"
    }
  }
  
  output$result = result
  correct = subset(test_class, result == "correct")
  incorrect = subset(test_class, result == "incorrect")
  output$correct = correct
  output$incorrect = incorrect
  output$accuracy = length(correct)/length(result)
  output$corprop = summary(correct)/length(correct)
  output$est = kest
  output$percents = summary(correct)/summary(test_class)

  return(output)
  
}


qumap$pred = knn.reg(train = refmap[,1:2], test = qumap[,1:2], 
    y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = query$arcsin_ages
qumap$diff = abs(qumap$pred - qumap$age)

#add samples - and remove non-specific samples
qumap$sample = query$sample
qumap$sample = fct_reorder(qumap$sample, qumap$age)
qumap = subset(qumap, !(sample %in% 
    c('5893_PFC', '5538_PFC_Nova', '5879_PFC_Nova','5936_PFC_Nova', 
      '5408_PFC_Nova', '5144_PFC', '5278_PFC', '5403_PFC', 
      '5419_PFC', '5945_PFC')))
smolquery = query[, !(query$sample %in% c('5893_PFC', '5538_PFC_Nova', 
    '5879_PFC_Nova','5936_PFC_Nova', '5408_PFC_Nova', '5144_PFC', 
    '5278_PFC', '5403_PFC', '5419_PFC', '5945_PFC'))]
out = myknn(ref_embed = refmap[,1:2], test_embed = qumap[,1:2], 
    test_class = smolquery$cluster, ref_class = reference$cluster, k =10)
qumap$result = out$result
qumap$truecluster = smolquery$cluster

g1 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(plot=g1, 'supp_figures/krieg_umap.png')

g2 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap[sample(1:nrow(qumap)),], alpha=0.5, size=1, aes(col=truecluster))  +
  scale_colour_manual(values = c('Astro' ='#FFC857','CGE_dev'='#C6D5C0',
        'ID2'='#558140','L2/3_CUX2'='#6E0614',
         'L4_RORB' ='#8B3843', 'L5/6_THEMIS_TLE4' = '#A86A72',
         'MGE_dev' ='#ECD1C8','Micro' ='#484848',
         'Oligo'= '#255F85','OPC' ='#92AFC2',
         'PN_dev' ='#E2CDD0', 'Poor-Quality' ='#FFFFFF',
         'PV'= '#C77459', 'SST' ='#B44622',
         'VIP'= '#1C5701', 'Vas' ='#A3A3A3', 
        'Inhib' ='pink', 'Mix' = 'light blue'))

ggsave(plot=g2, 'supp_figures/krieg_umap_clusters.png')

allmap = qumap

#add ga39 sample

latent = readRDS('snRNAseq/processed_data/transfer_learning/outs/4999_latent.RDS')
qlatent = readRDS('snRNAseq/processed_data/transfer_learning/outs/ga39_qlatent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
    init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
    min_dist = 0.5, verbose = TRUE, n_epochs = 0)
qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'ga39'

qumap$pred = knn.reg(train = refmap[,1:2], 
    test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = rep(0, 4852)
qumap$diff = abs(qumap$pred - qumap$age)

qumap$sample = rep('ga39', 4852)
qumap$result = rep('correct', 4852)          
qumap$truecluster = qumap$cluster

g3 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(plot=g3, 'main_figures/ga39_umap.png')


cols = intersect(colnames(allmap), colnames(qumap))
allmap = rbind(allmap[,cols], qumap[,cols])

#add RL2986 - 261 day sample

latent = readRDS('snRNAseq/processed_data/transfer_learning/outs/RL2986_latent.RDS')
qlatent = readRDS('snRNAseq/processed_data/transfer_learning/outs/RL2986_qlatent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
                      init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
                      min_dist = 0.5, verbose = TRUE, n_epochs = 0)
qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'RL2986'

qumap$pred = knn.reg(train = refmap[,1:2], 
                     test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
                    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = rep(0.715, nrow(qumap))
qumap$diff = abs(qumap$pred - qumap$age)

qumap$sample = rep('261d', nrow(qumap))
qumap$result = rep('correct', nrow(qumap))          
qumap$truecluster = qumap$cluster

g4 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(plot=g4, 'main_figures/RL2986_umap.png')

cols = intersect(colnames(allmap), colnames(qumap))
allmap = rbind(allmap[,cols], qumap[,cols])

#add RL2987 - 443 day sample

latent = readRDS('snRNAseq/processed_data/transfer_learning/outs/RL2987_latent.RDS')
qlatent = readRDS('snRNAseq/processed_data/transfer_learning/outs/RL2987_qlatent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
                      init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
                      min_dist = 0.5, verbose = TRUE, n_epochs = 0)
qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'RL2986'

qumap$pred = knn.reg(train = refmap[,1:2], 
                     test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
                    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = rep(0.715, nrow(qumap))
qumap$diff = abs(qumap$pred - qumap$age)

qumap$sample = rep('261d', nrow(qumap))
qumap$result = rep('correct', nrow(qumap))          
qumap$truecluster = qumap$cluster

g5 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(plot=g5, 'main_figures/RL2987_umap.png')

cols = intersect(colnames(allmap), colnames(qumap))
allmap = rbind(allmap[,cols], qumap[,cols])

#plot age prediction density plots

astro = subset(qumap, cluster == 'Astro' & result == 'correct')
l23 = subset(qumap, cluster == 'L2/3_CUX2' & result == 'correct')

plotastro = ggplot(astro, aes(x = pred, y = sample, fill = ..x..))+ geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted astrocyte ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2,3, 4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
                     label = c('ga22', 
                               '34d', 
                               '301d', 
                               '3yr', '10yr',
                               '20yr', '40yr'))+
  scale_fill_viridis_c()



plotl23 = ggplot(l23, aes(x = pred, y = sample, fill = ..x..)) + geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted L2/3 ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2,3, 4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
                     label = c('ga22', 
                               '34d', 
                               '301d', 
                               '3yr', '10yr',
                               '20yr', '40yr'))+
  scale_fill_viridis_c()


ggsave(plotastro, 'main_figures/astro_age_prediction.png')
ggsave(plotl23, 'main_figures/l23_age_prediction.png')

### 1b. Organoids 

reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS')
initembedding = reducedDim(reference, "UMAP")

latent = readRDS('snRNAseq/processed_data/transfer_learning/outs/4999_latent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
    init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
    min_dist = 0.5, verbose = TRUE, n_epochs = 0)
refmap = as.data.frame(reference_umap$embedding)
refmap$dataset = rep('reference', 154748)

sce = list('Luc9228' = 'Luc9228', 
           'RL2290' = 'RL2290', 
           'RL2432' = 'RL2432')

allqumaps = list()

for(i in names(sce)){
    
  qlatent = readRDS(paste0('snRNAseq/processed_data/transfer_learning/outs/', i, '_qlatent.RDS'))
  qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)
  qumap = as.data.frame(qumap)
  qumap$dataset= rep(paste0(i), nrow(qumap))
  qumap$sample= rep(paste0(i), nrow(qumap))
  qumap$pred = knn.reg(train = refmap[,1:2], test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
  qumap$cluster = knn(train = refmap[,1:2], test = qumap[,1:2], cl = reference$cluster, k = 10)
  qumap$stage = knn(train = refmap[,1:2], test = qumap[,1:2], cl = reference$stage_ids, k =10)
  g6 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
    geom_point(col='black', size=1.5) +
    geom_point(col='white', size=0.8) +
    theme_classic() +
    geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
    geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                    adjust=1/2)
  
  ggsave(plot=g6, paste0('main_figures/', i, '_umap.png'))
  saveRDS(qumap, paste0('snRNAseq/processed_data/', 'i', '_umap.RDS'))
  allqumaps = rbind(allqumaps, qumap)
}

l23 = subset(allqumaps, allqumaps$cluster == 'L2/3_CUX2')
l23$dataset = factor(l23$dataset, levels = c('RL2432', 'Luc9228', 'RL2290'))

plotl23 = ggplot(l23, aes(x = pred, y = dataset, fill = ..x..)) + geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted L2/3 ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2,3, 4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
                     label = c('ga22', 
                               '34d', 
                               '301d', 
                               '3yr', '10yr',
                               '20yr', '40yr'))+
  scale_fill_viridis_c()


ggsave(plotl23, 'supp_figures/organoid_l23_age_prediction.png')


### 1c. GBM
latent = readRDS('snRNAseq/processed_data/transfer_learning/outs/GBM_latent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
    init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
    min_dist = 0.5, verbose = TRUE, n_epochs = 0)
refmap = as.data.frame(reference_umap$embedding)
refmap$dataset = rep('reference', 154748)
    
  qlatent = readRDS(paste0('snRNAseq/processed_data/transfer_learning/outs/GBM_qlatent.RDS'))
  qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)
  qumap = as.data.frame(qumap)
  qumap$dataset= rep('GBM', nrow(qumap))
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
  
  ggsave(plot=g1, paste0('main_figures/GBM_umap.png'))
  
saveRDS(qumap, 'snRNAseq/processed_data/GBM_umap.RDS')


sce <- readRDS("snRNAseq/processed_data/gbm_primary_dataset.RDS")
sce <- logNormCounts(sce)

stem_cell_core <- read.table("annotation/stem_cell_core.txt",
                             sep="\t", skip=1, fill=T)
stem_cell_core <- stem_cell_core[1:140,]

g2m.genes <- cc.genes$g2m.genes
ind.g2m.genes <- match(g2m.genes, rowData(sce)$Symbol)
ind.g2m.genes <- ind.g2m.genes[!is.na(ind.g2m.genes)]

ind <- match(stem_cell_core[,3], rowData(sce)$Symbol)
ind <- ind[!is.na(ind)]
eigen <- prcomp(as.matrix(logcounts(sce[ind,])))

df <- reducedDim(sce, "UMAP")
colnames(df) <- c("UMAP1", "UMAP2")
df <- as.data.frame(df)
df$eigen <- eigen$rotation[,1]
df$mean <- colMeans(as.matrix(logcounts(sce[ind,])))
df$g2m_genes <- colMeans(as.matrix(logcounts(sce[ind.g2m.genes,])))
set.seed(2)
df <- df[sample(1:dim(df)[1]),]


logo_file <- "snRNAseq/processed_data/background.png"
img <- png::readPNG(logo_file)

p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = g2m_genes, fill = stat(mean)),
     funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title="G2M Phase"))

p2 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = eigen, fill = stat(mean)),
    funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title="Stem Cell Core Network"))
        

ggsave(p1, file="supp_figures/gbm_g2m.png",
       height=7, width=7)
ggsave(p2, file="supp_figures/gbm_stem_cell_score.png",
       height=7, width=7)


# 2. Independent organoid analysis

organoids <- c("snRNAseq/processed_data/orgpred/Luc9228.RDS",
               "snRNAseq/processed_data/orgpred/org2290.RDS",
               "snRNAseq/processed_data/orgpred/org2432.RDS")
samples <- c(`9_months`="LucH9228", `12_months`="RL2290", `5_months`="RL2432")

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

org_sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

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
        geom_point(size=1, alpha=0.75) + theme_classic() +
        scale_color_manual(values=col) + 
        guides(color= FALSE)
    
    g2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=stage)) + 
        geom_point(size=1, alpha=0.75) + theme_classic() +
        scale_color_manual(values=col_stage) +
        guides(color= FALSE)
    
    p1 <- plot_grid(g1, g2)
    ggsave(paste0("supp_figures/organoids_single_",
        names(samples)[i], ".png"), p1,  height=5, width=10, device="png")
    
}
    
    
# Combined analysis   
set.seed(20)
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
                   BSPARAM=BiocSingular::RandomParam(), ncomponents=50)
percent.var <- attr(reducedDim(corrected), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")

nn.clusters <- clusterCells(corrected, use.dimred="PCA",  
      BLUSPARAM=NNGraphParam(k=20))
table(nn.clusters)
corrected$cluster <- nn.clusters

corrected <- runUMAP(corrected, dimred="PCA", spread=0.2, min_dist=0.075)
corrected$predcelltype <- uncorrected$predcelltype
corrected$predstage <- uncorrected$predstage

df <- reducedDim(corrected, "UMAP")
colnames(df) <- c("UMAP1", "UMAP2")
df <- as.data.frame(df)
df$cell_type <- corrected$predcelltype
df$stage <- corrected$predstage
df$batch <- corrected$batch
df$cluster <- corrected$cluster

make_marker_plot <- function(gene, df){
  
  ind <- which(rowData(rescaled)$Symbol==gene)
  if(length(ind)==0){
    stop("Gene not found!")
  }
  df$gene <- as.matrix(t(logcounts(rescaled[ind,])))
  set.seed(30)
  df <- df[sample(1:dim(df)[1]),]
  ggplot(df, aes(x=UMAP1, y=UMAP2, col=gene)) + 
    geom_point(size=0.5, alpha=0.5) + theme_classic() +
    scale_color_viridis_c(name = gene) +
    theme(legend.title = element_text(face = "italic"))
  
}

g_gad1 <- make_marker_plot("GAD1", df)
g_cux2 <- make_marker_plot("CUX2", df)
g_slc1a3 <- make_marker_plot("SLC1A3", df)
g_pdgfrb <- make_marker_plot("PDGFRB", df)
g_meis2 <- make_marker_plot("MEIS2", df)
g_unc5d <- make_marker_plot("UNC5D", df)

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
    geom_point(size=1, alpha=0.75) + theme_classic() +
    scale_color_manual(values=col) + 
    guides(color= FALSE) 

g1a <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=cluster)) + 
    geom_point(size=1, alpha=0.75) + theme_classic() + 
     guides(color= FALSE)


leg <- get_legend( ggplot(df, aes(x=UMAP1, y=UMAP2, col=cell_type)) + 
    geom_point(size=1, alpha=0.75) + theme_classic() +
    scale_color_manual(values=col) + theme(legend.position="bottom",
        legend.title = element_blank())+
      guides(color=guide_legend(nrow=10,ncol=2,byrow=TRUE)))

lega <- get_legend( ggplot(df, aes(x=UMAP1, y=UMAP2, col=cluster)) + 
    geom_point(size=1, alpha=0.75) + theme_classic() + 
      theme(legend.position="bottom",  legend.title = element_blank())) 


p1 <- plot_grid(lega, NULL, ncol=1)

plot_grid(g1, g1a, leg, p1, labels=c("A", "B", NULL, NULL))
plot_grid(g_gad1, g_cux2,  g_slc1a3, g_pdgfrb, g_meis2,
g_unc5d,  ncol=2, labels=c("C", "D", "E", "F", "G", "H"))



g2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=stage)) + 
    geom_point(si ze=1, alpha=0.75) + theme_classic() +
    scale_color_manual(values=col_stage) + 
    guides(color= FALSE)

g3 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=batch)) + 
    geom_point(size=1, alpha=0.75) + theme_classic()  + 
    guides(color= FALSE)

p1 <- plot_grid(g1,g1a,g2,g3, ncol=2, nrow=2)
ggsave(p1, file="supp_figures/organoid_combined.png", height=7, width=7)


saveRDS(list(rescaled=rescaled, corrected=corrected), 
        file="snRNAseq/processed_data/organoid_independent.rds")

# 3. Differentially expressed genes

### 3a. number overlap

#### mature PN vs immature PN in organoids

markers_PN <- read.table(file="snRNAseq/processed_data/organoids_PN_markers.csv",
             sep=",", header = T)
gene_trends <-  read.csv("snRNAseq/processed_data/gene_cluster_ids.csv")
gene_trends <- gene_trends[gene_trends$major_clust %in% 
  c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4"),]

list_gene <- split(gene_trends$gene_name, gene_trends$major_trend)
list_organoids <- list()
list_organoids[["mature"]] <- markers_PN$Symbol[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC>0]
list_organoids[["immature"]] <- markers_PN$Symbol[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC<0]

tmp_overlap <- sapply(list_organoids, function(x) sapply(list_gene, function(y)
  length(intersect(x,y))))
tmp_overlap <- apply(tmp_overlap, 1, function(x) x/lengths(list_organoids))

gg_tmp <- reshape2::melt(tmp_overlap)
colnames(gg_tmp) <- c("Organoid PN Gene", "Major Gene Trend", "Overlap")
gg_tmp$`Major Gene Trend` <- as.character(gg_tmp$`Major Gene Trend`)
gg_tmp$`Major Gene Trend`[gg_tmp$`Major Gene Trend`=="interup"] <- "trans up"
gg_tmp$`Major Gene Trend`[gg_tmp$`Major Gene Trend`=="interdown"] <- "trans down"

g1 <- ggplot(gg_tmp, aes(x=`Organoid PN Gene`, y=`Major Gene Trend`, fill=Overlap)) +
  geom_tile()  + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours =colorspace::diverge_hcl(7) , 
  na.value = "transparent") 
  
ggsave(g1, file="supp_figures/compare_organoid_PN_markers_DEGs.svg",
       width=3, height=3)
  
list_gene <- split(gene_trends$gene_name, gene_trends$gene_trend)
list_organoids <- list()
list_organoids[["mature"]] <- markers_PN$Symbol[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC>0]
list_organoids[["immature"]] <- markers_PN$Symbol[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC<0]

genes_g4 <- intersect(list_organoids$mature, list_gene$G4)
write.csv(genes_g4, "snRNAseq/processed_data/Genes_mature_PNs_organoids.csv",
    quote=FALSE, row.names=FALSE)

tmp_overlap <- sapply(list_organoids, function(x) sapply(list_gene, function(y)
  length(intersect(x,y))))
tmp_overlap <- apply(tmp_overlap, 1, function(x) x/lengths(list_organoids))

gg_tmp <- reshape2::melt(tmp_overlap)
colnames(gg_tmp) <- c("Organoid PN Gene", "Gene Trend", "Overlap")
gg_tmp$`Gene Trend` <- factor(gg_tmp$`Gene Trend`, levels=
      rev(paste0("G", 1:14)))


g1a <- ggplot(gg_tmp, aes(x=`Organoid PN Gene`, y=`Gene Trend`, fill=Overlap)) +
  geom_tile()  + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours =colorspace::diverge_hcl(7) , 
  na.value = "transparent") 
  
ggsave(g1a, file="supp_figures/compare_organoid_PN_markers_DEGs_fine.svg",
       width=3, height=3)


plot_grid(g1a, g1, labels = c('A', 'B'),  hjust = -0.1,
          label_size = 12)


#### mature PN in brain vs mature PN in organoids

markers_PN <- read.table(file="snRNAseq/processed_data/Organoids_Brain_markers_PN.csv",
             sep=",", header = T)
gene_trends <-  read.csv("snRNAseq/processed_data/gene_cluster_ids.csv")
gene_trends <- gene_trends[gene_trends$major_clust %in% 
  c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4"),]

list_gene <- split(gene_trends$gene_name, gene_trends$major_trend)
list_organoids <- list()
list_organoids[["organoid up"]] <- rownames(markers_PN)[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC>0]
list_organoids[["brain up"]] <- rownames(markers_PN)[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC<0]

tmp_overlap <- sapply(list_organoids, function(x) sapply(list_gene, function(y)
  length(intersect(x,y))))
tmp_overlap <- apply(tmp_overlap, 1, function(x) x/lengths(list_organoids))

gg_tmp <- reshape2::melt(tmp_overlap)
colnames(gg_tmp) <- c("Organoid PN Gene", "Major Gene Trend", "Overlap")

g1 <- ggplot(gg_tmp, aes(x=`Organoid PN Gene`, y=`Major Gene Trend`, fill=Overlap)) +
  geom_tile()  + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours =colorspace::diverge_hcl(7) , 
  na.value = "transparent") 
  
ggsave(g1, file="supp_figures/compare_organoid_brain_markers_DEGs.svg",
       width=5, height=7)
  
list_gene <- split(gene_trends$gene_name, gene_trends$gene_trend)
list_organoids <- list()
list_organoids[["organoids up"]] <- rownames(markers_PN)[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC>0]
list_organoids[["brain up"]] <- rownames(markers_PN)[markers_PN$FDR<0.05 &
    markers_PN$summary.logFC<0]

tmp_overlap <- sapply(list_organoids, function(x) sapply(list_gene, function(y)
  length(intersect(x,y))))
tmp_overlap <- apply(tmp_overlap, 1, function(x) x/lengths(list_organoids))

gg_tmp <- reshape2::melt(tmp_overlap)
colnames(gg_tmp) <- c("Organoid PN Gene", "Gene Trend", "Overlap")
gg_tmp$`Gene Trend` <- factor(gg_tmp$`Gene Trend`, levels=
      paste0("G", 1:14))

g1 <- ggplot(gg_tmp, aes(x=`Organoid PN Gene`, y=`Gene Trend`, fill=Overlap)) +
  geom_tile()  + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours =colorspace::diverge_hcl(7) , 
  na.value = "transparent") 
  
ggsave(g1, file="supp_figures/compare_organoid_brain_markers_DEGs_fine.svg",
       width=5, height=7)

### 3B. Volcano plots

#### Mature versus immature organoid PNs

markers_PN <- read.table(file="snRNAseq/processed_data/Organoids_PN_markers.csv",
             sep=",", header = T)

genes_of_interest <- c("KCNB2", "KCNMA1", "GABRB1", "GRIA4", "LRRTM4", "PCDH9")

markers_PN$p.adjust <- -log10(markers_PN$FDR)
g1 <- ggplot(markers_PN, aes(x=summary.logFC, y=p.adjust)) + 
    geom_point(col="black", alpha=0.5) +
    geom_point(data=markers_PN[markers_PN$summary.logFC > 0 &
      markers_PN$FDR < 0.05,], col="#990099", alpha=0.5) +
    geom_point(data=markers_PN[markers_PN$summary.logFC < 0 &
        markers_PN$FDR < 0.05,], col="#009900", alpha=0.5) +
    theme_classic() +
    geom_text_repel(data=markers_PN[markers_PN$Symbol %in% genes_of_interest,],
        col="grey", aes(label=Symbol),
        size=4, max.overlaps = 20) +
    geom_point(data=markers_PN[markers_PN$Symbol %in% genes_of_interest,], 
               col="red", alpha=0.5) 

ggsave(g1, file="main_figures/volcano_pN.png", width=7,
       height=7)

##### organoid versus brain mature PNs   
  
markers_organoid <- read.table(file="snRNAseq/processed_data/Organoids_Brain_markers_PN.csv",
             sep=",", header = T)
markers_organoid$summary.logFC <- -markers_organoid$summary.logFC
markers_organoid$Symbol <- rownames(markers_organoid)
genes_of_interest <- c("MEG3", "GRIN2B", "GRIA1", "KCND2", "KCNQ5")

markers_organoid$p.adjust <- -log10(markers_organoid$FDR)
g1 <- ggplot(markers_organoid, aes(x=summary.logFC, y=p.adjust)) + 
    geom_point(col="black", alpha=0.5) +
    geom_point(data=markers_organoid[markers_organoid$summary.logFC > 0 &
      markers_organoid$FDR < 0.05,], col="cornflowerblue", alpha=0.5) +
    geom_point(data=markers_organoid[markers_organoid$summary.logFC < 0 &
        markers_organoid$FDR < 0.05,], col="#990099", alpha=0.5) +
    theme_classic() +
    #geom_text_repel(data=markers_organoid[markers_organoid$Symbol %in% genes_of_interest,],
    #    col="grey", aes(label=Symbol),
    #    size=4, max.overlaps = 20) +
    geom_point(data=markers_organoid[markers_organoid$Symbol %in% genes_of_interest,], 
               col="red", alpha=0.5) 

ggsave(g1, file="main_figures/volcano_organoid_brain.png", width=7,
       height=7)
  

### 3c. GO plots

#### mature versus immature PN organoids

go_pn_up <- read.table("snRNAseq/processed_data/GO_PN_up.tsv", sep="\t")
go_pn_down <- read.table("snRNAseq/processed_data/GO_PN_down.tsv", sep="\t")
go_pn <- read.table("snRNAseq/processed_data/GO_PN.tsv", sep="\t")

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

up_paths <- c("GO:0042391", "GO:0050808", "GO:0007215", "GO:0006836",
              "GO:0001508", "GO:0048167", "GO:0007214", "GO:2001257",
              "GO:0007611", "GO:0001764", "GO:0036465", "GO:1990138")
down_paths <- c("GO:0000819", "GO:0006302", "GO:0006260", "GO:0048285",
                "GO:0061351", "GO:0098727", "GO:0019827", "GO:0050768",
                "GO:0051961", "GO:0006338", "GO:0007405", "GO:0042063")

go_pn_up1 <- go_pn_up[up_paths,]
go_pn_up1$PValue <- 1

g1 <- ggplot(go_pn_up1, aes(y=Description, x=PValue, colour=p.adjust, size=Count)) + 
  geom_point() + theme_classic() +  scale_colour_gradientn(limits=c(0, 0.05), 
    colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
  scale_size(limits=c(0, 80))
ggsave(g1, file="main_figures/organoids_GO_PN_up.svg", width=7, height=7)

go_pn_down1 <- go_pn_down[down_paths,]
go_pn_down1$PValue <- 1
g1 <- ggplot(go_pn_down1, aes(y=Description, x=PValue, colour=p.adjust, size=Count)) + 
  geom_point() + theme_classic() +  scale_colour_gradientn(limits=c(0, 0.05), 
    colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    scale_size(limits=c(0, 80))
ggsave(g1, file="main_figures/organoids_GO_PN_down.svg", width=7, height=7)

#### mature organoid versus mature brain PNs

go_pn_up <- read.table("snRNAseq/processed_data/GO_brain_organoid_down.tsv", sep="\t")
go_pn_up <- go_pn_up[, -8]
go_pn_down <- read.table("snRNAseq/processed_data/GO_brain_organoid_up.tsv", sep="\t")

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

up_paths <- go_pn_up$ID[1:10]
down_paths <-  go_pn_down$ID[1:10]

go_pn_up1 <- go_pn_up[up_paths,]
go_pn_up1$PValue <- 1
g1 <- ggplot(go_pn_up1, aes(y=Description, x=PValue, colour=p.adjust, size=Count)) + 
  geom_point() + theme_classic() +  scale_colour_gradientn(limits=c(0, 0.05), 
    colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    scale_size(limits=c(0, 300))
ggsave(g1, file="main_figures/brain_vs_organoids_GO_PN_up.svg", width=7, height=7)

go_pn_down1 <- go_pn_down[down_paths,]
g1 <- ggplot(go_pn_down1, aes(y=Description, x=Count, fill=p.adjust)) +
  geom_col() +
  scale_fill_gradientn(limits=c(0, 0.05), 
   colours = rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
  theme_minimal() + theme(axis.title.y = element_blank())
ggsave(g1, file="supp_figures/brain_vs_organoids_GO_PN_down.svg", width=7, height=7)


### 3c UMAP plots with PNs

#### mature versus immature PNs in organoids

path1 <- "snRNAseq/processed_data/orgpred/"
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
ggsave(g1, file="main_figures/ogranoids_umap_PN.png", width=7, height=7)     


##### brain versus organoid mature PNs 

umap_org <- readRDS(umap_org, file="snRNAseq/processed_data/umap_orginial.RDS")
umap_organoid <- readRDS(file="snRNAseq/processed_data/umap_organoid.RDS")

g1 <- ggplot(umap_org[,1:2], aes(x=UMAP1, y=UMAP2)) +
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) +
    theme_classic() +
    geom_point(data=umap_org[ind,1:2], 
        alpha=0.5, size=1, col="cornflowerblue") +
    geom_point(data=umap_organoid[,1:2], alpha=0.5, size=0.5, col="#990099")


ggsave(g1, file="main_figures/ogranoids_brain_umap_PN.png", width=7, height=7)
