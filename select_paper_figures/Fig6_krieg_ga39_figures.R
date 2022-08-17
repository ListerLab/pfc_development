#plotting figures

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

#velmeshev plotting

query = readRDS('snRNAseq/processed_data/kriegstein_dataset.RDS')
# this data can be downloaded from https://www.science.org/doi/abs/10.1126/science.aav8130

qumap <- readRDS("snRNAseq/processed_data/scArches_umap_kriegstein.RDS")
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

ggsave(plot=g1, 'paper_figures/krieg_umap.png')

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

ggsave(plot=g2, 'paper_figures/krieg_clusters.png')


#add ga39 sample

ga39 = readRDS('snRNAseq//ga39_dataset.RDS')
latent = readRDS('snRNA/processed_data/4999_latent.RDS')
qlatent = readRDS('snRNA/processed_data/ga39_qlatent.RDS')

reference_umap = umap(latent, n_neighbors = 25, ret_model = TRUE, 
    init = initembedding, a = 0.58303, b = 1.334167, spread = 1, 
    min_dist = 0.5, verbose = TRUE, n_epochs = 0)
refmap = as.data.frame(reference_umap$embedding)
qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
refmap$dataset = character(154748)
refmap$dataset[1:154748] = 'reference'
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'ga39'

qumap$pred = knn.reg(train = refmap[,1:2], 
    test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = rep(0, 4852)
qumap$diff = abs(qumap$pred - qumap$age)

#add samples - and remove non-specific samples
qumap$sample = rep('ga39', 4852)
qumap$result = rep('correct', 4852)          
qumap$truecluster = qumap$cluster

g1 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(plot=g1, file='paper_figures/ga39_umap.png') 

qumap = rbind(qumap, kriegmap)
qumap$sample = factor(qumap$sample, levels = names(qumap$sample))

astro = subset(qumap, cluster == 'Astro' & result == 'correct')
l23 = subset(qumap, cluster == 'L2/3_CUX2' & result == 'correct')

plotastro = ggplot(astro, aes(x = pred, y = sample, fill = ..x..)) + 
    geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted astrocyte ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = 
        c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2, 3, 
        4,6,7,8,10,11, 12,14,15,16,18, 19, 20,21,23)]),
        label = c('ga22','34d', '301d', '3yr', '10yr','20yr', '40yr')) +
  scale_fill_viridis_c()

ggsave(plot=g2, 'paper_figures/krieg_ga39_astro.png')

plotl23 = ggplot(l23, aes(x = pred, y = sample, fill = ..x..)) + 
    geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted L2/3 ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(
        levels(as.factor(reference$arcsin_ages)))[-c(2,3, 
        4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
        label = c('ga22','34d', '301d', '3yr', '10yr', '20yr', '40yr')) +
  scale_fill_viridis_c()

ggsave(plot=g2, 'paper_figures/krieg_ga39_l23.png')








