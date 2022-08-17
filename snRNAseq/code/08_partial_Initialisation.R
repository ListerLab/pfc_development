#Partial Initialisation: Checking the suitability of datasets for integration

#Step 1: load packages
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
library(dplyr)

sc = import('scanpy', convert = FALSE)
scvi = import('scvi', convert = FALSE)
numpy = import('numpy', convert = FALSE)
mp = import('matplotlib', convert = FALSE)
mp$use('Agg')


#Step 2: load and combine data: select relevant latent representations for reference and query
reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS')
qlatent <- readRDS('snRNAseq/processed_data/transfer_learning/outs/RL2290_qlatent.RDS')
reflatent <- readRDS('snRNA/processed_data/transfer_learning/outs/4999_latent.RDS')
latent = rbind(reflatent, qlatent)

#Step 3: Initialise data - reference nuclei initialised to positions in UMAP, query nuclei randomly distributed over range of reference

set.seed(100)
initembedding <- reducedDim(reference, "UMAP")
initembedding = rbind(initembedding, matrix(0, nrow = nrow(qlatent), ncol = 2))
initembedding[(nrow(reflatent)+1):nrow(initembedding),1] <-  runif(nrow(qlatent), -6, 21.5)
initembedding[(nrow(reflatent)+1):nrow(initembedding),2] <-  runif(nrow(qlatent), -9.9, 15)

#Step 4: visualise partial initialisation with UMAP
#Thinking to remove this section since we don't show the UMAPs in figures/supps

reference_umap <-  umap(latent, n_neighbors  <-  25, ret_model = TRUE, init = initembedding, a = 0.58303, b = 1.334167, spread = 1, min_dist = 0.5, verbose = TRUE, n_epochs = 50)
fullmap <- as.data.frame(reference_umap$embedding)
refmap <-  fullmap[1:nrow(reflatent),]
qumap <- fullmap[(nrow(reflatent)+1):nrow(fullmap),]
refmap$dataset[1:nrow(refmap)]<- 'reference'
qumap$dataset <- character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] <- 'query'
refmap$cluster <- reference$cluster

g1 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2)) +
  geom_point(col='black', size=1.5) +
  geom_point(col='white', size=0.8) +
  theme_classic() +
  geom_point(data=qumap, alpha=0.5, size=1, col='#990000')  +
  geom_density_2d(data=qumap, color='#FF7B00', size=0.5, bins=10,
                  adjust=1/2)

ggsave(g1, file = 'snRNAseq/figures/RL2290_partial_init.png')


g2 <- ggplot(refmap[sample(1:nrow(refmap)), ], aes(x=V1, y=V2, col = cluster)) +
  geom_point( size=1.5) +
  theme_classic() +
  scale_colour_manual(values = c('Astro' ='#FFC857','CGE_dev'='#C6D5C0',
                                 'ID2'='#558140','L2/3_CUX2'='#6E0614',
                                 'L4_RORB' ='#8B3843', 'L5/6_THEMIS_TLE4' = '#A86A72',
                                 'MGE_dev' ='#ECD1C8','Micro' ='#484848',
                                 'Oligo'= '#255F85','OPC' ='#92AFC2',
                                 'PN_dev' ='#E2CDD0', 'Poor-Quality' ='#F2F2F2',
                                 'PV'= '#C77459', 'SST' ='#B44622',
                                 'VIP'= '#1C5701', 'Vas' ='#A3A3A3'))


ggsave(g2, file = 'snRNAseq/figures/RL2290_partial_init_ref.png')


#Step 5: Prepare Anndata object
m <-  Matrix(nrow = nrow(latent), sparse = TRUE)
rdata <-  sc$AnnData(X = m, obsm = list(X_scVI = latent))

#Step 6: Compute nearest neighbours, perform leiden clustering and calculate PAGA connectivities
sc$pp$neighbors(rdata, use_rep = 'X_scVI')
sc$tl$leiden(rdata)
sc$tl$paga(rdata)

connectivities <- py_to_r(rdata$uns[['paga']]['connectivities'])
groups <- py_to_r(numpy$array(rdata$obs[['leiden']]))
groups <-  factor(groups, levels = c(0:max(as.numeric(groups))))

rownames(connectivities)  <- c(0:(nrow(connectivities)-1))
colnames(connectivities) <-  c(0:(nrow(connectivities)-1))

refmap$groups <- groups[1:nrow(reflatent)]
qumap$groups <- groups[(nrow(reflatent)+1):length(groups)]
reference$paga_groups <-  groups[1:154748]
fullmap$groups <- groups


#Step 7: position leiden cluster centroids, and create table of PAGA connectivities (ignore warning message!)

group_centers <- fullmap %>%
  group_by(groups) %>%
  summarize(x = median(V1), y = median(V2)) %>%
  arrange(groups)


paga_position = group_centers
rownames(paga_position) <- c(0:(nrow(group_centers)-1))



paga_edges <- tibble(
  group1 = rownames(connectivities)[row(connectivities)[upper.tri(connectivities)]],
  group2 = colnames(connectivities)[col(connectivities)[upper.tri(connectivities)]],
  weight = connectivities[upper.tri(connectivities)]
) %>%
  mutate(
    x1 = paga_position$x[match(.$group1, rownames(paga_position))],
    y1 = paga_position$y[match(.$group1, rownames(paga_position))],
    x2 = paga_position$x[match(.$group2, rownames(paga_position))],
    y2 = paga_position$y[match(.$group2, rownames(paga_position))]
  ) %>%
  filter(weight > 0)


saveRDS(paga_edges, 'snRNAseq/processed_data/paga/RL2290_paga_edges.RDS')


#Step 8: Plot connection weights by cell type
paga_query <-  paga_edges[(paga_edges$group1 %in% q_clusters & !(paga_edges$group2 %in% q_clusters))  |(paga_edges$group2 %in% q_clusters & !(paga_edges$group1 %in% q_clusters)),]

paga_query$celltype <-  character(nrow(paga_query))
paga_query$cluster <-  character(nrow(paga_query))

for(i in 1:nrow(paga_query)){
  paga_query$cluster[i] <- paga_query[[i,which(!(paga_query[i,1:2] %in% q_clusters))]]}

celltypes = list()

for(i in 1:nrow(paga_query)){
  if(paga_query$group1[i] %in% q_clusters){
    celltypes[[i]] = names(which(summary(refmap[refmap$groups %in% paga_query$group2[i],]$cluster)>100))}
  
  if(!(paga_query$group1[i] %in% q_clusters)){
    celltypes[[i]] = names(which(summary(refmap[refmap$groups %in% paga_query$group1[i],]$cluster)>100))}
}

weights = c()

for(i in 1:length(celltypes)){
  weights = c(weights, rep(paga_query$weight[i], length(celltypes[[i]])))
}

df = data.frame(weight = weights, celltype = unlist(celltypes), dataset = rep('12 month', length(weights)))

g3 = ggplot(df, aes(x=dataset, y = weight, col = celltype)) + geom_quasirandom( cex = 16, width = 0.1, bandwidth = 0.05, dodge.width = 0.5)+ theme_classic()+
  scale_color_manual(name = 'Cell Type', values = c('Astro' ='#FFC857','CGE_dev'='#C6D5C0',
                                                    'ID2'='#558140','L2/3_CUX2'='#6E0614',
                                                    'L4_RORB' ='#8B3843', 'L5/6_THEMIS_TLE4' = '#A86A72',
                                                    'MGE_dev' ='#ECD1C8','Micro' ='#484848',
                                                    'Oligo'= '#255F85','OPC' ='#92AFC2',
                                                    'PN_dev' ='#E2CDD0', 'Poor-Quality' ='#F2F2F2',
                                                    'PV'= '#C77459', 'SST' ='#B44622',
                                                    'VIP'= '#1C5701', 'Vas' ='#A3A3A3'))+ 
  geom_hline(yintercept= 0.15, color='red', linetype = 'dashed') + ylim(0,1)+
  xlab('Sample')+ ylab('Connection Weight')+
  theme(axis.text = element_text(size =14), axis.title=element_text(size = 20))


ggsave(g3, file = 'supp_figures/RL2290_paga_points.png')



##For additional samples, follow the same procedure with the appropriate latent representations
