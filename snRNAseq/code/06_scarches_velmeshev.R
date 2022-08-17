#SCVI workthrough

#0:library loading
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

source("snATACseq/R/knn_function.R")

sc = import('scanpy', convert = FALSE)
scvi = import('scvi', convert = FALSE)

#1: reference downsampling

refmake = readRDS('snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS')
refmake = refmake[rowSums(counts(refmake)>0)>=5,]
bycol = 1000/colSums(counts(refmake))            #ensures downsample target is 1000
set.seed(100)
downsampled = downsampleMatrix(x = counts(refmake), prop = bycol, bycol = TRUE)
assay(refmake, 'downsample1000') = downsampled

oldref = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS') 
rowSubset(refmake, 'highly_variable') = subset(rowData(oldref)$index, 
        rowData(oldref)$highly_variable == TRUE)
refmake@colData = oldref@colData[c(76,71,69,68,66,62,59,57)]
refmake = refmake[subset(rowData(oldref)$index, 
        rowData(oldref)$highly_variable == TRUE),]

#2:  make query
ribogeneset = read.csv('annotation/ribo_geneset.txt', header = FALSE)$V1 

query = readRDS('snRNAseq/processed_data/transfer_learning/datasets/velmeshev_dataset.RDS') 
# this dataset can be downloaded from https://www.science.org/doi/abs/10.1126/science.aav8130


ribout = subset(rowData(query)$ID, rowData(query)$Symbol %in% ribogeneset)

query = query[!( rownames(query) %in% ribout),]
query = query[rowData(query)$Symbol != 'MALAT1',]
query =query[-which(str_detect(rowData(query)$Symbol, '^MT-')),]
query = query[rowSums(counts(query)>0)>=5,]

qbycol = 1000/colSums(counts(query))
set.seed(100)
qdownsampled = downsampleMatrix(x = counts(query), prop = qbycol, bycol = TRUE)
assay(query, 'downsample1000') = qdownsampled

smolquery = query[intersect(rownames(query), rowData(refmake)$gene_ids),]  

#3: create reference scvi - takes a while 

smolrefmake = refmake[subset(rowData(refmake)$index, 
    rowData(refmake)$gene_ids %in% intersect(rownames(query), 
    rowData(refmake)$gene_ids)),]
rownames(smolrefmake) = rowData(smolrefmake)$gene_ids

rdata = sc$AnnData(X = t(assay(smolrefmake, 'downsample1000')), 
    obs = as.data.frame(smolrefmake@colData))
scvi$data$setup_anndata(rdata)  
vae_ref = scvi$model$SCVI(rdata, use_layer_norm = 'both', 
    use_batch_norm = 'none', encode_covariates = TRUE, 
    dropout_rate = as.numeric(0.2), n_layers = as.integer(2), 
    n_latent = as.integer(365))
refstart = proc.time()
vae_ref$train(max_epochs = as.integer(150), early_stopping = TRUE)
refend = proc.time()
reftime = refend - refstart 
latent = vae_ref$get_latent_representation()
latent = as.matrix(latent)
rownames(latent) = colnames(smolrefmake)


#4 create query scvi

qdata = sc$AnnData(X = t(assay(smolquery, 'downsample1000')), 
                 obs = as.data.frame(smolquery@colData))
scvi$data$setup_anndata(qdata)
vae_q = scvi$model$SCVI$load_query_data(qdata, vae_ref)
querystart = proc.time()
vae_q$train(max_epochs = as.integer(150), 
    plan_kwargs=dict(weight_decay = as.integer(0)), 
    early_stopping = TRUE)
queryend = proc.time()
querytime = queryend - querystart
RunTime = querytime + reftime
qlatent = vae_q$get_latent_representation()
qlatent = as.matrix(qlatent)
rownames(qlatent) = colnames(smolquery)


#5 save latent representations for further analyses

saveRDS(latent, 'snRNAseq/processed_data/transfer_learning/outs/velmeshev_latent.RDS')
saveRDS(qlatent, 'snRNAseq/processed_data/transfer_learning/outs/velmeshev_qlatent.RDS')

#6: umap transform

reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS') 
initembedding = reducedDim(reference, "UMAP")

reference_umap = umap(latent, n_neighbors = 25, 
    ret_model = TRUE, init = initembedding, a = 0.58303, b = 1.334167, 
    spread = 1, min_dist = 0.5, verbose = TRUE, n_epochs = 0)
qumap = umap_transform(qlatent, reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'query'
commap = rbind(refmap, qumap)

qumap$pred = knn.reg(train = refmap[,1:2], test = qumap[,1:2], 
    y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], test = qumap[,1:2], 
    cl = reference$cluster, k = 10)
qumap$age = query$arcsin_ages
qumap$diff = abs(qumap$pred - qumap$age)

#add samples - different for Velmeshev
qumap$sample = query$sample
qumap$sample = fct_reorder(qumap$sample, qumap$age)
qumap = subset(qumap, !(sample %in% c('5893_PFC', '5538_PFC_Nova', 
    '5879_PFC_Nova','5936_PFC_Nova', '5408_PFC_Nova', '5144_PFC', 
    '5278_PFC', '5403_PFC', '5419_PFC', '5945_PFC')))
smolquery = query[, !(query$sample %in% c('5893_PFC', '5538_PFC_Nova', 
    '5879_PFC_Nova','5936_PFC_Nova', '5408_PFC_Nova', '5144_PFC', 
    '5278_PFC', '5403_PFC', '5419_PFC', '5945_PFC'))]
out = myknn(ref_embed = refmap[,1:2], test_embed = qumap[,1:2], 
    test_class = smolquery$cluster, ref_class = reference$cluster, k =10)
qumap$result = out$result
qumap$truecluster = smolquery$cluster


saveRDS(qumap, file="snRNAseq/processed_data/velmeshev_umap.RDS")

