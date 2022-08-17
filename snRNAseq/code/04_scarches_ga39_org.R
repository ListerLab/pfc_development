#organoid/ga39/GBM scArches


library(BiocManager)
library(scater)
library(scran)
library(scuttle)
library(reticulate)
library(Seurat)
library(uwot)
library(class)
library(DropletUtils)


sc = import('scanpy', convert = FALSE)
scvi = import('scvi', convert = FALSE)

#load data

sce = list(Luc9228 = readRDS('snRNAseq/processed_data/transfer_learning/datasets/LucH9228'), 
           RL2290 = readRDS('snRNAseq/processed_data/transfer_learning/datasets/RL2290'), 
           RL2432 =readRDS('snRNAseq/processed_data/transfer_learning/datasets/RL2432.RDS'), 
           ga39 = readRDS('snRNAseq/processed_data/transfer_learning/datasets/ga39_dataset.RDS'),
           GBM = readRDS('snRNAseq/processed_data/transfer_learning/datasets/gbm_primary_dataset.RDS'))

#reference downsampling

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


## make scvi reference - same for all samples in sce, as all contain the 4999 hvg genes
rdata = sc$AnnData(X = t(assay(refmake, 'downsample1000')), 
                   obs = as.data.frame(refmake@colData))
scvi$data$setup_anndata(rdata)
vae_ref = scvi$model$SCVI(rdata, use_layer_norm = 'both', 
                          use_batch_norm = 'none', encode_covariates = TRUE, 
                          dropout_rate = as.numeric(0.2), n_layers = as.integer(2), 
                          n_latent = as.integer(30))
refstart = proc.time()
vae_ref$train(max_epochs = as.integer(150), early_stopping = TRUE)
refend = proc.time()
reftime = refend - refstart
latent = vae_ref$get_latent_representation()
latent = as.matrix(latent)
rownames(latent) = colnames(smolrefmake)

saveRDS(latent, 'snRNAseq/processed_data/transfer_learning/outs/4999_latent.RDS')


## preprocessing for all query samples

for(i in names(sce)){

    bycol = 1000/colSums(counts(sce[[i]]))            #ensures downsample target is 1000
    set.seed(100)
    downsampled = downsampleMatrix(x = counts(sce[[i]]), 
        prop = bycol, bycol = TRUE)
    assay(sce[[i]], 'downsample1000') = downsampled

    sce[[i]] = subset(sce[[i]], rowData(sce[[i]])$ID %in% hvg)
    
}


## query scArches

for(i in names(sce)){
    
    qdata = sc$AnnData(X = t(assay(sce[[i]], 'downsample1000')), 
        obs = as.data.frame(sce[[i]]@colData))
    scvi$data$setup_anndata(qdata)
    vae_q = scvi$model$SCVI$load_query_data(qdata, vae_ref)
    vae_q$train(max_epochs = as.integer(150), 
        plan_kwargs=dict(weight_decay = as.integer(0)), 
        early_stopping = TRUE) #50 again?
    qlatent = vae_q$get_latent_representation()
    qlatent = as.matrix(qlatent)
    rownames(qlatent) = colnames(sce[[i]])
    saveRDS(qlatent, paste0('snRNAseq/processed_data/transfer_learning/outs/', i, '_qlatent.RDS'))

}

