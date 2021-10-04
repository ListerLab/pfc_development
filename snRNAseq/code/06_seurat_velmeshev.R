#Seurat sctransform attempt


library(BiocManager)
library(scater)
library(scran)
library(scuttle)
library(reticulate)
library(Seurat)
library(uwot)
library(FNN)
library(stringr)


rawref= readRDS('snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS')
rownames(rawref) = rowData(rawref)$gene_ids
reference = readRDS('snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS') 

#prepare reference - no downsampling

seuraw = CreateSeuratObject(counts = counts(rawref), 
    meta.data = as.data.frame(colData(reference)))
VariableFeatures(seuraw) = subset(rowData(reference)$gene_ids, 
    rowData(reference)$highly_variable == TRUE)
seuraw[['PCA']] = CreateDimReducObject(embeddings = 
    reducedDim(reference, 'PCA')[,1:365], key = 'PCA_')
seuraw[['UMAP']] = CreateDimReducObject(embeddings = 
    reducedDim(reference, 'UMAP'),key = 'UMAP_')
seuraw[['PCA_est']] = CreateDimReducObject(embeddings = 
    reducedDim(reference, 'PCA_est')[,1:365], key = 'PCA_est', 
    loadings = attributes(reducedDim(reference, 'PCA_est'))$rotation)
attributes(seuraw[['PCA_est']]@feature.loadings)$dimnames[[1]] = 
    VariableFeatures(seuraw)
seuraw = SCTransform(seuraw, verbose = TRUE, 
    do.correct.umi = FALSE, residual.features=VariableFeatures(seuraw))


query = readRDS('snRNAseq/processed_data/kriegstein_dataset.RDS') 
# this can be downloaded by https://www.science.org/doi/abs/10.1126/science.aav8130

ribogeneset = read.csv('annotation/ribo_geneset.txt', header = FALSE)$V1
ribout = subset(rowData(query)$ID, rowData(query)$Symbol %in% ribogeneset)

query = query[!( rownames(query) %in% ribout),]
query = query[rowData(query)$Symbol != 'MALAT1',]
query =query[-which(str_detect(rowData(query)$Symbol, '^MT-')),]
query = query[rowSums(counts(query)>0)>=5,]

seuq = CreateSeuratObject(counts = counts(query), 
    meta.data = as.data.frame(colData(query)))
VariableFeatures(seuq) = subset(rownames(seuq), 
    rownames(seuq) %in% VariableFeatures(seuraw))


qanchors = FindTransferAnchors(reference = seuraw, query = seuq, 
    normalization.method='SCT', reduction = 'pcaproject', 
    reference.reduction = "PCA_est", npcs = 365, dims = 1:365, l2.norm=FALSE)


qproj = qanchors@object.list[[1]]@reductions$pcaproject@cell.embeddings[154749:207304,]
saveRDS(qproj, "snRNAseq/processed_data/seurat_umap_kriegstein.RDS") 

