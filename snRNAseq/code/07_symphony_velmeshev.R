#Symphony

library(BiocManager)
library(scater)
library(scran)
library(scuttle)
library(reticulate)
library(Seurat)
library(uwot)
library(symphony)

reference = readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS") 
rawref =readRDS('snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS')
rownames(rawref) = rowData(rawref)$gene_ids
bycol = 1000/colSums(counts(rawref))            #ensures downsample target is 1000
set.seed(100)
ref_exp = downsampleMatrix(x = counts(rawref), prop = bycol, bycol = TRUE)
ref_metadata = reference@colData[c(71,69,68,66,62,59)]
initembedding = reducedDim(reference, "UMAP")

hvg = subset(rowData(reference)$gene_ids, 
    rowData(reference)$highly_variable == TRUE)
set.seed(0)

new_sym_ref = symphony::buildReference(
  ref_exp,
  ref_metadata,
  K = 100,
  verbose = TRUE,
  do_umap = TRUE,
  do_normalize = TRUE,
  vargenes_method = 'vst',
  topn = 0,
  d = 365,
  additional_genes = hvg,
  save_uwot_path = '.symphony_trial2'
)

krieg = readRDS('datasets/kriegstein_dataset.RDS') 
#this data is available https://www.science.org/doi/full/10.1126/science.aav8130

ribogeneset = read.csv('annotation/ribo_geneset.txt', header = FALSE)$V1 

ribout = subset(rowData(krieg)$ID, rowData(krieg)$Symbol %in% ribogeneset)

krieg = krieg[!( rownames(krieg) %in% ribout),]
krieg = krieg[rowData(krieg)$Symbol != 'MALAT1',]
krieg =krieg[-which(str_detect(rowData(krieg)$Symbol, '^MT-')),]
krieg = krieg[rowSums(counts(krieg)>0)>=5,]
bycol = 1000/colSums(counts(krieg))            #ensures downsample target is 1000
set.seed(100)
krieg_exp = counts(krieg)
krieg_metadata = krieg@colData[c(1,2,3,8,13,14,15)]

query = mapQuery(krieg_exp, krieg_metadata, new_sym_ref, 
    do_normalize = TRUE, do_umap = TRUE)

reference_umap = umap(t(new_sym_ref$Z_corr), 
        n_neighbors = 25, ret_model = TRUE, init = initembedding, 
        a = 0.58303, b = 1.334167, spread = 1, min_dist = 0.5, 
        verbose = TRUE, n_epochs = 0)
qumap = umap_transform(t(query$Z), reference_umap, n_epochs = 50)

saveRDS(qumap, file="snRNAseq/processed_data/symphony_umap_kriegstein.RDS")






