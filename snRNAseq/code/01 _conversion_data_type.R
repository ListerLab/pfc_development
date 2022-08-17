###############################################################################
#                                                                             #
#                                                                             #
#  Make SingleCellExperiment from anndata                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(SingleCellExperiment)
library(zellkonverter)

link <- "" #insert path to 
file_name <- "2020-12-18_whole-tissue_post-restaged-GABA-clustering.h5ad"

sce <- readH5AD(paste0(link, file_name))
rowData(sce)$index <- rownames(sce)
names(sce@int_colData$reducedDims@listData) <- c("PCA", "UMAP")

names(sce@assays@data@listData) <- "logcounts"

outlink <- 'snRNAseq/processed_data/'

saveRDS(sce, file=paste0(outlink, gsub(".h5ad",".RDS", file_name)))

file_name <- "2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.h5ad"

sce <- readH5AD(paste0(link, file_name))
rowData(sce)$index <- rownames(sce)
names(sce@int_colData$reducedDims@listData) <- c("PCA", "UMAP", "UMAT")

names(sce@assays@data@listData) <- "counts"

outlink <- 'snRNAseq/processed_data/'

saveRDS(sce, file=paste0(outlink, gsub(".h5ad",".RDS", file_name)))

