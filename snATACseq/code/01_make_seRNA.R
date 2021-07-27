###############################################################################
#                                                                             #
#                                                                             #
# Make SummarizedExperiment from snRNAseq                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(SummarizedExperiment)
library(rhdf5)
library(Matrix)
library(rtracklayer)

file <- "snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"
# read in gene locations
genes <- rtracklayer::import("annotation/hg19_genes.bed")

# read in file
sce <- readRDS(paste0(file))

# keep only genes with known locations 
genes_sce <- rowData(sce)$index
ind <- genes_sce %in% genes@elementMetadata$name 
sce <- sce[ind,]
genes_sce <- genes_sce [ind]
ind <- match(genes_sce, genes@elementMetadata$name)
genes <- genes[ind,]

# save as Summarized Experiment with locations
seRNA <- SummarizedExperiment(assays=SimpleList(counts=counts(sce)), 
     rowRanges=genes, colData=colData(sce))

saveRDS(seRNA, file="snRNAseq/processed_data/snRNAseq_SE.RDS")

