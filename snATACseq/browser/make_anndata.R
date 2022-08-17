library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

sce_atac <- readRDS("snATACseq/processed_data/umap_gene_score_rm.rds")
anndata <- import("anndata")

adata <- SCE2AnnData(sce_atac)

anndata$AnnData$write(adata, 
  "snATACseq/processed_data/neuro_mat_anndata_scatac.h5ad")
