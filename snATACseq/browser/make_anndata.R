library(SingleCellExperiment)
library(zellkonverter)
library(reticulate)

dat <- readRDS("snATACseq/processed_data/UMAP_GeneScore_rm.rds")
heanndata <- import("anndata")

col_data <- DataFrame(barcode=rownames(dat$coldata),
    sample=dat$coldata$Sample, stage=dat$coldata$Stage,
    cell_type=dat$coldata$Anno1, arcsin_age=dat$coldata$arcsin_ages)
rownames(col_data) <- rownames(dat$coldata)

aa <- as(dat$logcounts, "dgCMatrix")
aa@Dimnames <- list(NULL, NULL)

umap <- as.matrix(dat$df_umap)
colnames(umap) <- c("UMAP1", "UMAP2")

sce <- SingleCellExperiment(SimpleList(logcounts=aa),
            colData=col_data)
reducedDim(sce, "UMAP") <- umap
rowData(sce) <- DataFrame(Gene=rownames(dat$logcounts))
adata <- SCE2AnnData(sce)

anndata$AnnData$write(adata, "processed_data/neuro_mat_anndata_scatac.h5ad")

