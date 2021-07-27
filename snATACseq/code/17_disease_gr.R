###############################################################################
#                                                                             #
#                                                                             #
# Regions for disease gene trends motif enrichment                            #
#                                                                             #
#                                                                             #    
###############################################################################

library(GenomicRanges)
library(readxl)
library(ggplot2)
library(SingleCellExperiment)
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(disgenet2r)

# 1. read in data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data"

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")

for(i in 1:length(peaks)){
    
    names(peaks[[i]]) <- paste0(names(peaks)[i], ".", names(peaks[[i]]))
}

sce <- readRDS(paste0(paths,
    "/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))


# 2. remove all promoters regions

peaks <- lapply(peaks, function(x) {
    hits <- findOverlaps(x, promoters)
    ind <- setdiff(1:length(x), unique(hits@from))
    x[ind]
})

# 3. load disease genes
all_files <- list.files(paths, "/disease_enrichment_genes/")
all_diseases <- lapply(all_files, function(x) read.csv(paste0(paths,
    "/disease_enrichment_genes/", x)))
all_diseases <- lapply(all_diseases, function(x) x[,2])
names(all_diseases) <- gsub(".csv", "", all_files)

all_diseases <- lapply(all_diseases, function(x)
    rowData(sce)$gene_ids[match(x, rowData(sce)$index)])

cell_types <- matrix(c("Astro", "Astro",
                       "L2-3_CUX2", "L2_3",
                       "L5-6_THEMIS", "L5_6",
                       "L5-6_TLE4", "L5_6",
                       "L4_RORB", "L4",
                       "VIP", "CGE_der",
                       "LAMP5_CA1", "CGE_der",
                       "ID2", "CGE_der",
                       "SST", "MGE_der",
                       "PV", "MGE_der",
                       "PV_SCUBE3", "MGE_der",
                       "Oligo", "Oligo",
                       "Micro", "Micro"
), ncol=2, byrow=T)

# 4. get peaks for disease, gene trend, cell type combinations

obtain_peak_ids <- function(disease_tmp, disease_name_tmp, peaks, 
                            cell_types){
    
    ind <- which(str_count(disease_name_tmp, cell_types[,1])==1)
    tmp <- peaks[[cell_types[ind[1],2]]]
    ind <- unlist(sapply(disease_tmp, function(y) 
        grep(y, tmp$CRE)))
    tmp[ind]
}

cre_diseases <- mapply(function(X,Y) obtain_peak_ids(X,Y, peaks, cell_types),
       X=all_diseases, Y=names(all_diseases))

saveRDS(cre_diseases, file="snATACseq/processed_data/tf_diseases_gr.RDS")



