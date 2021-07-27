###############################################################################
#                                                                             #
#                                                                             #
# Regions for GO pathways and trend motif enrichment                          #
#                                                                             #
#                                                                             #    
###############################################################################

library(qusage)
library(rtracklayer)
library(SingleCellExperiment)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reshape2)
library(cba)
library(cowplot)
library(viridis)
library(ggplot2)
library(dplyr)
library(jaccard)
library(cowplot)

# 1. load data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
gene_trends$cluster_trend <- paste0(gene_trends$major_clust, ".", 
                                    gene_trends$gene_trend)

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")

sce <- readRDS(paste0(paths,
    "/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))
gene_trends$gene_name <- rowData(sce)$gene_ids[match(gene_trends$gene_name, 
                                                     rownames(sce))]

# 2. remove all promoters regions

peaks <- lapply(peaks, function(x) {
    hits <- findOverlaps(x, promoters)
    ind <- setdiff(1:length(x), unique(hits@from))
    x[ind]
})

# 3. get genes in GO terms

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
                       "OPC", "Oligo",
                       "Micro", "Micro"
), ncol=2, byrow=T)

go_terms <- read.gmt("annotation/hsapiens.GO:BP.name.gmt")
gos <- c(Ensheatment="GO:0007272", `Cell Cycle`="GO:0007049", 
    `Neuron Migration`="GO:0001764", `Synapse Organization`="GO:0050808",
    `Cytoskeletal Organization`="GO:0051493", Learning="GO:0007612",
    `Cell Death`="GO:0008219", `Ion Transport`="GO:0006811")
ind <- which(names(go_terms) %in% gos)
go_ids <- go_terms[ind]
go_ids <- lapply(go_ids, function(x)
    rowData(sce)$gene_ids[match(x, rowData(sce)$index)])
go_ids <- lapply(go_ids, function(x) x[!is.na(x)])

down_up <- c(`0`="up" , `1`="down", `2`="down", `3`="up", `4`="down", `5`="interup", 
    `6`="up", `7`="up", `8`="interdown", `9`="down", `10`="down", `11`="interup", 
    `12`="interdown", `13`="interup")
gene_trends$up_down <- down_up[match(gene_trends$gene_trend,names(down_up))]
gene_trends$cluster <- cell_types[match(gene_trends$major_clust, cell_types[,1]),2]
gene_trends$cluster_up_down <- paste0(gene_trends$cluster, ".",
    gene_trends$up_down)

genes_cluster_trends <- split(gene_trends$gene_name, gene_trends$cluster_up_down)
# remove any combinations with less than 10 genes
genes_cluster_trends <- genes_cluster_trends[!lengths(genes_cluster_trends)<10]

genes_cluster_trends <- lapply(go_ids, function(x) lapply(genes_cluster_trends,
    function(y) intersect(x,y)))
genes_cluster_trends <- lapply(genes_cluster_trends, function(x)
    x[!lengths(x)<10])

# 5. get regions

obtain_peak_ids <- function(peak_tmp, genes_cluster_trends, cluster_tmp){
    
    ind <- sapply(cluster_tmp, function(x) 
        grepl(paste0(x, "."), names(genes_cluster_trends)))
    tmp <- genes_cluster_trends[ind]
    ind <- lapply(tmp, function(x) unlist(sapply(x, function(y) 
        grep(y, peak_tmp$CRE))))
    split_tmp <- lapply(ind, function(x) peak_tmp[x])
    split_tmp
}

peaks <- peaks[-2]

cre_trends <- list()

for(i in 1:length(genes_cluster_trends)){

cre_trends[[i]] <- mapply(function(X, Y) obtain_peak_ids(X, 
    genes_cluster_trends[[i]], Y), X=peaks, Y=names(peaks))
cre_trends[[i]] <- unlist(cre_trends[[i]], recursive = FALSE) 

}

names(cre_trends) <- names(genes_cluster_trends)
cre_trends <- unlist(cre_trends, recursive = FALSE)

# 6. save regions

saveRDS(cre_trends, file="snATACseq/processed_data/go_gr.RDS")







