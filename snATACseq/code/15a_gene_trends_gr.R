###############################################################################
#                                                                             #
#                                                                             #
# Regions for gene trends motif enrichment                                    #
#                                                                             #
#                                                                             #    
###############################################################################

library(GenomicRanges)
library(readxl)
library(ggplot2)
library(SingleCellExperiment)
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# 1. read in data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data/"

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

# 3. find CREs linked to cell types and gene trends

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

genes_cluster_trends <- split(gene_trends$gene_name, gene_trends$cluster_trend)
# remove cell type and gene trends combinations with less than 10 CREs
genes_cluster_trends <- genes_cluster_trends[!lengths(genes_cluster_trends)<10]
genes_cluster_trends$`L2-3_CUX2.0` <- genes_cluster_trends$`L2-3_CUX2.0`[!is.na(
    genes_cluster_trends$`L2-3_CUX2.0`)]

obtain_peak_ids <- function(peak_tmp, genes_cluster_trends, cluster_tmp, 
                            cell_types){
    
    cluster_tmp <- cell_types[cell_types[,2] %in% cluster_tmp,1]
    ind <- sapply(cluster_tmp, function(x) 
        grepl(paste0(x, "."), names(genes_cluster_trends)))
    if(class(ind)[[1]]=="matrix"){
        ind <- apply(ind, 1, function(x) any(x))
    }
    tmp <- genes_cluster_trends[ind]
    ind <- lapply(tmp, function(x) unlist(sapply(x, function(y) 
        grep(y, peak_tmp$CRE))))
    split_tmp <- lapply(ind, function(x) peak_tmp[x])
    split_tmp
}

peaks <- peaks[-2]

cre_trends <- mapply(function(X, Y) obtain_peak_ids(X, genes_cluster_trends, 
    Y, cell_types), X=peaks, Y=names(peaks))

cre_trends <- lapply(cre_trends, function(x) x[!is.na(names(x))])
for(i in 1:length(cre_trends)){
    
    for(j in 1:length(cre_trends[[i]])){
        
        names(cre_trends[[i]][[j]]) <- paste0(names(cre_trends)[i], ".",
          names(cre_trends[[i]][[j]]))
    }
    names(cre_trends[[i]]) <- paste0(names(cre_trends)[i], ".", 
        names(cre_trends[[i]]))
}


for(i in 1:length(cre_trends)){
    
    print(i)
    tmp_names <- names(cre_trends[[i]])
    tmp_names <- sapply(tmp_names, function(x) 
        strsplit(x, ".", fixed=TRUE)[[1]][3]) 
    if(any(duplicated(tmp_names))){
        
        split_names <- split(names(tmp_names), tmp_names)
        tmp <- lapply(split_names, function(x) 
            unlist(as(cre_trends[[i]][x], "GRangesList")))
        names(tmp) <- paste0(names(cre_trends)[i], ".", names(cre_trends)[i],
            ".", names(tmp))
        cre_trends[[i]] <- tmp
    }
    
}

saveRDS(cre_trends, file="snATACseq/processed_data/tf_gene_trends_trajectories_gr.RDS")
cre_trends <- readRDS("snATACseq/processed_data/tf_gene_trends_trajectories_gr.RDS")
## 4. make downwards 

down_up <- c(`0`="up" , `1`="down", `2`="down", `3`="up", `4`="down", `5`="interup", 
             `6`="up", `7`="up", `8`="interdown", `9`="down", `10`="down", `11`="interup", 
             `12`="interdown", `13`="interup")

cre_general <- list()

for(i in 1:length(cre_trends)){
    
    tmp_names <- names(cre_trends[[i]])
    tmp_names <- sapply(tmp_names, function(x) 
        strsplit(x, ".", fixed=TRUE)[[1]][3])
    tmp_down <- down_up[tmp_names]
    
    split_names <- split(names(tmp_names), tmp_down)
    tmp <- lapply(split_names, function(x) 
        unlist(as(cre_trends[[i]][x], "GRangesList")))
    names(tmp) <- paste0(names(cre_trends)[i], ".", names(cre_trends)[i],
                         ".", names(tmp))
    cre_general[[i]] <- tmp
    
}

names(cre_general) <- names(cre_trends)
saveRDS(cre_general, file="snATACseq/processed_data/tf_general_trends_trajectories_gr.RDS")    
    