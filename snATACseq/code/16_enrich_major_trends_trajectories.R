###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment in CREs of major trends per trajectory                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(GenomicRanges)
library(readxl)
library(ggplot2)
library(SingleCellExperiment)
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ArchR)

source("snATACseq/R/functions_TF.R")

# 1. read in data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data/"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
gene_trends$cluster_trend <- paste0(gene_trends$major_clust, ".", 
    gene_trends$major_trend)

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

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

peaks <- peaks[-2]

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

peaks_of_interest <- cre_trends

saveRDS(peaks_of_interest, 
  file="snATACseq/processed_data/major_trends_trajectories_gr.rds")

# 4. load ArchR project

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("clustering_final")
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")
peaks <- peaks[-2]

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
          "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

# 2. enrichment for each cell type

enriched_list <- list()

for(i in 1:length(peaks)){
    
    traj <- names(peaks)[i]
    sc_sub <- sc[sc$Anno1==names(peaks)[i],]
    
    sc_sub <- addPeakSet(sc_sub, peaks[[i]], force=TRUE)
    sc_sub <- addPeakMatrix(sc_sub, force=TRUE)
    sc_sub <- addMotifAnnotations(ArchRProj = sc_sub, motifSet = "JASPAR2018", 
                          name = "Motif", force=TRUE)
    set.seed(10)
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc_sub, method="ArchR"))
    matches <- getMatches(sc_sub, NULL)
    
    enriched_list[[i]] <- lapply(peaks_of_interest[[traj]], function(x) 
        compute_enrich_traj(x, bgdPeaks, matches))
}

# 3. save motifs

names(enriched_list) <- names(peaks)

for(i in 1:length(enriched_list)){
    
    for(j in 1:length(enriched_list[[i]])){
        
        enriched_list[[i]][[j]]$trajectory <- names(enriched_list[[i]])[j]
        enriched_list[[i]][[j]]$pvalue <- 10^(-enriched_list[[i]][[j]]$mlog10p)
    }
    
    enriched_list[[i]] <- do.call(rbind, enriched_list[[i]])
    enriched_list[[i]]$adjusted_pvalue <- p.adjust(enriched_list[[i]]$pvalue, 
      "fdr")
}

test <- lapply(enriched_list, function(x) x)
test <- do.call(rbind, test)

test$feature <- sapply(test$feature, function(x) strsplit(x, "_")[[1]][1])
test$general_trend <- sapply(test$trajectory, function(x) 
    strsplit(x, ".", fixed=TRUE)[[1]][3])
test$trajectory <- sapply(test$trajectory, function(x) 
    strsplit(x, ".", fixed=TRUE)[[1]][1])

write.csv(test, 
  file="snATACseq/processed_data/tf_general_trends_trajectories.csv", 
    row.names=FALSE, quote=FALSE)


