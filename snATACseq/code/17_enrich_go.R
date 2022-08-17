###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment for GO terms per trajectory                                #
#                                                                             #
#                                                                             #    
###############################################################################

library(qusage)
library(rtracklayer)
library(SingleCellExperiment)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reshape2)
library(dplyr)
library(cowplot)
library(ArchR)
library(chromVAR)

source("snATACseq/R/functions_TF.R")

# 1. load data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data/"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
gene_trends$cluster_trend <- paste0(gene_trends$major_clust, ".", 
                                    gene_trends$gene_trend)

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
gos <- c(`Ion Transport`="GO:0006811")
ind <- which(names(go_terms) %in% gos)
go_ids <- go_terms[[ind]]
go_ids <- rowData(sce)$gene_ids[match(go_ids, rowData(sce)$index)]
go_ids <- go_ids[!is.na(go_ids)]

gene_trends$cluster_trend <- paste0(gene_trends$major_clust, ".",
    gene_trends$major_trend)

genes_cluster_trends <- split(gene_trends$gene_name, gene_trends$cluster_trend)
# remove any combinations with less than 10 genes
genes_cluster_trends <- genes_cluster_trends[!lengths(genes_cluster_trends)<10]

genes_cluster_trends <- lapply(genes_cluster_trends,
    function(y) intersect(go_ids,y))
genes_cluster_trends <- genes_cluster_trends[!lengths(genes_cluster_trends)<10]

# 5. get regions

cre_trends <- mapply(function(X, Y) obtain_peak_ids(X, genes_cluster_trends, 
    Y, cell_types), X=peaks, Y=names(peaks))
cre_trends <- cre_trends[!lengths(cre_trends)==0]

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


saveRDS(cre_trends, file="snATACseq/processed_data/go_ion_gr.RDS")

# 6. load ArchR project

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")
peaks <- peaks[-2] # remove vasculature

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
          "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

# 7. enrichment for each cell type

enriched_gos <- list()
for(i in 1:length(all_gos)){
    
    sc_sub <- sc[sc$Anno1==names(all_gos)[i],]
    
    sc_sub <- addPeakSet(sc_sub, peaks[[names(all_gos)[i]]], force=TRUE)
    sc_sub <- addPeakMatrix(sc_sub, force=TRUE)
    sc_sub <- addMotifAnnotations(ArchRProj = sc_sub, motifSet = "JASPAR2018", 
                          name = "Motif", force=TRUE)
    set.seed(10)
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc_sub, method="ArchR"))
    matches <- getMatches(sc_sub, NULL)
    
    enriched_gos[[i]] <- lapply(all_gos[[i]], function(x)
      compute_enrich_traj(x, bgdPeaks, matches, 
                      names(all_gos)[i]))
}

names(enriched_gos) <- names(all_gos)

enriched_gos <- readRDS("snATACseq/processed_data/go_tf.rds")

for(i in 1:length(enriched_gos)){
  
  for(j in 1:length(enriched_gos[[i]])){
    
    enriched_gos[[i]][[j]]$general_trend <- strsplit(names(enriched_gos[[i]])[j],
      ".", fixed=T)[[1]][3]
  }
  
  tmp <- do.call(rbind, enriched_gos[[i]])
  tmp$pvalue <- 10^(-(tmp$`mlog10p`))
  tmp$adjusted_pvalue <- p.adjust(tmp$pvalue, "fdr")
  
  enriched_gos[[i]] <- tmp
  enriched_gos[[i]]$trajectory <- names(enriched_gos)[i]
  
}

enriched_gos <- do.call(rbind, enriched_gos)

write.csv(enriched_gos, 
  file="snATACseq/processed_data/tf_ion_trends.csv", 
    row.names=FALSE, quote=FALSE)









