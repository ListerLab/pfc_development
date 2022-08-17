###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment in CREs of major trends                                    #
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

# 1. read in data

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

for(i in 1:length(cre_trends)){
  
  for(j in 1:length(cre_trends[[i]])){
    cell_type <- colnames(cre_trends[[i]][[j]]@elementMetadata)[1]
    cell_type <- strsplit(cell_type, "_RL")[[1]][1]
    cre_trends[[i]][[j]]@elementMetadata <- DataFrame(new_name=
        paste0(cell_type, ".", names(cre_trends[[i]][[j]])))
    
  }
         
}

cre_trends <- unlist(cre_trends, recursive = F, use.names = T)

# 4. save regions for each gene trend

peaks_of_interest <- lapply(unique(gene_trends$gene_trend), function(x) 
    unlist(cre_trends[grepl(paste0(x, "$"), names(cre_trends))]))
peaks_of_interest <- lapply(peaks_of_interest, function(x) 
    unlist(as(x, "GRangesList"), use.names = F))
names(peaks_of_interest) <- unique(gene_trends$gene_trend)

if(!grepl(".", names(peaks_of_interest[[1]])[1], fixed=T)){
  for(i in 1:length(peaks_of_interest)){
    
    names(peaks_of_interest[[i]]) <- paste0(
    names(peaks_of_interest)[[i]], ".", 
    names(peaks_of_interest[[i]]))
      
  }
}

saveRDS(peaks_of_interest, file="snATACseq/processed_data/gene_trends_gr.rds")

peaks_of_interest <- readRDS("snATACseq/processed_data/gene_trends_gr.rds")
peaks_of_interest <- lapply(peaks_of_interest, function(x)
  unique(x))
 
# 5. load ArchR project

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")
peaks <- peaks[-2]

# 6. annotate with peaks and motifs, make peak matrix

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
        "nearestGene", "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

peaks <- unlist(as(peaks, "GRangesList"))

sc <- addPeakSet(sc, peaks, force=TRUE)
sc <- addPeakMatrix(sc, force=TRUE)
sc <- addMotifAnnotations(ArchRProj = sc, motifSet = "JASPAR2018", 
    name = "Motif", force=TRUE)

rm(peaks)
gc()
gc()

# 7. find background peaks

set.seed(10)
bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc, method="ArchR"))
matches <- getMatches(sc, NULL)

# 8. get motif enrichment and save
    
enriched <- lapply(peaks_of_interest, function(x) 
  compute_enrich(x, bgdPeaks, matches))
enriched <- enriched[!sapply(enriched, function(x) dim(x)[1]==0)]
for(i in 1:length(enriched)){
        
  enriched[[i]]$cell_type <- names(enriched)[i]
}
enriched <- do.call(rbind, enriched)
enriched$pval <- 10^-(enriched$mlog10p)
enriched$padjust <- p.adjust(enriched$pval, method="fdr")
write.table(enriched, file=
    paste0("snATACseq/processed_data/gene_trends_tf.tsv"), quote=F, sep="\t")


