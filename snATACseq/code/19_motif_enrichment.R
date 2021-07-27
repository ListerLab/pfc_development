###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment                                                            #
#                                                                             #
#                                                                             #    
###############################################################################

library(ArchR)
library(SummarizedExperiment)
library(ggplot.multistats)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(chromVAR)

# 1. load data

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_annotation")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")

# 2. annotate with peaks and motifs, make peak matrix

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

saveArchRProject(sc,  outputDirectory = "clustering_annotation")

# 3. find background peaks

set.seed(10)
bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc, method="ArchR"))
matches <- getMatches(sc, NULL)

# 4. get motif enrichment and save

compute_enrich <- function(tmp_interest, bgdPeaks, matches){
    
    idx <- match(names(tmp_interest), rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
        as.vector(bgdPeaks[idx,])))
    return(res)

}

tf_regions <- c("tf_gene_trends_gr.RDS", "tf_diseases_gr.RDS")

for(j in 1:length(tf_regions)){

    peaks_of_interest <- readRDS(paste0("processed_data/", tf_regions[j]))
    if(!grepl(".", names(peaks_of_interest[[1]])[1], fixed=T)){
        for(i in 1:length(peaks_of_interest)){
    
            names(peaks_of_interest[[i]]) <- paste0(
                names(peaks_of_interest)[[i]], ".", 
                names(peaks_of_interest[[i]]))
        }
    }
    
    enriched <- lapply(peaks_of_interest, function(x) 
        compute_enrich(x, bgdPeaks, matches))
    enriched <- enriched[!sapply(enriched, function(x) dim(x)[1]==0)]
    for(i in 1:length(enriched)){
        
        enriched[[i]]$cell_type <- names(enriched)[i]
    }
    enriched <- do.call(rbind, enriched)
    enriched$pval <- 10^-(enriched$mlog10p)
    enriched$padjust <- p.adjust(enriched$pval, method="fdr")
    write.table(enriched, file=paste0("snATACseq/processed_data/", 
        gsub("_gr.RDS", ".tsv", tf_regions[j])), quote=F, sep="\t")
}

