###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment NMF regions                                                #
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
         "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

peaks <- unlist(as(peaks, "GRangesList"))

sc <- addPeakSet(sc, peaks, force=TRUE)
sc <- addPeakMatrix(sc, force=TRUE)
sc <- addMotifAnnotations(ArchRProj = sc, motifSet = "JASPAR2018", 
            name = "Motif", force=TRUE)

saveArchRProject(sc,  outputDirectory = "snATACseq/clustering_annotation")

# 3. find background peaks

set.seed(10)
bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc, method="ArchR"))
matches <- getMatches(sc, NULL)

# 4. motif enrichment distinct

peaks_nmf <- readRDS("snATACseq/processed_data/peaks_nmf_distinct.RDS")
peaks_nmf <- peaks_nmf[!lengths(peaks_nmf)<20]

compute_enrich <- function(tmp_interest, bgdPeaks, matches){
    
    idx <- match(tmp_interest, rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
        as.vector(bgdPeaks[idx,])))
    return(res)
    
}

enriched <- lapply(peaks_nmf, function(x) 
    compute_enrich(x, bgdPeaks, matches))
enriched <- enriched[!sapply(enriched, function(x) dim(x)[1]==0)]
for(i in 1:length(enriched)){
    
    enriched[[i]]$cell_type <- names(enriched)[i]
}
enriched <- do.call(rbind, enriched)
enriched$pval <- 10^-(enriched$mlog10p)
enriched$padjust <- p.adjust(enriched$pval, method="fdr")
write.table(enriched, file=paste0("snATACseq/processed_data/TF_nmf.tsv"), 
    quote=F, sep="\t")

# 5. motif enrichment all

peaks_max_mod <- readRDS("snATACseq/processed_data/peaks_max_mod.RDS")

enriched <- lapply(peaks_max_mod, function(x) 
    compute_enrich(x, bgdPeaks, matches))
enriched <- enriched[!sapply(enriched, function(x) dim(x)[1]==0)]
for(i in 1:length(enriched)){
    
    enriched[[i]]$cell_type <- names(enriched)[i]
}
enriched <- do.call(rbind, enriched)
enriched$pval <- 10^-(enriched$mlog10p)
enriched$padjust <- p.adjust(enriched$pval, method="fdr")
write.table(enriched, file=paste0("snATACseq/processed_data/TF_nmf_max_mod.tsv"), 
            quote=F, sep="\t")
