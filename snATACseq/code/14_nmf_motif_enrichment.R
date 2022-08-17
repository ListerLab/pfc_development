###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment NMF regions                                                #
#                                                                             #
#                                                                             #    
###############################################################################

library(ArchR)
library(SummarizedExperiment)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(chromVAR)

# 1. load data

source("snATACseq/R/functions_TF.R")

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

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

saveArchRProject(sc,  outputDirectory = "snATACseq/clustering_final")

# 3. find background peaks

set.seed(10)
bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc, method="ArchR"))
matches <- getMatches(sc, NULL)

# 4. motif enrichment distinct

peaks_nmf <- readRDS("snATACseq/processed_data/peaks_nmf_distinct.RDS")
peaks_nmf <- peaks_nmf[!lengths(peaks_nmf)<20]

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
