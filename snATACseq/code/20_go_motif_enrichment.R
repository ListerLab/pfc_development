###############################################################################
#                                                                             #
#                                                                             #
# Motif enrichment GO pathways                                                #
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
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
peaks <- peaks[-2] # remove vasculature

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
          "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

all_gos <- readRDS("snATACseq/processed_data/go_gr_neurotransmitter.RDS")

# 2. enrichment for each cell type

compute_enrich <- function(tmp_interest, bgdPeaks, matches){
    
    idx <- match(names(tmp_interest), rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
          as.vector(bgdPeaks[idx,])))
    return(res)
    
}


enriched_gos <- list()
for(i in 1:length(peaks)){
    
    sc_sub <- sc[sc$Anno1==names(peaks)[i],]
    
    sc_sub <- addPeakSet(sc_sub, peaks[[i]], force=TRUE)
    sc_sub <- addPeakMatrix(sc_sub, force=TRUE)
    sc_sub <- addMotifAnnotations(ArchRProj = sc_sub, motifSet = "JASPAR2018", 
                          name = "Motif", force=TRUE)
    set.seed(10)
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc_sub, method="ArchR"))
    matches <- getMatches(sc_sub, NULL)
    
    gos_tmp <- all_gos[grepl(names(peaks)[i], names(all_gos))]
    enriched_gos[[i]] <- lapply(gos_tmp, function(x) 
        compute_enrich(x, bgdPeaks, matches))
}

# 3. save motifs

names(enriched_gos) <- names(peaks)
saveRDS(enriched_gos, file="processed_data/Gos_tf_new_no_proms.RDS")

