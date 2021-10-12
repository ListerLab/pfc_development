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

sc <- loadArchRProject("clustering_annotation")
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "processed_data/cell_type_atac_peaks_filtered_anno_")
peaks <- peaks[-2] # remove vasculature

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
          "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

all_regions <- readRDS("processed_data/tf_general_trends_trajectories_gr.RDS")

# 2. enrichment for each cell type

compute_enrich <- function(tmp_interest, bgdPeaks, matches, traj){
    
    names(tmp_interest) <- sapply(names(tmp_interest), function(x)
        strsplit(x, "id_", fixed=TRUE)[[1]][2])
    names(tmp_interest) <- paste0("id_", names(tmp_interest))
    idx <- match(names(tmp_interest), rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
          as.vector(bgdPeaks[idx,])))
    return(res)
    
}


enriched_gos <- list()
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
    
    enriched_gos[[i]] <- lapply(all_regions[[traj]], function(x) 
        compute_enrich(x, bgdPeaks, matches))
}

# 3. save motifs

names(enriched_gos) <- names(peaks)

for(i in 1:length(enriched_gos)){
    
    for(j in 1:length(enriched_gos[[i]])){
        
        enriched_gos[[i]][[j]]$trajectory <- names(enriched_gos[[i]])[j]
        enriched_gos[[i]][[j]]$pvalue <- 10^(-enriched_gos[[i]][[j]]$mlog10p)
    }
    
    enriched_gos[[i]] <- do.call(rbind, enriched_gos[[i]])
    enriched_gos[[i]]$adjusted_pvalue <- p.adjust(enriched_gos[[i]]$pvalue, "fdr")
}

test <- lapply(enriched_gos, function(x) x[x$adjusted_pvalue<0.05,])
test <- do.call(rbind, test)

test$feature <- sapply(test$feature, function(x) strsplit(x, "_")[[1]][1])
test$general_trend <- sapply(test$trajectory, function(x) 
    strsplit(x, ".", fixed=TRUE)[[1]][3])
test$trajectory <- sapply(test$trajectory, function(x) 
    strsplit(x, ".", fixed=TRUE)[[1]][1])

write.csv(test, file="tf_general_trends_trajectories_sig.csv", row.names=FALSE,
          quote=FALSE)

saveRDS(enriched_gos, file="processed_data/tf_general_trends_trajectories.RDS")

