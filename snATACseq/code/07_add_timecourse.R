###############################################################################
#                                                                             #
#                                                                             #
# Peak annotation and temporal timecourse                                     #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(rtracklayer)
library(parallel)
library(readxl)
library(data.table)
library(stringr)
library(ArchR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

source("snATACseq/R/functions_timecourse.R")

# 1. read in data and annotations

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

cell_types <- c("Astro", "Vas", "L5_6", "L4", "L2_3", "Micro", "Oligo",
                "MGE_der", "CGE_der")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peak_set_all <- readRDS("snATACseq/processed_data/cell_type_atac_peaks_filtered_gr.rds")
names(peak_set_all) <- gsub("MGE", "MGE_der", names(peak_set_all))
names(peak_set_all) <- gsub("CGE", "CGE_der", names(peak_set_all))

all_files <- list.files("snATACseq/peaks")
all_files <- all_files[grepl("_peaks_final_merge.bed", all_files)]
all_files <- paste0("snATACseq/peaks/", all_files)

# 2. add information of sample-specific accessibility

make_access_matrix <- function(cell_type){
    
    peak_set <- peak_set_all[[cell_type]]
    names(peak_set) <- paste0("id_", 1:length(peak_set))
        
    files_for_import <- all_files[grepl(cell_type, all_files)]
    files_for_import <- setdiff(files_for_import, paste0("peaks/", cell_type, 
            "_peaks_final_merge.bed"))
    files_for_import <- files_for_import[grepl("RL", files_for_import)]
    
    all_peaks <- lapply(files_for_import, function(x) 
        import(x, format="bed"))
    all_peaks <- GRangesList(all_peaks)
    names(all_peaks) <- str_extract(files_for_import, "RL[0-9]{4}")
    all_peaks <- all_peaks[!lengths(all_peaks)==0]
    
    peak_set <- annotate_peaks(peak_set, txdb)
    overlap_index <- lapply(all_peaks, function(x) 
        findOverlaps(peak_set, x))
    # find width of overlap over all intersecting peaks
    overlap_length <- lapply(1:length(overlap_index), function(x)
        data.frame(index_query=overlap_index[[x]]@from, 
            width=overlapsRanges(ranges(peak_set), 
            ranges(all_peaks[[x]]), overlap_index[[x]])@width))
    overlap_length <- lapply(overlap_length, function(x)
        split(x, x$index_query)
    )
    overlap_length <- lapply(overlap_length, function(x) t(sapply(x, function(y)
        data.frame(index_query=unique(y$index_query), width=sum(y[,2]))
    )))
    # see whether overlap is more than half the peak
    overlap_length <- lapply(overlap_length, function(y)
        data.frame(index_query=unlist(y[,1]), 
            width=unlist(y[,2])/ranges(peak_set[unlist(y[,1])])@width>0.5)
    )
    access_mat <- matrix(FALSE, nrow=length(peak_set), ncol=length(all_peaks)) 
    for(j in 1:length(overlap_length)){
        access_mat[overlap_length[[j]][,1], j] <- overlap_length[[j]][,2]
    }
    print(sum(rowSums(access_mat)==FALSE)/dim(access_mat)[1])
    colnames(access_mat) <- names(all_peaks)
    peak_set@elementMetadata <- cbind(peak_set@elementMetadata, access_mat)
    return(peak_set)
}

processed_data <- lapply(cell_types, function(x) make_access_matrix(x))
names(processed_data) <- cell_types

# 3. add stage accessibility 

for(i in 1:length(processed_data)){

    tmp <- processed_data[[i]]@elementMetadata
    tmp <- tmp[,colnames(tmp) %in% atac_samples$Sample]
    tmp <- as.matrix(tmp)
    tmp <- apply(tmp, 2, function(x) as.numeric(x))
    stages_tmp <- atac_samples$Stage[atac_samples$Sample %in% colnames(tmp)]
    tmp <- sapply(stages, function(x) 
    rowSums(tmp[, stages_tmp==x, drop=FALSE])>0)
    processed_data[[i]]@elementMetadata <- cbind(
        processed_data[[i]]@elementMetadata, tmp)
  
}

# 3. save data
saveRDS(processed_data, 
    file="snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

# 4. add peaks to scArchR object

sc <- loadArchRProject("snATACseq/clustering_final")
sc$Anno1 <- gsub("/", "_", sc$Anno1)
sc$Anno1 <- gsub(" ", "_", sc$Anno1)

peaks <- processed_data

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(
        "nearestGene", "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

peaks <- unlist(as(peaks, "GRangesList"))
peaks <- reduce(peaks)

sce_atac <- readRDS("snATACseq/processed_data/umap_gene_score_rm.rds")

sc <- addPeakSet(sc, peaks, force=TRUE)
sc <- addPeakMatrix(sc)
peak_mat <- getMatrixFromProject(sc, "PeakMatrix")

index <- match(colnames(sce_atac), colnames(peak_mat))

sce_atac$FRIP <- sc$FRIP
sce_atac$ReadsInPeaks <- sc$ReadsInPeaks

saveRDS(sce_atac, file="snATACseq/processed_data/umap_gene_score_rm.rds")



  
  

