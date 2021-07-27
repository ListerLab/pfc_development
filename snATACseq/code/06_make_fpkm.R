###############################################################################
#                                                                             #
#                                                                             #
# Make ATAC peak normalied insertion matricies                                #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(bsseq)
library(GenomicRanges)
library(stringr)
library(magrittr)
library(data.table)
library(GenomicFeatures)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(rtracklayer)
library(limma)
library(edgeR)
library(parallel)
#library(readxl)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

# 1. list of files
peak_list <- list.files(path = "peaks", pattern = ".bed",
                        full.names = TRUE)
peak_list <- peak_list[!grepl(pattern = "_RL", x = peak_list)]

cell_types <- c("Astro", "L2_3", "L5_6", "Micro",
                "CGE", "L4", "MGE", "Oligo", "Vas")

sample_beds <- list.files(path = "snATACseq/bedfiles/",
                          pattern = "_blacklistrm.bed", full.names = TRUE)
sample_beds <- sample_beds[grepl(pattern = "_RL", x = sample_beds)]

# 2. make normalized peak overlap

count_peak_olaps <- function(cell_type){
    
    cell_peak_fl <- peak_list[grepl(pattern = cell_type, x = peak_list)]
    stopifnot(length(cell_peak_fl) == 1)
    
    cell_bed_files <- sample_beds[grepl(pattern = cell_type, x = sample_beds)]
    
    bed_to_gr <- function(bed){
        peaks <- fread(bed)
        gr <- GRanges(seqnames = peaks$V1,
                      ranges = IRanges(start = peaks$V2, end = peaks$V3))
        return(gr)
    }
    
    peaks_gr <- bed_to_gr(cell_peak_fl)
    
    # Subset for autosomes
    peaks_gr <- peaks_gr[seqnames(peaks_gr) %in% str_c("chr", 1:22)]
    
    # Count overlaps for each sample.
    # Minoverlap set to 24 to ensure insertion is within peak
    count_peak_olaps <- function(x){
        sample_gr <- bed_to_gr(cell_bed_files[x])
        
        sample_id <- basename(cell_bed_files[x]) %>%
            str_remove("_blacklistrm.bed")
        
        olaps <- GenomicRanges::countOverlaps(query = peaks_gr,
                                              subject = sample_gr,
                                              minoverlap = 24,
                                              ignore.strand=TRUE)
        df <- data.frame(olaps)
        colnames(df) <- sample_id
        return(df)
    }
    
    peak_counts <- mclapply(1:length(cell_bed_files),
                            count_peak_olaps, mc.cores=24)
    peak_counts <- do.call(cbind, peak_counts)
    
    # Add the peak loci to matrix
    loci <- str_c(seqnames(peaks_gr), ":", start(peaks_gr), "-", end(peaks_gr))
    rownames(peak_counts) <- loci
    
    # Get peak widths for length normalisation
    peak_widths <- width(peaks_gr)
    
    peak_fpkm <- edgeR::rpkm(peak_counts, gene.length=peak_widths)
    
    out <- list(peak_counts=peak_counts, peak_fpkm=peak_fpkm)
    
    return(out)
}

peak_mat_list <- lapply(cell_types, count_peak_olaps)
names(peak_mat_list) <- cell_types

# 3. save peak counts as Rds object
saveRDS(peak_mat_list, "snATACseq/processed_data/atac_peak_counts.Rds")