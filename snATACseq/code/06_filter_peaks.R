###############################################################################
#                                                                             #
#                                                                             #
# Filter peaks                                                                #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(GenomicRanges)
library(stringr)
library(magrittr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(rtracklayer)
library(limma)
library(edgeR)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg19)

# 1. read in data

peak_mat_list <- readRDS("snATACseq//processed_data/atac_peak_counts.rds")
cell_types <- names(peak_mat_list)

# 2. filter peaks by widths <300 or >30000 and FPKM

filter_peaks <- function(cell_type, fpkm_min=2,
    sample_min=2, peak_width_range=c(300, 30000)){
    
    cell_dat <- peak_mat_list[[cell_type]]
    
    df <- cell_dat$peak_fpkm
    
    # Filter FPKM
    keep_fpkm <- rowSums(df >= fpkm_min) >= sample_min
    
    widths <- GRanges(rownames(df)) %>% width()
    keep_width <- (widths >= min(peak_width_range)) & 
        (widths <= max(peak_width_range))    
    
    keep <- keep_fpkm & keep_width
    
    pc_keep <- round((sum(keep) / nrow(df)) * 100, digits = 2)
    
    message(str_c("Retained peaks: ", pc_keep, "%"))
    
    df <- df[keep, ]
    
    gr_peaks <- GRanges(rownames(df))
    mcols(gr_peaks) <- df
    
    return(gr_peaks)
}

peak_list <- lapply(cell_types, filter_peaks)

names(peak_list) <- cell_types

#3. save peaks
saveRDS(object = peak_list, 
        file = "snATACseq//processed_data/cell_type_atac_peaks_filtered_gr.rds")
