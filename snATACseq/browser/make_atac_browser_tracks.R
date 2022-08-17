library(GenomicRanges)
library(stringr)
library(magrittr)
library(data.table)
library(ggplot2)
library(reshape2)
library(rtracklayer)
library(edgeR)
library(BSgenome.Hsapiens.UCSC.hg19)


setwd("snATACseq")

### Make the Tn5 insertion bigwig files


# List the Tn5 insertion bed files
bed_fls <- list.files(path = "bedfiles/stages/", pattern = ".bed",
                      full.names = TRUE)

bed_to_cov <- function(fl, normalise=TRUE){
    
    dat <- fread(fl, select = c(1:3)) 
    
    # Make GRanges for each insertion point
    gr <- GRanges(seqnames = dat$V1,
                  ranges = IRanges(start = dat$V2, end = dat$V3))
    rm(dat)

    gr <- sort(gr) # Coverage runs faster with sorted ranges
    
    cov <- GenomicRanges::coverage(x = gr) # Get the coverage into RLE

    gr_cov <- as(cov, "GRanges") # Coverage RLE to GRanges

    gr_cov <- gr_cov[gr_cov$score > 0] # Drop 0-score rows

    # Get the seq info needed for bigwig files
    info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    info <- info[names(info)[names(info) %in% seqnames(gr_cov)]] 
    #info <- sort(info)
    
    seqlevels(gr_cov) <- seqlevels(info)
    seqinfo(gr_cov) <- info

    # Normalise the counts (counts per million insertions)
    million_counts <- length(gr) / 1e6

    if(normalise == TRUE){
        gr_cov$score <- gr_cov$score / million_counts
    }
    
    
    # Write the bigwig file
    out <- str_replace(string = basename(fl), pattern = ".bed", replacement = "")
    out <- str_c("bigwigs/", out, "_atac_insertions.bigwig")
    
    rtracklayer::export.bw(object = gr_cov, con = out)
    
}

mclapply(bed_fls, bed_to_cov, mc.cores = 1)


lapply(bed_fls, function(x){
  
  out <- str_replace(string = basename(x), pattern = ".bed", replacement = "")
  out <- str_c("bigwigs/", out, "_atac_insertions.bigwig")
  gr_cov <- rtracklayer::import.bw(out)
  print(length(gr_cov))
  gr_cov <- trim(gr_cov)
  print(length(gr_cov))
  out1 <- gsub(".bigwig", "_trim.bigwig", out)
  
  rtracklayer::export.bw(object = gr_cov, con = out1)
  
})

all_peaks <- readRDS("processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

for(i in 1:length(all_peaks)) {
  
  gr <- all_peaks[[i]]
  gr@elementMetadata <- DataFrame(gr@elementMetadata$peakType)
  
  rtracklayer::export(gr, con=paste0("bigwigs/", names(all_peaks)[i],
    "_peaks_final_merge.bed"))
  
  
}




  
  
  
  
