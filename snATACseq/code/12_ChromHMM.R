###############################################################################
#                                                                             #
#                                                                             #
# Overlap ChromHMM                                                            #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(rtracklayer)
library(parallel)

# 1. overlap with CREs

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

cis <- lapply(peaks, function(x) x[!is.na(x$CRE)])
cis <- as(cis, "GRangesList")
cis <- unlist(cis)
cis <- reduce(cis)

encode_anno <- read.table("annotation/encode/EIDlegend.txt", sep="\t")

ind <- grepl("Liver$|Brain|Esophagus$|Lung$|Ovary$|Placenta$|Lung$|Pancreas$|Spleen$|Thymus$|Left Ventricle$", 
             encode_anno$V2)
files <- paste0("annotation/encode/", 
    encode_anno$V1[ind], "_25_imputed12marks_dense.bed.gz")


compare_encode <- function(file, cis) {
    
    chromHMM <- import(file)
    chromHMM_enh <- chromHMM[grepl("Enh", chromHMM$name)] 
    chromHMM_het <- chromHMM[grepl("Het", chromHMM$name)] 
    chromHMM_prom<- chromHMM[grepl("Prom", chromHMM$name)]
    enh_prop <- length(unique(findOverlaps(cis, chromHMM_enh)@from))/length(cis)
    het_prop <- length(unique(findOverlaps(cis, chromHMM_het)@from))/length(cis)
    prom_prop <- length(unique(findOverlaps(cis, chromHMM_prom)@from))/length(cis)
    return(c(Enhancer=enh_prop, Heterochromatin=het_prop, Promoters=prom_prop))
    
}

all_encode_comparisons <- mclapply(files, function(x) compare_encode(x, cis),
    mc.cores=5)

all_encode_comparisons <- do.call(rbind, all_encode_comparisons)
rownames(all_encode_comparisons) <- encode_anno$V2[ind]  

write.table(all_encode_comparisons, 
  file="snATACseq/processed_data/chrom_hmm_cre.txt", quote=FALSE, sep="\t")

# 2. cell type specific overlap

all_encode_comparisons <- list()

for(i in 1:length(peaks)){
    
    tmp <- mclapply(files, function(x) 
        compare_encode(x, peaks[[i]]),
         mc.cores=5)
    
    tmp <- do.call(rbind, tmp)
    tmp <- as.data.frame(tmp)
    tmp$cell_type <- names(peaks)[i] 
    tmp$tissue<- encode_anno$V2[ind] 
    all_encode_comparisons[[i]] <- tmp
    
}

all_encode_comparisons  <- do.call(rbind, all_encode_comparisons)
write.table(all_encode_comparisons, 
  file="snATACseq/processed_data/chrom_hmm_peaks_cell_types.txt", quote=FALSE, 
  sep="\t")


# 3. overlap with all peaks

all_peaks <- unlist(as(peaks, "GRangesList"), recursive = T, use.names=T)

all_encode_comparisons <- mclapply(files, function(x) 
    compare_encode(x, all_peaks), mc.cores=5)
all_encode_comparisons <- do.call(rbind, all_encode_comparisons)
rownames(all_encode_comparisons) <- encode_anno$V2[ind]  
write.table(all_encode_comparisons, 
  file="snATACseq/processed_data/chrom_hmm_all_peaks.txt", quote=FALSE, sep="\t")




