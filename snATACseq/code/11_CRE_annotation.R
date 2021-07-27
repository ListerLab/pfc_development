###############################################################################
#                                                                             #
#                                                                             #
# Add correlated gene to annotation                                           #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(rtracklayer)
library(parallel)

# 1. read in data

processed_data <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_gr.Rds")
pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")


#2. annotate with CREs

for(i in 1:length(processed_data)){
    
    peak_names <- paste0(names(processed_data)[i], ".", 
        names(processed_data[[i]]))
    ind <- mclapply(peak_names, function(x)
        which(pearson_corr_sig$peak_name %in% x), mc.cores=10)
    all_genes <- sapply(ind, function(x) {
        if(length(x)==0) return(NA)
        paste0(pearson_corr_sig$Gene[x], collapse="|")
    })
    processed_data[[i]]$CRE <- unlist(all_genes)
    
    
}

# 3. save peaks

saveRDS(processed_data, file=
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
