# save for GEO upload

library(rtracklayer)
library(scater)
library(DropletUtils)

peaks <- readRDS("snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
for(i in 1:length(peaks)){
    
    tmp <- as.data.frame(peaks[[i]])
    write.table(tmp, file=paste0("snATACseq/GEO_upload/",
        names(peaks)[i], ".tsv"), quote=FALSE, sep="\t", row.names = F,
        col.names = T)
    
}

sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

