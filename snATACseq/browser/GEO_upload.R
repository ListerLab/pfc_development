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


peaks <- as(peaks, "GRangesList")
peaks <- unlist(peaks, recursive=TRUE)
peaks@elementMetadata <- DataFrame(OCR=rep("OCR", length(peaks)),
    id=names(peaks))
write.table(peaks, file=paste0("snATACseq/GEO_upload/all_peaks.tsv"), 
            quote=FALSE, sep="\t", row.names = T,
            col.names = T)

paths <- "snRNAseq/processed_data"
sce <- readRDS(paste0(paths, 
                      "/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"))
write.table(rowData(sce), file="snATACseq/GEO_upload/gene_info.txt",
            row.names=TRUE, quote=FALSE, col.names=TRUE)
