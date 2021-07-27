###############################################################################
#                                                                             #
#                                                                             #
# Add identifier for stages                                                   #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(rtracklayer)
library(reshape2)
library(dplyr)
library(readxl)
library(SingleCellExperiment)

# 1. read in data

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

processed_data <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.Rds")

# 2. add stage accessibility 

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

saveRDS(processed_data, file=
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_gr.Rds")
