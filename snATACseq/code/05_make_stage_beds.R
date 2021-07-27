###############################################################################
#                                                                             #
#                                                                             #
# Making normalized bigWigs                                                   #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(ArchR)
library(rtracklayer)
library(parallel)
library(readxl)
addArchRGenome("hg19")
addArchRThreads(threads = 1)

source("snATACseq/R/functions_peak_calling.R")

# 1. read in data
sc <- loadArchRProject("snATACseq/clustering_annotation")

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)

# 2. make list of cells for each trajectory and stage
ArrowFiles <- getArrowFiles(sc)
Groups <- getCellColData(ArchRProj = sc, 
    select = "Anno1", drop = TRUE)

Cells <- sc$cellNames

cellGroups <- split(Cells, Groups)
cellGroups_new <- list()
for(i in 1:length(cellGroups)){
    if(names(cellGroups)[i] %in% c("L5_6", "L4", "L2_3")){
        a <- c(i, which(names(cellGroups)=="PN_dev"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("MGE_der", "CGE_der")){
        a <- c(i, which(names(cellGroups)=="IN_dev"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("Oligo")){
        a <- c(i, which(names(cellGroups)=="OPC"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("Astro", "Micro", "Vas")){
        cellGroups_new[[names(cellGroups)[i]]] <- cellGroups[[names(cellGroups)[i]]]
    }
    
}
input <- lapply(1:length(cellGroups_new), function(x) 
    lapply(names(ArrowFiles), function(y)
    if(sum(grepl(paste0(y, "#"), cellGroups_new[[x]]))>=40) 
    c(y, names(cellGroups_new)[x]))) 
input <- unlist(input, recursive = FALSE)
input <- input[!sapply(input, is.null)]
input <- Reduce(rbind, input)
input <- as.data.frame(input)
input <- split(input$V1, input$V2)

# 3. make bed files for all stage and trajectory combinations
setwd("snATACseq/bedfiles")
system("mkdir stages")

make_bed_stages <- function(cellGroup, input, atac_samples){
    

    stages <- unique(atac_samples$Stage[match(input[[cellGroup]],
            atac_samples$Sample)])
    for(j in 1:length(stages)){
        
        all_samples <- intersect(input[[cellGroup]],
            atac_samples$Sample[atac_samples$Stage==stages[j]])
        all_samples <- paste0(cellGroup, "_",
              all_samples, "_blacklistrm.bed")
        bed_files <- list.files()
        all_samples <- intersect(all_samples, bed_files)
        if(length(all_samples)==0) {
            return(NULL)
        } else {
            all_beds <- concateBeds(all_samples)
            export(all_beds, con=paste0("stages/", 
                cellGroup, "_", stages[j], ".bed"), format="bed")

        }
    }
}


mclapply(names(input), function(x) make_bed_stages(x, input, 
    atac_samples), mc.cores=1)

