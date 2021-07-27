###############################################################################
#                                                                             #
#                                                                             #
# Peak calling in pseudobulk                                                  #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(ArchR)
library(rtracklayer)
library(parallel)
library(readxl)
library(data.table)
library(GenomicRanges)
addArchRGenome("hg19")
addArchRThreads(threads = 1)
set.seed(10)

source("snATACseq/R/functions_peak_calling.R")

system("mkdir snATACseq/bedfiles_pseudo")
system("mkdir snATACseq/bedfiles")
system("mkdir snATACseq/peaks")

# 1. read in data and blacklist
sc <- loadArchRProject("snATACseq/clustering_annotation")
blacklist <- import("annotation/ENCFF001TDO.bed")

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

# 2. make bedfiles for pseudo-pseudobulks for cluster and sample

sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)

ArrowFiles <- getArrowFiles(sc)
Groups <- getCellColData(ArchRProj = sc, 
            select = "Anno1", drop = TRUE)

Cells <- sc$cellNames
cellGroups <- split(Cells, Groups)

availableChr <- ArchR:::.availableSeqnames(head(getArrowFiles(sc)))
chromLengths <- getChromLengths(sc)
chromSizes <- getChromSizes(sc)

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

cell_types <- names(cellGroups_new)

input <- lapply(1:length(cellGroups_new), function(x) lapply(names(ArrowFiles), 
     function(y) {
        if(sum(grepl(paste0(y, "#"), cellGroups_new[[x]]))>=40) 
        c(names(cellGroups_new)[x], y)
      }
)) 
input <- unlist(input, recursive = FALSE)
input <- input[!sapply(input, is.null)]

# get fragments from ArrowFiles for every sample and every celltype
# split into two pseudo-pseudobulks of equal size

make_beds_pseudo <- function(ArrowFile, ArrowFileName, 
                             cellGroups_new=cellGroups_new){
    for(i in 1:length(cellGroups_new)){
        
        cellGroupi <- cellGroups_new[[i]]
        cellGroupi <- cellGroupi[grepl(paste0(ArrowFileName, "#"), 
                        cellGroupi)]
        if(length(cellGroupi)<40) next
        
        ind_groupi <- sample(1:length(cellGroupi), floor(length(cellGroupi)/2))
        ind_groupi <- list(pseudo_1=ind_groupi, 
            pseudo_2=setdiff(1:length(cellGroupi), ind_groupi))
        
        for(k in 1:2){
            
            cellNames <- cellGroupi[ind_groupi[[k]]]
        
            covList <- getFrags(ArrowFile, cellNames)
            
            # remove blacklist and set strand to positive
            covList <- subsetByOverlaps(covList, blacklist, invert=TRUE)
            strand(covList) <- "+"
        
            rtracklayer::export(covList, con=paste0("snATACseq/bedfiles_pseudo/",
                names(cellGroups_new)[i],"_", ArrowFileName, "_pseudo", k,
                "_blacklistrm.bed"), format="bed")
        }
    }
}

mcmapply(function(X,Y) make_beds_pseudo(X,Y, cellGroups_new=cellGroups_new), 
    X=ArrowFiles, Y=names(ArrowFiles),
    mc.cores=length(ArrowFiles))

# 3. call peaks in pseudo-pseudobulk

call_pseudo <- function(cellGroup, ArrowFile, rep=2){
    
    for(k in 1:rep){
        macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cellGroup, "_", ArrowFile, 
        "_", "pseudo", k, "_blacklistrm.bed"), 
        paste0("snATACseq/bedfiles_pseudo/", cellGroup, "_", ArrowFile, "_pseudo", k))
    }
    
}

mclapply(input, function(x) call_pseudo(x[1], x[2], rep=2), mc.cores=10)
system("rm snATACseq/bedfiles_pseudo/*peaks.xls")

# 4. make bedfiles and call peaks for pseudo-pseudobulks for cell types

cell_types <- names(cellGroups_new)
bed_files <- list.files("snATACseq/bedfiles_pseudo")

call_pseudo_cluster <- function(cell_type, bed_files){

    pseudo_files <- bed_files[grepl(cell_type, bed_files)]
    pseudo_files <- pseudo_files[grepl(cell_type, pseudo_files)]
    pseudo_files <- pseudo_files[grepl("_blacklistrm.bed", pseudo_files)]

    pseudo_files1 <- paste0("snATACseq/bedfiles_pseudo/", 
        pseudo_files[grepl("_pseudo1", pseudo_files)])
    pseudo_files2 <- paste0("snATACseq/bedfiles_pseudo/",
        pseudo_files[grepl("_pseudo2", pseudo_files)])
    
    all_pseudo_files_1 <- concateBeds(pseudo_files1)
    export(all_pseudo_files_1, con=paste0("snATACseq/bedfiles_pseudo/", 
        cell_type, "_pseudo1.bed"), format="bed")
    all_pseudo_files_2 <- concateBeds(pseudo_files2)
    export(all_pseudo_files_2, con=paste0("snATACseq/bedfiles_pseudo/",
        cell_type, "_pseudo2.bed"), format="bed")
    
    macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo1.bed"), 
        paste0("bedfiles_pseudo/", cell_type, "_pseudo1"))
    macsPeaks(paste0("snATACseq/bedfiles_pseudo/", cell_type, "_pseudo2.bed"), 
        paste0("bedfiles_pseudo/", cell_type, "_pseudo2"))
    
}

mclapply(cell_types, function(x) call_pseudo_cluster(x, bed_files), mc.cores=10)
 
# 5. make bedfiles for pseudobulk for cluster and sample

make_beds<- function(ArrowFile, ArrowFileName, 
                             cellGroups_new=cellGroups_new){
    for(i in 1:length(cellGroups_new)){
        cellGroupi <- cellGroups_new[[i]]
        cellGroupi <- cellGroupi[grepl(paste0(ArrowFileName, "#"), 
                                   cellGroupi)]
        if(length(cellGroupi)<40) next
        
        covList <- getFrags(ArrowFile, cellGroupi)
        covList <- subsetByOverlaps(covList, blacklist, invert=TRUE)
        
        strand(covList) <- "+"
        
        rtracklayer::export(covList, con=paste0("snATACseq/bedfiles/",
            names(cellGroups_new)[i], "_", ArrowFileName, 
            "_blacklistrm.bed"), format="bed")
    }
}

mcmapply(function(X,Y) make_beds(X, Y, cellGroups_new), 
         X=ArrowFiles, Y=names(ArrowFiles),
         mc.cores=length(ArrowFiles))

# 6. call peaks for pseudobulk

call_peaks <- function(cellGroup, ArrowFile){
    
    macsPeaks(paste0("snATACseq/bedfiles/", cellGroup, "_", ArrowFile, "_blacklistrm.bed"), 
            paste0("snATACseq/bedfiles/", cellGroup, "_", ArrowFile))

}


mclapply(input, function(x) call_peaks(x[1], x[2]), mc.cores=10)

# 7. make bedfiles and call peaks for pseudobulks for clusters

bed_files <- list.files("snATACseq/bedfiles")
    
call_cluster <- function(cell_type, bed_files){
        
    pseudo_files <- bed_files[grepl(cell_type, bed_files)]
    pseudo_files <- pseudo_files[grepl(cell_type, pseudo_files)]
    pseudo_files <- pseudo_files[grepl("_blacklistrm.bed", pseudo_files)]
        
    all_pseudo_files <- concateBeds(paste0("snATACseq/bedfiles/", pseudo_files))
    export(all_pseudo_files, con=paste0("snATACseq/bedfiles/", 
        cell_type, "_blacklistrm.bed"), format="bed")
    
    macsPeaks(paste0("snATACseq/bedfiles/", cell_type, "_blacklistrm.bed"), 
        paste0("snATACseq/bedfiles/", cell_type))
    
}

mclapply(cell_types, function(x) call_cluster(x, bed_files), mc.cores=10)

# 8. call IDR and merge

link1 <- "snATACseq/bedfiles_pseudo/"
link2 <- "snATACseq/bedfiles/"
link3 <- "snATACseq/peaks/"

idr_call <- function(input1, link1=link1, link2=link2, link3=link3){
    
    system(paste0("source activate idr_filtering
                  idr --samples ", link1, input1, "_pseudo1_peaks.narrowPeak ",  
                  link1, input1,
                  "_pseudo2_peaks.narrowPeak --peak-list ", link2, 
                  input1,
                  "_peaks.narrowPeak --input-file-type ",
                  "narrowPeak --output-file ", link3, input1,
                  "_IDR0.05.narrowPeak ",
                  "--rank p.value --soft-idr-threshold 0.05 ",
                  "--use-best-multisummit-IDR"))
    idr_trans <- -log(0.05)/log(10)
    dat <- fread(paste0(link3, input1, "_IDR0.05.narrowPeak"), data.table = FALSE)
    ind <- (dat[,9] >= idr_trans & dat[,12] >= idr_trans)
    dat <- dat[ind,]
    g_dat <- GRanges(seqnames=dat[,1], IRanges(start=dat[,2], end=dat[,3]),
            score=dat[,7])
    if(length(g_dat)>0){
        g_dat <- ArchR::nonOverlappingGR(g_dat, by="score", decreasing=TRUE)
        export(g_dat, con=paste0(link3, 
                input1, "_peaks_final_merge.bed"), format="bed")
    }
    
}


mclapply(input, function(x) idr_call(paste0(x[1], "_", x[2]), link1=link1,
    link2=link2, link3=link3), mc.cores=10)
mclapply(cell_types, function(x) idr_call(x, link1=link1,
    link2=link2, link3=link3), mc.cores=10)

