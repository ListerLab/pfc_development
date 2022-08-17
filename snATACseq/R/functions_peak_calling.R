## Functions for peak calling

getFrags <- function(ArrowFile, cellNames){ 

    covList <- lapply(availableChr, function(kk){
    
        ArchR:::.getFragsFromArrow(ArrowFile, 
             chr = kk, out = "GRanges", 
             cellNames = cellNames)
    })
    covList <- do.call(c, covList)

    covList <- c(GRanges(seqnames = seqnames(covList),
        ranges = IRanges(start = start(covList), end = start(covList))),
        GRanges(seqnames = seqnames(covList),
        ranges = IRanges(start = end(covList), end = end(covList))))

    covList <- covList + 25 
    covList <- sort(covList)
    
}

concateBeds <- function(x){
    
    all_beds <- lapply(x, function(y) rtracklayer::import(y))
    all_beds <- as(all_beds, "GRangesList")
    all_beds <- unlist(all_beds)
    return(all_beds)
    
}

macsPeaks <- function(input, output){
    
    system(paste0("macs2 callpeak -t ", input,
       " -f BED -n ", output, 
       " -g hs -p 0.01 ",
       "--nomodel --keep-dup all --call-summits --nolambda"))
}