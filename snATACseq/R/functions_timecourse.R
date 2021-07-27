annotate_peaks <- function(peak_set, txdb){
    
    gens <- genes(txdb)
    nearest_gene <- names(gens)[nearest(peak_set, gens, 
            select=c("arbitrary"), ignore.strand=TRUE)]
    dist_gene <- distanceToNearest(peak_set, gens, 
            ignore.strand=TRUE)
    ind <- dist_gene@elementMetadata$distance>=250000
    nearest_gene[ind] <- ""
    peak_set@elementMetadata$nearestGene <- nearest_gene
    
    peak_set@elementMetadata$peakType <- "distal"
    introns <- intronsByTranscript(txdb)
    introns <- unlist(introns)
    strand(introns) <- "*"
    ind <- findOverlaps(peak_set, introns)
    peak_set@elementMetadata$peakType[unique(ind@from)] <- "introns"
    exons <- exons(txdb)
    strand(exons) <- "*"
    ind <- findOverlaps(peak_set, exons)
    peak_set@elementMetadata$peakType[unique(ind@from)] <- "exons"
    promoters <- promoters(txdb, upstream = 500, downstream = 100)
    strand(promoters) <- "*"
    ind <- findOverlaps(peak_set, promoters)
    peak_set@elementMetadata$peakType[unique(ind@from)] <- "promoters"
    
    return(peak_set)
    
    
}