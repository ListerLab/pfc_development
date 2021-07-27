

# Bed file to GRanges
bed_to_gr <- function(bed){
    peaks <- fread(bed)
    gr <- GRanges(seqnames = peaks$V1,
                  ranges = IRanges(start = peaks$V2, end = peaks$V3))
    return(gr)
}
