library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(gridExtra)
library(viridis)
library(org.Hs.eg.db)
library(cba)
library(clusterProfiler)

pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")

peaks_per_gene <- split(pearson_corr_sig$peak_name, pearson_corr_sig$Gene)

find_number_unique_peaks <- function(x, pearson_corr_sig) {
    
    ind <- match(x, pearson_corr_sig$peak_name)
    gr <- GRanges(seqnames = pearson_corr_sig$seqnames[ind],
                  IRanges(start=pearson_corr_sig$start[ind], 
                          end=pearson_corr_sig$end[ind]))
    gr <- reduce(gr)
    return(length(gr))
    
}

no_peaks_genes <- mclapply(peaks_per_gene, function(x) find_number_unique_peaks(x, 
     pearson_corr_sig), mc.cores=3)

median(unlist(no_peaks_genes))

hist(unlist(no_peaks_genes), breaks=30)
dat <- data.frame(Freq=table(unlist(no_peaks_genes)))
colnames(dat) <- c("No CREs", "Freq")
dat$`No CREs` <- as.numeric(dat$`No CREs`)
gg <- ggplot(dat, aes(x=`No CREs`, y=Freq)) + geom_col() + theme_classic()

ggsave(gg, file="supp_figures/SuppFig5_Hist_CRE.svg",
       height=7, width=7)

