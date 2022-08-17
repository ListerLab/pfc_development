################################################################################
#                                                                              #
#                                                                              #
# Footprinting plots                                                           #
#                                                                              #
#                                                                              #    
################################################################################

library(ArchR)
library(SummarizedExperiment)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

cell_types <- list(Astro=c("Astro"), L5_6=c("PN dev", "L5/6"),
    L4=c("PN dev", "L4"), L2_3=c("PN dev", "L2/3"), 
    MGE_der=c("IN dev", "MGE der"), CGE_der=c("IN dev", "CGE der"),
    Mirco=c("Micro"), Oligo=c("OPC", "Oligo"))

gr_trends <-readRDS("snATACseq/processed_data/gene_trends_gr.rds")
tf_trends <-read.table("snATACseq/processed_data/gene_trends_tf.tsv", sep="\t",
      stringsAsFactors = FALSE)
rownames(tf_trends) <- NULL
tf_trends$id <- paste0(tf_trends$cell_type, "_", tf_trends$feature)
tf_trends <- tf_trends[!duplicated(tf_trends$id),]
tf_trends <- tf_trends[tf_trends$padjust<0.05,]
tab <- table(tf_trends$feature)
tf_all <- names(tab)[tab>7]

all_features <- lapply(tf_all, function(x) 
  tf_trends$cell_type[tf_trends$feature==x])
names(all_features) <- tf_all

all_features <- lapply(all_features, function(x) "G7")


peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

# 2. annotate with peaks and motifs, make peak matrix

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
         "nearestGene", "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

peaks <- unlist(as(peaks, "GRangesList"))

sc <- addPeakSet(sc, peaks, force=TRUE)
sc <- addPeakMatrix(sc, force=TRUE)
sc <- addMotifAnnotations(ArchRProj = sc, motifSet = "JASPAR2018", 
                          name = "Motif", force=TRUE)


motifs <- split(tf_trends$feature, tf_trends$cell_type)
motifs <- motifs[!lengths(motifs)==0]

sc <- addGroupCoverages(ArchRProj = sc, groupBy = "Stage")
motifPositions <- getPositions(sc)

xxx <- list(SOX10_405=c("G12"), TFEC_340=c("G1"))
xxx <- all_features

xxx <- list(POU2F1_256=c("G11"), TEF_318=c("G3"),
      SOX10_405=c("G2"), SOX10_405=c("G8"), SOX15_445=c("G8"),
      SOX15_445=c("G2"), STAT3_68=c("G7"))

xxx <- list(THAP1_75=c("G7"), SREBF1_73=c("G7"),
    SREBF2_74=c("G7"), EBF1_90=c("G7"), HIC2_204=c("G7"))

for(i in 5:length(xxx)){

    mm <- motifPositions[names(xxx)[i]]
    gene_trends_tmp <- gr_trends[names(gr_trends) %in% xxx[[i]]]
    gene_trends_tmp <- lapply(gene_trends_tmp, function(x)
        reduce(x))
    gene_trends_tmp <- unlist(as(gene_trends_tmp, "GRangesList"))
    mm_interest <- lapply(mm, function(x) x[findOverlaps(x, 
        gene_trends_tmp)@from])
    mm_interest <- lapply(mm_interest, function(x) 
      x[!duplicated(x)])
    
    seFoot <- getFootprints(
        ArchRProj = sc, 
        positions = mm_interest,
        groupBy = "Stage"
    )
    
    cn_tmp <- tf_trends$cluster_name[match(names(xxx)[i], tf_trends$feature)[1]]
    cn_tmp <- gsub("/", "", cn_tmp, fixed=T)
        
    
    ggsave(plot_footprint(seFoot, names(xxx)[i]), 
           file=paste0("snATACseq/region_figures/footprints/cc_", names(xxx)[i], 
            "_", xxx[i] , "_", cn_tmp, ".pdf"), height=10, width=10)
}

plot_footprint <- function(seFoot, 
                           #seFoot_ctrl, 
                           motif){

    rowDF <- SummarizedExperiment::rowData(seFoot)
    footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
    biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

    footMat <- ArchR:::.getAssay(seFoot[
        BiocGenerics::which(rowDF[,2]=="footprint"),], motif)
    biasMat <- ArchR:::.getAssay(seFoot[
        BiocGenerics::which(rowDF[,2]=="bias"),], motif)

    footMat <- footMat - biasMat

    footMatMean <- apply(footMat, 1, function(x) mean(x)) 
    footMatSd <- apply(footMat, 1, function(x) sd(x))

    plotFootDF <- data.frame(
        x = footDF$x, 
        mean = footMatMean, 
        sd = footMatSd
    )

    plotFootDF$mean1 <- c(NA, NA, NA, sapply(1:(length(plotFootDF$mean)-5), 
        function(x) mean(plotFootDF$mean[x:(x+5)])), NA, NA) 

    gg <- ggplot(plotFootDF, aes(x = x, y = mean1)) + 
    geom_rect(mapping=aes(xmin=-10, xmax=10,
     ymin=0, ymax=max(plotFootDF$mean1, na.rm=T)+10), 
     alpha=0.2, fill="#ffffba") +
    geom_line() + 
    xlab("Distance to motif center (bp)") +
    ylab(paste0("Normalized Insertions")) +
    theme_classic() + ggtitle(motif)
}



