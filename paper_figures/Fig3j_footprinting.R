## Footprinting

library(ArchR)
library(SummarizedExperiment)
library(ggplot.multistats)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_annotation")
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

cell_types <- list(Astro=c("Astro"), L5_6=c("PN dev", "L5/6"),
    L4=c("PN dev", "L4"), L2_3=c("PN dev", "L2/3"), 
    MGE_der=c("IN dev", "MGE der"), CGE_der=c("IN dev", "CGE der"),
    Mirco=c("Micro"), Oligo=c("OPC", "Oligo"))

gr_trends <-readRDS("snATACseq/processed_data/tf_gene_trends_gr.RDS")
tf_trends <-readRDS("snATACseq/processed_data/tf_gene_trends_sub.RDS")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")

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

xxx <- list(SOX10_405=c(0,2,4,10), POU3F3_259=c(8),
    HSF4_239=c(3,4,7))

for(i in 1:length(xxx)){

    mm <- motifPositions[names(xxx)[i]]
    gene_trends_tmp <- gr_trends[names(gr_trends) %in% xxx[[i]]]
    gene_trends_tmp <- lapply(gene_trends_tmp, function(x)
        reduce(x))
    gene_trends_tmp <- unlist(as(gene_trends_tmp, "GRangesList"))
    mm_interest <- lapply(mm, function(x) x[findOverlaps(x, 
        gene_trends_tmp)@from])
    #mm_control <- lapply(mm, function(x) x[-findOverlaps(x, 
    #    gene_trends_tmp)@from])
    #set.seed(10)
    #mm_control <- mapply(function(X,Y) X[sample(1:length(X), Y)],
    #    X=mm_control, Y=lengths(mm_interest))
    
    seFoot <- getFootprints(
        ArchRProj = sc, 
        positions = mm_interest,
        groupBy = "Stage"
    )
    
    seFoot_ctrl <- getFootprints(
        ArchRProj = sc, 
        positions = mm_control,
        groupBy = "Stage"
    )
    
    cn_tmp <- tf_trends$cluster_name[match(names(xxx)[i], tf_trends$feature)[1]]
    cn_tmp <- gsub("/", "", cn_tmp, fixed=T)
        
    #ggsave(plot_footprint(seFoot, seFoot_ctrl, names(xxx)[i]), 
    #    file=paste0("Fig3j_Footprints_", names(xxx)[i], 
    #    "_", cn_tmp, ".pdf"), height=10, width=10)
    
    ggsave(plot_footprint(seFoot, names(xxx)[i]), 
           file=paste0("Fig3j_Footprints_", names(xxx)[i], 
                       "_", cn_tmp, ".pdf"), height=10, width=10)
    
    
}

plot_footprint <- function(seFoot, 
                           #seFoot_ctrl, 
                           motif){

    rowDF <- SummarizedExperiment::rowData(seFoot)
    footDF <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
    biasDF <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

    #rowDF_ctrl <- SummarizedExperiment::rowData(seFoot_ctrl)
    #footDF_ctrl <- rowDF[BiocGenerics::which(rowDF[,2]=="footprint"),]
    #biasDF_ctrl <- rowDF[BiocGenerics::which(rowDF[,2]=="bias"),]

    footMat <- ArchR:::.getAssay(seFoot[
        BiocGenerics::which(rowDF[,2]=="footprint"),], motif)
    biasMat <- ArchR:::.getAssay(seFoot[
        BiocGenerics::which(rowDF[,2]=="bias"),], motif)
    #footMat_ctrl <- ArchR:::.getAssay(seFoot_ctrl[
    #    BiocGenerics::which(rowDF_ctrl[,2]=="footprint"),], motif)
    #biasMat_ctrl <- ArchR:::.getAssay(seFoot_ctrl[
    #    BiocGenerics::which(rowDF_ctrl[,2]=="bias"),], motif)

    footMat <- footMat - biasMat
    #footMat_ctrl <- footMat_ctrl - biasMat_ctrl

    footMatMean <- apply(footMat, 1, function(x) mean(x)) 
    footMatSd <- apply(footMat, 1, function(x) sd(x))
    #footMatMean_ctrl <- apply(footMat_ctrl, 1, function(x) mean(x)) 
    #footMatSd_ctrl <- apply(footMat_ctrl, 1, function(x) sd(x))

    plotFootDF <- data.frame(
        x = footDF$x, 
        mean = footMatMean, 
        sd = footMatSd
    )
    #plotFootDF_ctrl <- data.frame(
    #    x = footDF_ctrl$x, 
    #    mean = footMatMean_ctrl,
    #    sd = footMatSd_ctrl
    #)

    plotFootDF$mean1 <- c(NA, NA, NA, sapply(1:(length(plotFootDF$mean)-5), 
        function(x) mean(plotFootDF$mean[x:(x+5)])), NA, NA) 
    #plotFootDF_ctrl$mean1 <- c(NA, NA, NA, sapply(1:(length(plotFootDF_ctrl$mean)-5), 
    #       function(x) mean(plotFootDF_ctrl$mean[x:(x+5)])), NA, NA)

    gg <- ggplot(plotFootDF, aes(x = x, y = mean1)) + 
    geom_rect(mapping=aes(xmin=-10, xmax=10,
     ymin=0, ymax=max(plotFootDF$mean1, na.rm=T)+10), 
     alpha=0.2, fill="#ffffba") +
    geom_line() + 
    #geom_line(data=plotFootDF_ctrl,
    #    aes(x = x, y = mean1), col="grey") +
    xlab("Distance to motif center (bp)") +
    ylab(paste0("Normalized Insertions")) +
    theme_classic() + ggtitle(motif)
}



