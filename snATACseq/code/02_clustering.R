###############################################################################
#                                                                             #
#                                                                             #
# Clustering of snATACseq                                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(ArchR)
library(ggplot.multistats)
library(readxl)
library(scuttle)


addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
# subsample only to atac_samples used
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

# 1. make arrow files
link <- "" #insert path to raw data

system("mkdir snATACseq/arrows_good")
setwd("snATACseq/arrows_good")

for(i in 1:length(atac_samples$Sample)){
    
    # use only valid 10x barcodes 
    barcodes <- getValidBarcodes(csvFiles = paste0(link, atac_samples$Sample[i], 
        "/outs/singlecell.csv"), sampleNames = 
        paste0(link, atac_samples$Sample[i], 
        "/outs/filtered_peak_bc_matrix/barcodes.tsv"))
    
    ArrowFiles <- createArrowFiles(inputFiles = paste0(link, 
        atac_samples$Sample[i], "/outs/fragments.tsv.gz"),
        sampleNames = atac_samples$Sample[i], minTSS = 0, minFrags = 0, 
        validBarcodes = barcodes[[1]], addTileMat = TRUE, 
        addGeneScoreMat = TRUE, force=TRUE, offsetPlus = 0, offsetMinus = 0,
        excludeChr = c("chrM", "chrY", "chrX"))
}

# 2. find doublets

setwd("../..")

ArrowFiles <- paste0("snATACseq/arrows_good/", atac_samples$Sample, ".arrow")

set.seed(1)
doubScores <- addDoubletScores(input = ArrowFiles, k = 20, 
    knnMethod = "UMAP", LSIMethod = 1, force=TRUE)

# 3. make ArchR project

sc <- ArchRProject(ArrowFiles = ArrowFiles, 
    outputDirectory = "snATACseq/clustering", copyArrows = TRUE)

sc <- loadArchRProject("snATACseq/clustering")

# 4. quality control plots

p1 <- plotGroups(ArchRProj = sc, groupBy = "Sample", 
    colorBy = "cellColData", name = "TSSEnrichment",plotAs = "ridges")
p2 <- plotGroups(ArchRProj = sc, groupBy = "Sample", 
    colorBy = "cellColData", name = "TSSEnrichment",
    plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
p3 <- plotGroups( ArchRProj = sc, groupBy = "Sample", 
    colorBy = "cellColData", name = "log10(nFrags)",
    plotAs = "ridges")
p4 <- plotGroups(ArchRProj = sc, groupBy = "Sample", 
    colorBy = "cellColData", name = "log10(nFrags)",
    plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
plotPDF(p1, p2, p3, p4, name = "QC-Sample-Statistics.pdf", 
        ArchRProj = sc, addDOC = FALSE, width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = sc)
p2 <- plotTSSEnrichment(ArchRProj = sc)
plotPDF(p1, p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = sc, addDOC = FALSE, width = 5, height = 5)

# 5. add meta information

sc$Stage <- atac_samples$Stage[match(sc$Sample, atac_samples$Sample)]
sc$Age <- atac_samples$Age[match(sc$Sample, atac_samples$Sample)]
sc$Sex <- atac_samples$Sex[match(sc$Sample, atac_samples$Sample)]
sc$arcsin_ages <-as.numeric(
    atac_samples$`ArcSin Age`[match(sc$Sample, atac_samples$Sample)])

# 6. filter doublets

sc <- filterDoublets(sc)

#Filtering 8162 cells from ArchRProject!
#RL1784 : 511 of 7154 (7.1%)
#RL1785 : 491 of 7011 (7%)
#RL1914 : 212 of 4606 (4.6%)
#RL1994 : 19 of 1383 (1.4%)
#RL2085 : 412 of 6425 (6.4%)
#RL2207 : 421 of 6490 (6.5%)
#RL2208 : 348 of 5904 (5.9%)
#RL2209 : 565 of 7517 (7.5%)
#RL2210 : 542 of 7364 (7.4%)
#RL2366 : 46 of 2162 (2.1%)
#RL2367 : 449 of 6708 (6.7%)
#RL2368 : 466 of 6833 (6.8%)
#RL2369 : 2217 of 14890 (14.9%)
#RL2371 : 226 of 4759 (4.7%)
#RL2372 : 67 of 2601 (2.6%)
#RL2373 : 1020 of 10103 (10.1%)
#RL2364 : 150 of 3883 (3.9%)

#7. filter TSS 

tss_outliers <- list()
tss_outliers_names <- list()

for(i in 1:length(atac_samples$Sample)){
    
    sample_i <- atac_samples$Sample[i]
    tss_enrich <- sc$TSSEnrichment[sc$Sample==sample_i]
    tss_outliers[[sample_i]] <- isOutlier(tss_enrich, nmads=1,  type="lower")
    tss_outliers_names[[sample_i]] <- 
        rownames(sc)[sc$Sample==sample_i][tss_outliers[[sample_i]]]
}

lower_threshold <- sapply(tss_outliers, function(x) attr(x, "thresholds")[1])
#RL1784.lower RL1785.lower RL1914.lower RL1994.lower RL2085.lower RL2207.lower 
#2.265997     5.124557     3.832662     2.522640     4.082562     4.126503 
#RL2208.lower RL2209.lower RL2210.lower RL2366.lower RL2367.lower RL2368.lower 
#2.051703     2.774853     1.599835     4.536740     7.234283     3.720166 
#RL2369.lower RL2371.lower RL2372.lower RL2373.lower RL2364.lower 
#2.327747     3.236841     1.113877     1.751072     1.806042 
filter_tss <- sapply(tss_outliers, function(x) sum(x))
#RL1784 RL1785 RL1914 RL1994 RL2085 RL2207 RL2208 RL2209 RL2210 RL2366 RL2367 
#34    276    420      0    283    761    190    355    121    252    909 
#RL2368 RL2369 RL2371 RL2372 RL2373 RL2364 
#426    211    553    290    565    365 

saveRDS(sc@cellColData, 
    file="snATACseq/processed_data/scatac_cell_data_prior_filtering.rds")

tss_outliers_names <- unlist(tss_outliers_names)
sc <- sc[!rownames(sc) %in% tss_outliers_names,]

# 8. dimension reduction 

sc <- addIterativeLSI(ArchRProj = sc, useMatrix = "TileMatrix", 
    name = "IterativeLSI", iterations = 4, clusterParams = list( 
    resolution = 4, sampleCells = 10000, n.start = 10), 
    varFeatures = 50000, dimsToUse = 1:20, force=TRUE, seed=11)

# 9. add batch correction

sc <- addHarmony(ArchRProj = sc, reducedDims = "IterativeLSI",
    name = "Harmony", groupBy = "Sample", force=TRUE)

# 8. clustering

sc <- addClusters(input = sc, reducedDims = "Harmony", method = "Seurat",
    name = "Clusters", resolution = 4, force=TRUE, seed = 12,
    maxClusters = 100)

# 9. UMAP

sc <- addUMAP(ArchRProj = sc, reducedDims = "Harmony", name = "UMAP", 
    nNeighbors = 3, minDist = 0.08, metric = "cosine",
    force=TRUE, seed = 13)

saveArchRProject(sc, outputDirectory="snATACseq/clustering")

# 10. plots

p2 <- plotEmbedding(ArchRProj = sc, colorBy = "cellColData", name = "Clusters", 
    embedding = "UMAP")
plotPDF(p2, name = "snATACseq/clustering/Plots/Plot-UMAP-Sample-Clusters.pdf", ArchRProj = sc, 
    addDOC = FALSE, width = 5, height = 5)

# make confusion matrix
cM <- confusionMatrix(paste0(sc$Clusters), paste0(sc$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), border_color = "black")
ggsave(p, file=paste0("snATACseq/clustering/Plots/Heatmap_clustering.pdf"))
write.csv(as.matrix(cM), file="snATACseq/clustering/Plots/Cluster_sample.csv", quote=F)

# make special UMAPs for sample location, etc.
df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample[match(rownames(df), sc$cellNames)]
df$Cluster <- sc$Clusters[match(rownames(df), sc$cellNames)]
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Sample)) + geom_point(size=0.2) +
    theme_classic()
g2 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Cluster)) + geom_point(size=0.2) + 
    theme_classic()

system("mkdir snATACseq/clustering/Plots/samples")

for(i in 1:length(atac_samples$Sample)){
    
    df$SampleInterest <- (df$Sample==atac_samples$Sample[i])*1
    
    gg <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
        stat_summaries_hex( aes(z = SampleInterest, fill = stat(mean)),
        funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
        theme_minimal() +
        ggtitle(paste0(atac_samples$Sample[i], "-", atac_samples$Age[i]))
    
    ggsave(gg, file=paste0("snATACseq/clustering/Plots/",
         "samples/UMAP_", atac_samples$Sample[i], ".pdf"))
    
}

ggsave(g1, file="snATACseq/clustering/Plots/UMAP_Sample.pdf")
ggsave(g2, file="snATACseq/clustering/Plots/UMAP_Cluster.pdf")


