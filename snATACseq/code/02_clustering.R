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

addArchRGenome("hg19")
addArchRThreads(threads = 20)

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
# subsample only to atac_samples used
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

# 1. make arrow files
link <- "path_to_scatac_raw"

for(i in 1:length(atac_samples$Sample)){
    
    # use only valid 10x barcodes 
    barcodes <- getValidBarcodes(csvFiles = paste0(link, atac_samples$Sample[i], 
        "/outs/singlecell.csv"), sampleNames = 
        paste0(link, atac_samples$Sample[i], 
        "/outs/filtered_peak_bc_matrix/barcodes.tsv"))
    
    ArrowFiles <- createArrowFiles(inputFiles = paste0(link, 
        atac_samples$Sample[i], "/outs/fragments.tsv.gz"),
        sampleNames = atac_samples$Sample[i], filterTSS = 0, filterFrags = 0, 
        validBarcodes = barcodes[[1]], addTileMat = TRUE, 
        addGeneScoreMat = TRUE, force=TRUE, offsetPlus = 0, offsetMinus = 0)
}

# 2. find doublets

ArrowFiles <- paste0("arrows_good/", atac_samples$Sample, ".arrow")
system("mv *.arrow arrows_good/")

set.seed(1)
doubScores <- addDoubletScores(input = ArrowFiles, k = 20, 
    knnMethod = "UMAP", LSIMethod = 1, force=TRUE)

# 3. make ArchR project

sc <- ArchRProject(ArrowFiles = paste0(ArrowFiles), 
    outputDirectory = "snATACseq/clustering", copyArrows = FALSE)

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

# 5. filter doublets

sc <- filterDoublets(sc)

# 6. dimension reduction 

sc <- addIterativeLSI(ArchRProj = sc, useMatrix = "TileMatrix", 
    name = "IterativeLSI", iterations = 4, clusterParams = list( 
    resolution = 4, sampleCells = 10000, n.start = 10), 
    varFeatures = 50000, dimsToUse = 1:20, force=TRUE, seed=11)

# 7. add batch correction

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

# 10. plots

p2 <- plotEmbedding(ArchRProj = sc, colorBy = "cellColData", name = "Clusters", 
    embedding = "UMAP")
plotPDF(p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = sc, 
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

# 13. add meta information

sc$Stage <- atac_samples$Stage[match(sc$Sample, atac_samples$Sample)]
sc$Age <- atac_samples$Age[match(sc$Sample, atac_samples$Sample)]
sc$Sex <- atac_samples$Sex[match(sc$Sample, atac_samples$Sample)]
sc$arcsin_ages <-as.numeric(
    atac_samples$`ArcSin Age`[match(sc$Sample, atac_samples$Sample)])

saveArchRProject(ArchRProj = sc, load = FALSE)

# 14. save UMAP and colData 

umap <- sc@embeddings@listData$UMAP$df
coldata <- sc@colCellData

saveRDS(list(coldata=coldata, df_umap=umap), 
        file="snATACseq/processed_data/UMAP_all.rds")
