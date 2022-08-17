###############################################################################
#                                                                             #
#                                                                             #
# Annotation of snATACseq                                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(ArchR)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(ggplot.multistats)
library(readxl)

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

# 1. load data
sc <- loadArchRProject("snATACseq/clustering")
serna <- readRDS("snRNAseq/processed_data/scrna_sum_exp.rds")

# 2. integrate snRNA-seq data

set.seed(10)
sc <- addGeneIntegrationMatrix(ArchRProj = sc, useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix", reducedDims = "Harmony",
    seRNA = serna, addToArrow = TRUE, groupRNA = "sub_clust",
    nameCell = "predictedCell", nameGroup = "predictedGroup",
    nameScore = "predictedScore", force=TRUE)

# 3. plots

pal <- paletteDiscrete(values = colData(serna)$sub_clust)

p1 <- plotEmbedding(sc, colorBy = "cellColData", name = "predictedGroup")
plotPDF(p1, name = "Plot-UMAP-RNA-Unconstrained-Integration-Sub.pdf", 
        ArchRProj = sc, addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(sc$Clusters, sc$predictedGroup)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)

p <- pheatmap::pheatmap(mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), border_color = "black")

ggsave(p, file=paste0("snATACseq/clustering/Plots/Heatmap_annotation_sub.pdf"))
write.csv(as.matrix(cM), file="snATACseq/clustering/Plots/Cluster_anno.csv", 
    quote=F)

sc$Clusters_sub <- mapLabs[match(sc$Clusters, mapLabs[,1]),2]
table(sc$Clusters_sub, sc$Sample)

p1 <- plotEmbedding(sc, colorBy = "cellColData", name = "Clusters_sub",
    embedding = "UMAP")
plotPDF(p1, name = "Plot-UMAP-RNA-Clusters_sub.pdf", 
        ArchRProj = sc, addDOC = FALSE, width = 5, height = 5)

system("mkdir snATACseq/clustering/Plots/annotation")

cell_types <- unique(sc$predictedGroup)
df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample[match(rownames(df), sc$cellNames)]
df$Cluster <- sc$Clusters[match(rownames(df), sc$cellNames)]
df$Group <- sc$predictedGroup[match(rownames(df), sc$cellNames)]
df$Age <- sc$arcsin_ages[match(rownames(df), sc$cellNames)]
df$Age <- as.numeric(as.character(df$Age))
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

for(i in 1:length(cell_types)){
    
    df$GroupInterest <- (df$Group==cell_types[i])*1
    
    gg <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
        stat_summaries_hex(aes(z = GroupInterest, fill = stat(mean)),
        funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
        theme_minimal() + ggtitle(paste0(cell_types[i]))
    
    ggsave(gg, file=paste0("snATACseq/clustering/Plots/",
        "annotation/UMAP_", gsub("/", "-", cell_types[i]), ".pdf"))
    
}

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Age)) + 
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() + scale_colour_viridis_c()
ggsave(g1, file="snATACseq/clustering/Plots/UMAP_ages_arcsin.pdf")

system("mkdir snATACseq/clustering/Plots/clusters")

all_clusters <- unique(df$Cluster)

for(i in 1:length(all_clusters)){
    
    df$ClusterInterest <- (df$Cluster==all_clusters[i])*1
    
    gg <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    stat_summaries_hex(aes(z = ClusterInterest, fill = stat(mean)),
    funs = c('mean'), bins = 100) + scale_fill_viridis_c() + theme_minimal() +
    ggtitle(paste0("Cluster ", all_clusters[i]))
    
    ggsave(gg, file=paste0("snATACseq/clustering/Plots/",
        "clusters/UMAP_Cluster_", all_clusters[i], ".pdf"))
    
}

# 4. return gene matrix 

gene_sc <- getMatrixFromProject(sc, useMatrix = "GeneScoreMatrix")

sce_atac <- SingleCellExperiment(
    list(counts=gene_sc@data$counts@assays@data$GeneScoreMatrix,
    colData=sc@cellColData, rowData=gene_sc@assays@data$counts@elementMetadata))

reducedDim(sce_atac, "UMAP") <- sc@embeddings@listData$UMAP$df
saveRDS(sce_atac, file="snATACseq/processed_data/umap_gene_score_all.rds")

# 5. annotate clusters 
# (rational found in annotation/scatac_cluster_annotation.xlsx)

annotation_cluster <- read_xlsx("annotation/scatac_clustering_annotation.xlsx", 
             sheet="first_clustering")
annotation_cluster$Cluster <- paste0("C", annotation_cluster$Cluster)
sc$Anno <- annotation_cluster$Final_Annotation[match(sc$Clusters, 
   annotation_cluster$Cluster)]

# 4281 low quality cells 

df$Anno <- sc$Anno

col <- c(Astro= '#ffc857', `IN dev`='#c6d5c0', `L2/3`= '#6e0614', L4= '#8b3843', 
    `L5/6`='#a86a72', Micro='#484848',OPC='#92afc2', Oligo='#255f85', 
    `PN dev`='#e2cdd0', `Low quality`='lightgrey', `MGE der`='#b44622',
    `CGE der`='#1c5701', Vas='#a3a3a3')

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Anno)) + geom_point(size=0.2) +
    theme_classic() + scale_color_manual(values=col)
ggsave(g1, file="snATACseq/clustering/Plots/UMAP_annotation_new.pdf")

# save project
saveArchRProject(ArchRProj = sc, load = FALSE)

# 5. remove low quality clusters and redo dimension reduction

sc <- sc[!sc$Anno=="Low quality",]

sc <- addIterativeLSI(ArchRProj = sc, useMatrix = "TileMatrix", 
    name = "IterativeLSI", iterations = 4, clusterParams = list( 
    resolution = 4, sampleCells = 10000, n.start = 10), 
    varFeatures = 50000, dimsToUse = 1:20, force=TRUE, seed=11)

# 6. add batch correction

set.seed(1)
sc <- addHarmony(ArchRProj = sc, reducedDims = "IterativeLSI",
    name = "Harmony", groupBy = "Sample", force=TRUE)

# 7. clustering

sc <- addClusters(input = sc, reducedDims = "Harmony", method = "Seurat",
    name = "Clusters", resolution = 4, force=TRUE, seed = 12, maxClusters = 100)

# 8. UMAP

sc <- addUMAP(ArchRProj = sc, reducedDims = "Harmony", name = "UMAP", 
    nNeighbors = 5, minDist = 0.25, metric = "cosine", force=TRUE,
    seed = 87)

set.seed(10)
sc <- addGeneIntegrationMatrix(ArchRProj = sc, useMatrix = "GeneScoreMatrix", 
    matrixName = "GeneIntegrationMatrix", reducedDims = "Harmony",
    seRNA = serna, addToArrow = TRUE, groupRNA = "sub_clust",
    nameCell = "predictedCell", nameGroup = "predictedGroup",
    nameScore = "predictedScore", force=TRUE)

# save project
saveArchRProject(ArchRProj = sc, load = TRUE, 
    outputDirectory = "snATACseq/clustering_final")

# 9. plots

sc <- loadArchRProject("snATACseq/clustering_final")

df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample[match(rownames(df), sc$cellNames)]
df$Cluster <- sc$Clusters[match(rownames(df), sc$cellNames)]
df$Group <- sc$predictedGroup[match(rownames(df), sc$cellNames)]
df$Age <- sc$arcsin_ages[match(rownames(df), sc$cellNames)]
df$Age <- as.numeric(as.character(df$Age))

colnames(df)[1:2] <- c("UMAP1", "UMAP2")

all_clusters <- unique(df$Cluster)

system("rm -r snATACseq/clustering_final/Plots")
system("mkdir snATACseq/clustering_final/Plots")
system("mkdir snATACseq/clustering_final/Plots/clusters")

p1 <- plotEmbedding(sc, colorBy = "cellColData", name = "Clusters")
plotPDF(p1, name = "Plot-UMAP-RNA-Clusters.pdf", 
        ArchRProj = sc, addDOC = FALSE, width = 5, height = 5)

for(i in 1:length(all_clusters)){
    
    df$ClusterInterest <- (df$Cluster==all_clusters[i])*1
    
    gg <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    stat_summaries_hex( aes(z = ClusterInterest, fill = stat(mean)),
    funs = c('mean'), bins = 100 ) + scale_fill_viridis_c() + theme_minimal() +
    ggtitle(paste0("Cluster ", all_clusters[i]))
    
    ggsave(gg, file=paste0("snATACseq/clustering_final/Plots/",
        "clusters/UMAP_Cluster_", all_clusters[i], ".pdf"))
    
}

cM <- confusionMatrix(sc$Clusters, sc$Anno)
cM <- cM / Matrix::rowSums(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
mapLabs <- cbind(rownames(cM), labelNew)

p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
ggsave(p, file=paste0("snATACseq/clustering_final/Plots/Heatmap_annotation_sub.pdf"))
write.csv(as.matrix(cM), file="snATACseq/clustering_final/Plots/Cluster_anno.csv", 
    quote=F)

cM <- confusionMatrix(sc$Clusters, sc$Sample)
cM <- cM / Matrix::rowSums(cM)
write.csv(as.matrix(cM), file="snATACseq/clustering_final/Plots/Cluster_sample.csv", 
          quote=F)

# 10. reannotate

annotation_cluster <- read_xlsx("annotation/scatac_clustering_annotation.xlsx", 
     sheet="second_clustering")
sc$Anno1 <- annotation_cluster$Final_Annotation[match(sc$Clusters,
    paste0("C", annotation_cluster$Cluster))]

# 11. plots 

df$Anno1 <- sc$Anno1
df <- df[sample(1:nrow(df), nrow(df)),]
g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Anno1)) + geom_point(size=0.2) +
    theme_classic() + scale_color_manual(values=col)
ggsave(g1, file="snATACseq/clustering_final/Plots/UMAP_annotation_new.pdf")


g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2, col=Age)) + 
    geom_point(size=0.2, alpha=0.5) + theme_classic() + 
    scale_colour_viridis_c()
ggsave(g1, file="snATACseq/clustering_final/Plots/UMAP_ages_arcsin.pdf")

cM <- confusionMatrix(sc$Clusters, sc$predictedGroup)
cM <- cM / Matrix::rowSums(cM)
write.csv(as.matrix(cM), 
    file="snATACseq/clustering_final/Plots/Cluster_predicted.csv", quote=F)

# 12. save project
saveArchRProject(ArchRProj = sc, load = FALSE)

# 13. save UMAP, colData and gene logcounts

sc <- loadArchRProject("snATACseq/clustering_final")

gene_sc <- getMatrixFromProject(sc, useMatrix = "GeneScoreMatrix")

sce_atac <- SingleCellExperiment(
    list(counts=gene_sc@data$counts@assays@data$GeneScoreMatrix,
    colData=sc@cellColData, rowData=gene_sc@assays@data$counts@elementMetadata))

reducedDim(sce_atac, "UMAP") <- sc@embeddings@listData$UMAP$df
saveRDS(sce_atac, file="snATACseq/processed_data/umap_gene_score_rm.rds")
