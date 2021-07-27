###############################################################################
#                                                                             #
#                                                                             #
# Integrate organoids                                                         #
#                                                                             #
#                                                                             #    
###############################################################################

library(ArchR)
library(uwot)
library(readxl)
library(FNN)

addArchRGenome("hg19")
addArchRThreads(threads = 10)

# 1. read in data make ArrowFile

link <- "path_to_organoid_scatac"
atac_samples_org <- c("RL1738")

for(i in 1:length(atac_samples_org)){
    
    barcodes <- getValidBarcodes(csvFiles = paste0(link, atac_samples_org[i], 
        "/outs/singlecell.csv"), 
        sampleNames = paste0(link, atac_samples_org[i], 
         "/outs/filtered_peak_bc_matrix/barcodes.tsv"))
    
    ArrowFiles <- createArrowFiles(
        inputFiles = paste0(link, atac_samples_org[i], 
            "/outs/fragments.tsv.gz"),
        sampleNames = atac_samples_org[i],
        filterTSS = 0,
        filterFrags = 0, 
        validBarcodes = barcodes[[1]],
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
        force=TRUE,
        offsetPlus = 0,
        offsetMinus = 0
    )
}

# 2. Find doublets

ArrowFiles <- paste0(atac_samples_org, ".arrow")

set.seed(1)
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 20, 
    knnMethod = "UMAP", 
    LSIMethod = 1,
    force=TRUE
)

# 3. Make ArchR project

atac_samples <- read_xlsx("../neuronal_mat/annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

ArrowFiles <- c(paste0(atac_samples_org, ".arrow"),
    paste0("../neuronal_mat/arrows_good/", atac_samples$Sample, ".arrow"))

sc <- ArchRProject(
    ArrowFiles = paste0(ArrowFiles), 
    outputDirectory = "combine_RL1738",
    copyArrows = FALSE 
)

sc <- filterDoublets(sc)

# 4. Load reference UMAP

umap_org <- readRDS("snATACseq/processed_data/UMAP_GeneScore_all.rds")
colnames(umap_org[[3]]) <- c("UMAP1", "UMAP2")
umap_org[[3]]$dataset <- "ref"

# 5. Dimension reduction 

sc <- addIterativeLSI(
    ArchRProj = sc,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( 
        resolution = 4, 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:20,
    force=TRUE,
    seed=11
)

# 6. Add batch correction

sc <- addHarmony(
    ArchRProj = sc,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force=TRUE,
    max.iter.harmony = 100
)

# 7. UMAP

sc <- addUMAP(ArchRProj = sc, reducedDims = "Harmony", name = "UMAP", 
              nNeighbors = 3, minDist = 0.08, metric = "cosine",
              force=TRUE, seed = 13)

df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample[match(rownames(df), sc$cellNames)]
df$dataset <- "ref"
df$dataset[df$Sample=="RL1738"] <- "query"
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

g1 <- ggplot(df[rev(1:nrow(df)),], aes(x=UMAP1, y=UMAP2, col=dataset)) + 
    geom_point(size=0.2, alpha=0.5) + theme_classic()
ggsave(g1, file="UMAP_all.png", height=7, width=7)

# 8. Project into reference UMAP

ind <- match(rownames(umap_org[[3]]), rownames(sc))
reference_umap = umap(sc@reducedDims$Harmony$matDR[ind,], n_neighbors = 5, 
   ret_model = TRUE, a = NULL, b = NULL, init=as.matrix(umap_org[[3]][,1:2]),
   spread = 1, min_dist = 0.08, verbose = TRUE, n_epochs = 0)

qumap = umap_transform(sc@reducedDims$Harmony$matDR[sc$Sample %in% 
    atac_samples_org,], reference_umap)

# 9. Save and find low quality cells

source("snATACseq/R/predict_functions.R")

qumap = as.data.frame(qumap)
colnames(qumap) <- c("UMAP1", "UMAP2")
qumap$dataset <- "query"
rownames(qumap) <- rownames(sc[sc$Sample %in% atac_samples_org,])
commap = rbind(umap_org[[3]], qumap)
commap$Sample <- sapply(rownames(commap), function(x) 
    strsplit(x, "#", fixed=T)[[1]][1])
ind <- match(rownames(commap[commap$dataset=="ref",]), rownames(umap_org[[3]]))
commap$Cluster <- NA
commap$Cluster[commap$dataset=="ref"] <- umap_org[[2]]$Clusters[ind]

cluster_annotation <- read_xlsx("annotation/clustering_annotation.xlsx",
    sheet="First clustering")
commap$Anno <- NA
commap$Anno[commap$dataset=="ref"] <- cluster_annotation$Annotation[match(
    commap$Cluster[commap$dataset=="ref"], cluster_annotation$Cluster)]

commap$Anno[commap$dataset=="query"]  <- predict_cells(commap, 
    train=commap$dataset=="ref", test=commap$dataset=="query", knn=1, 
    col="Anno")

saveRDS(commap, file="snATACseq/processed_data/Organoids_Integrated_All.RDS")


# 10. Delete low quality cells

delete_lq <- rownames(commap)[commap$Anno=="Low quality"]
sc <- sc[-match(delete_lq, rownames(sc)),]

# 11. Dimension Reduction

sc <- addIterativeLSI(
    ArchRProj = sc,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( 
        resolution = 4, 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 50000, 
    dimsToUse = 1:20,
    force=TRUE,
    seed=11
)

# 12. Add batch correction

sc <- addHarmony(
    ArchRProj = sc,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force=TRUE,
    max.iter.harmony = 100
)

# 13. UMAP

sc <- addUMAP(ArchRProj = sc, reducedDims = "Harmony", name = "UMAP", 
              nNeighbors = 3, minDist = 0.08, metric = "cosine",
              force=TRUE, seed = 13)

df <- sc@embeddings@listData$UMAP$df
df$Sample <- sc$Sample[match(rownames(df), sc$cellNames)]
df$dataset <- "ref"
df$dataset[df$Sample=="RL1738"] <- "query"
colnames(df)[1:2] <- c("UMAP1", "UMAP2")

g1 <- ggplot(df[rev(1:nrow(df)),], aes(x=UMAP1, y=UMAP2, col=dataset)) + 
    geom_point(size=0.2, alpha=0.5) + theme_classic()
ggsave(g1, file="UMAP_rm.png", height=7, width=7)

# 14. Integrate into reference UMAP without low quality cells

umap_org <- readRDS("snATACseq/processed_data/UMAP_GeneScore_rm.rds")
colnames(umap_org[[3]]) <- c("UMAP1", "UMAP2")
umap_org[[3]]$dataset <- "ref"

ind <- match(rownames(umap_org[[3]]), rownames(sc))
reference_umap = umap(sc@reducedDims$Harmony$matDR[ind,], n_neighbors = 5, 
    ret_model = TRUE, a = NULL, b = NULL, init=as.matrix(umap_org[[3]][,1:2]),
     spread = 1, min_dist = 0.08, verbose = TRUE, n_epochs = 0)

qumap = umap_transform(sc@reducedDims$Harmony$matDR[sc$Sample %in% 
  atac_samples_org,], reference_umap)

# 15. Save data

qumap = as.data.frame(qumap)
colnames(qumap) <- c("UMAP1", "UMAP2")
qumap$dataset <- "query"
rownames(qumap) <- rownames(sc[sc$Sample %in% atac_samples_org,])
commap = rbind(umap_org[[3]], qumap)
commap$Sample <- sapply(rownames(commap), function(x) 
    strsplit(x, "#", fixed=T)[[1]][1])
ind <- match(rownames(commap[commap$dataset=="ref",]), rownames(umap_org[[3]]))
commap$Cluster <- NA
commap$Cluster[commap$dataset=="ref"] <- umap_org[[2]]$Clusters[ind]

cluster_annotation <- read_xlsx("annotation/clustering_annotation.xlsx",
     sheet="Second clustering")
commap$Anno <- NA
commap$Anno[commap$dataset=="ref"] <- cluster_annotation$Annotation[match(
    commap$Cluster[commap$dataset=="ref"], cluster_annotation$Cluster)]

saveRDS(commap, file="snATACseqprocessed_data/Organoids_Integrated_RL1738.RDS")

