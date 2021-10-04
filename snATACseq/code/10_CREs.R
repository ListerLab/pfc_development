###############################################################################
#                                                                             #
#                                                                             #
# Finding cis regulatory elements                                             #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(ArchR)
library(SingleCellExperiment)
library(reticulate)
library(FNN)
library(parallel)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(edgeR)
library(readxl)

# 1. load ArchR project and all meta data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

addArchRGenome("hg19")
addArchRThreads(threads = 1)

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_annotation")
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)
all_cell_types <- unique(sc$Anno1)

sc$Stages_Types <- paste0(sc$Anno1, "_", sc$Stage)

# 2. find cells for each cell type

ArrowFiles <- getArrowFiles(sc)
Groups <- getCellColData(ArchRProj = sc, 
     select = "Anno1", drop = TRUE)

Cells <- sc$cellNames

cellGroups <- split(Cells, Groups)
cellGroups_new <- list()
for(i in 1:length(cellGroups)){
    if(names(cellGroups)[i] %in% c("L5_6", "L4", "L2_3")){
        a <- c(i, which(names(cellGroups)=="PN_dev"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("MGE_der", "CGE_der")){
        a <- c(i, which(names(cellGroups)=="IN_dev"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("Oligo")){
        a <- c(i, which(names(cellGroups)=="OPC"))
        cellGroups_new[[names(cellGroups)[i]]] <- unlist(cellGroups[a])
    }
    if(names(cellGroups)[i] %in% c("Astro", "Micro", "Vas")){
        cellGroups_new[[names(cellGroups)[i]]] <- cellGroups[[names(cellGroups)[i]]]
    }
    
}
input <- lapply(1:length(cellGroups_new), function(x) 
    lapply(names(ArrowFiles), function(y)
    if(sum(grepl(paste0(y, "#"), cellGroups_new[[x]]))>=40) 
        c(y, names(cellGroups_new)[x]))) 
input <- unlist(input, recursive = FALSE)
input <- input[!sapply(input, is.null)]
input <- lapply(input, function(x) c(x[2], 
    atac_samples$Stage[atac_samples$Sample==x[1]]))
all_stage_types <- unique(sapply(input, function(x) paste0(x[1], "_", x[2])))


# 2. find random representative sample 

set.seed(10)
geosketch <- import('geosketch')
random <- import('random')
random$seed(10)
umap.mat <- as.matrix(sc@embeddings@listData$UMAP$df)
umap.mat <- umap.mat[sc$Stages_Types %in% all_stage_types, ]
sketch.indices <- geosketch$gs(umap.mat, as.integer(400))
cells.select <- rownames(umap.mat)[unlist(sketch.indices)+1]

# 3. find 50 nearest neighbors of every picked atac-seq cell

pca_atac <- sc@reducedDims@listData$Harmony$matDR
query <- pca_atac[cells.select,]
pca_atac <- pca_atac[-match(cells.select, rownames(pca_atac)),]

neighbors <- get.knnx(pca_atac, query, k=50, algorithm="kd_tree")
neighbors <- lapply(1:length(cells.select), function(x) neighbors$nn.index[x,])
names(neighbors) <- cells.select

colInfo.select <- data.frame(cells=cells.select, 
    stages=sc@cellColData[match(cells.select,
             rownames(sc@cellColData)), "Stage"],
    type=sc@cellColData[match(cells.select,
         rownames(sc@cellColData)), "Anno1"]
)

# 4. extract counts in each of the peaks from genes

availableChr <- ArchR:::.availableSeqnames(head(getArrowFiles(sc)))
chromLengths <- getChromLengths(sc)
chromSizes <- getChromSizes(sc)

cellGroups <- lapply(neighbors, function(x) rownames(pca_atac)[x])
cellGroups <- lapply(1:length(cellGroups), function(x) c(names(cellGroups)[x],
        cellGroups[[x]]))


getFrags_cellGroup <- function(cellGroupi) {    
    
    covList <- mclapply(getArrowFiles(sc), function(kkk) {
        tmp <- lapply(availableChr, function(kk){
            ArchR:::.getFragsFromArrow(kkk, 
              chr = kk, out = "GRanges", 
              cellNames = cellGroupi)
        })
        
        tmp <- do.call(c, tmp)
        
    }, mc.cores=5)
    
    covList <- covList[!sapply(covList, function(x) length(x)==0)]
    covList <- as(covList, "GRangesList")
    covList <- unlist(covList, recursive = TRUE, use.names = TRUE)
    
    covList <- c(GRanges(seqnames = seqnames(covList),
       ranges = IRanges(start = start(covList), end = start(covList))),
       GRanges(seqnames = seqnames(covList),
       ranges = IRanges(start = end(covList), end = end(covList))))
    
    covList <- covList + 25 
    covList <- sort(covList)
    
    return(covList)
}

frags_cellGroups <- lapply(cellGroups, function(x) getFrags_cellGroup(x))
saveRDS(frags_cellGroups, file="snATACseq/processed_data/frags_cellGroups.RDS")
saveRDS(cellGroups, file="snATACseq/processed_data/cellGroups.RDS")

# 5. load RNA-seq data 
paths <- "snRNAseq/processed_data"
sce <- readRDS(paste0(paths, 
    "/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"))

# 6. find matching RNA-seq cells and get expression (cpm)

cellGroupsRNA <- lapply(cellGroups, function(x) 
    sc$predictedCell[match(x, rownames(sc@cellColData))])

all_stages_rna <- lapply(cellGroupsRNA, function(x) 
    sce$stage_ids[match(x, colnames(sce))])
all_stages_rna <- lapply(all_stages_rna, function(x) {
    tmp <- sort(table(x), decreasing=TRUE)
    names(tmp)[1]
})
all_stages_rna <- unlist(all_stages_rna)
colInfo.select$stages_rna <- all_stages_rna

saveRDS(colInfo.select, file="snATACseq/processed_data/colInfo.select.RDS")

genes_per_groups <- sapply(cellGroupsRNA, function(x) rowSums(counts(sce[,
            match(x, colnames(sce))])))
rownames(genes_per_groups) <- rowData(sce)$gene_ids
genes_per_groups <- edgeR::cpm(genes_per_groups)

rm(sce)
gc()

# 7. get locations for genes and 250kb window around TSS
all_peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_gr.Rds")
all_peaks <- as(all_peaks, "GRangesList")
all_peaks <- unlist(all_peaks, recursive = TRUE)

genes <- promoters(txdb, filter=list(gene_id=rownames(genes_per_groups)),
   upstream=250000, downstream=250000, columns="GENEID")
strand(genes) <- "*"

genes <- promoters(txdb,
                   upstream=250000, downstream=250000, columns="GENEID")

ind <- split(1:length(genes), as.character(genes@elementMetadata$GENEID))
genes <- mclapply(ind, function(x) genes[x], mc.cores=5)

peaks_in_genes <- mclapply(genes, function(x) {
    
    unique(all_peaks[findOverlaps(x, all_peaks)@to])
    
}, mc.cores=10)


# 8. find background peak genes (not on the same chromosome)
set.seed(10)
peaks_in_background <- mclapply(genes, function(x){
    
    tmp <- all_peaks[!seqnames(all_peaks)==as.character(seqnames(x[1]))]
    sample(tmp, 10)
}, mc.cores=10)


# 9. find overlaps of peaks with gene windows and find fragments in peaks (cpm)

frags_cellGroups <- readRDS("snATACseq/processed_data/frags_cellGroups.RDS")
lib_sizes <- lengths(frags_cellGroups)
     
getFrags_genes <- function(frags_cellGroups, peaks_in_genesi, lib_sizes) {
        
    frags_in_peaks <- sapply(frags_cellGroups, function(xx)
        countOverlaps(peaks_in_genesi, xx))
    if(class(frags_in_peaks)[[1]]!="matrix"){
        frags_in_peaks <- matrix(frags_in_peaks, 
            nrow=length(peaks_in_genesi), ncol=length(frags_cellGroups))
    }
        
    return(edgeR::cpm(frags_in_peaks, lib.size=lib_sizes))
}

frags_in_peaks <- mclapply(peaks_in_genes, function(a)
     getFrags_genes(frags_cellGroups, a, lib_sizes), mc.cores=30)
frags_in_background <- mclapply(peaks_in_background, function(a)
    getFrags_genes(frags_cellGroups, a, lib_sizes), mc.cores=15)
saveRDS(frags_in_peaks, file="snATACseq/processed_data/frags_in_peaks.RDS")
saveRDS(frags_in_background, file="snATACseq/processed_data/frags_in_background.RDS")

#10. find correlation between 
rna_genes_pseudo_bulk <- genes_per_groups[names(frags_in_peaks),]
rna_genes_pseudo_bulk <- split(rna_genes_pseudo_bulk, 
                seq(nrow(rna_genes_pseudo_bulk)))
names(rna_genes_pseudo_bulk) <- names(frags_in_peaks)

ind <- sapply(frags_in_peaks, function(x) dim(x)[1])==0
frags_in_peaks <- frags_in_peaks[!ind]
frags_in_background <- frags_in_background[!ind]
rna_genes_pseudo_bulk <- rna_genes_pseudo_bulk[!ind]

ind1 <- unlist(sapply(frags_in_background, function(x) is.null(x)))
frags_in_background <- frags_in_background[!ind1]

pearson_corr <- mcmapply(function(X,Y) {
    t(apply(X, 1, function(x) c(cor(x,Y), cor.test(x, Y)$p.value)))
}, X=frags_in_peaks, Y=rna_genes_pseudo_bulk, mc.cores=10)

pearson_corr_background <- mcmapply(function(X,Y) {
    t(apply(X, 1, function(x) c(cor(x,Y), cor.test(x, Y)$p.value)))
}, X=frags_in_background, Y=rna_genes_pseudo_bulk[!ind1],
mc.cores=10, SIMPLIFY = FALSE)

# 11. put it together
clean_up <- function(pearson_corr, all_peaks){
    
    for(i in 1:length(pearson_corr)){
        
        pearson_corr[[i]] <- cbind(names(pearson_corr)[i], pearson_corr[[i]])
        
    }
    
    pearson_corr <- do.call(rbind, pearson_corr)
    colnames(pearson_corr) <- c("Gene", "Correlation", "PValue")
    pearson_corr <- as.data.frame(pearson_corr)
    pearson_corr$PValue <- as.numeric(as.character(pearson_corr$PValue))
    pearson_corr$Correlation <- as.numeric(as.character(pearson_corr$Correlation))
    ind <- is.na(pearson_corr$PValue)
    pearson_corr <- pearson_corr[!ind,]
    pearson_corr$FDR <- p.adjust(pearson_corr$PValue, method="fdr")
    pearson_corr$peak_name <- sapply(rownames(pearson_corr), function(y) 
        paste0(strsplit(y, ".", fixed=T)[[1]][1:2], collapse="."))
    ind <- is.na(match(pearson_corr$peak_name, names(all_peaks)))
    pearson_corr <- pearson_corr[!ind,]
    pearson_corr$seqnames <- as.character(seqnames(all_peaks)[
        match(pearson_corr$peak_name, names(all_peaks))])
    pearson_corr$start <- start(all_peaks)[
        match(pearson_corr$peak_name, names(all_peaks))]
    pearson_corr$end <- end(all_peaks)[
        match(pearson_corr$peak_name, names(all_peaks))]
    return(pearson_corr)
}

pearson_corr <- clean_up(pearson_corr, all_peaks)
pearson_corr$name <- paste0(pearson_corr$Gene, "_", pearson_corr$peak_name) 
pearson_corr <- pearson_corr[!duplicated(pearson_corr$name),]
pearson_corr_background <- clean_up(pearson_corr_background, all_peaks)

saveRDS(pearson_corr, file="snATACseq/processed_data/pearson_corr.RDS")
saveRDS(pearson_corr_background, 
    file="snATACseq/processed_data/pearson_corr_background.RDS")
pearson_corr <- readRDS("snATACseq/processed_data/pearson_corr.RDS")
pearson_corr_background <- readRDS("snATACseq/processed_data/pearson_corr_background.RDS")

# 12. find cutoff

nn <- 0.025*nrow(pearson_corr_background)
cutoff <-  sort(abs(pearson_corr_background$Correlation), decreasing = T)[nn+1]
rm(pearson_corr_background)
gc()

# 13. filter for significant correlations

pearson_corr_sig <- pearson_corr[(pearson_corr$FDR<0.05 &
                 abs(pearson_corr$Correlation)>cutoff), ]

tss <- promoters(txdb, filter=list(gene_id=pearson_corr_sig$Gene),
                 upstream=0, downstream=0, columns="GENEID")
strand(tss) <- "*"

tmp <- GRanges(seqnames=pearson_corr_sig$seqnames, IRanges(start=
    as.numeric(pearson_corr_sig$start), end=as.numeric(pearson_corr_sig$end)), 
    Gene=pearson_corr_sig$Gene)
strand(tmp) <- "*"

# 13. find distance to the nearest associated gene

find_distance <- function(tmptmp, tss){
    
    genes <- tss[unlist(tss$GENEID==tmptmp$Gene)]
    dist <- distanceToNearest(genes, tmptmp)
    min(dist@elementMetadata$distance)
    
}

min_dist <- mclapply(1:length(tmp), function(x) find_distance(tmp[x], tss), 
        mc.cores=15)

pearson_corr_sig$dist <- unlist(min_dist)

saveRDS(pearson_corr_sig, file="snATACseq/processed_data/pearson_corr_sig.RDS")

# 14. get out data for plotting 

frags_in_peaks <- readRDS("snATACseq/proccessed_data/frags_in_peaks.RDS")
saveRDS(rna_genes_pseudo_bulk, "snATACseq/processed_data/rna_pseudo_bulk.RDS")
rna_genes_pseudo_bulk <- readRDS("snATACseq/processed_data/rna_pseudo_bulk.RDS")

pearson_corr_sig$Gene <- as.character(pearson_corr_sig$Gene)
genes_sam <- unique(pearson_corr_sig$Gene)
peaks_sam <- split(rownames(pearson_corr_sig), pearson_corr_sig$Gene)
peaks_sam <- lapply(peaks_sam, function(x) sapply(x, function(y) 
    paste0(strsplit(y, ".", fixed=T)[[1]][1:2], collapse=".")))
peaks_sam <- lapply(peaks_sam, function(x) unique(x))

sam_rna <- do.call(rbind, 
    rna_genes_pseudo_bulk[match(genes_sam, names(rna_genes_pseudo_bulk))])

sam_atac <- lapply(1:length(genes_sam), function(i)
    colSums(frags_in_peaks[[genes_sam[i]]][peaks_sam[[i]],, drop=F]))
sam_atac <- do.call(rbind, sam_atac)
rownames(sam_atac) <- rownames(sam_rna) <- genes_sam

plot_list <- list(atac=sam_atac, rna=sam_rna, info=colInfo.select)
saveRDS(plot_list, file="snATACseq/processed_data/plot_atac_rna_corr.RDS")

