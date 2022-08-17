###############################################################################
#                                                                             #
#                                                                             #
# Regions for disease gene trends motif enrichment                            #
#                                                                             #
#                                                                             #    
###############################################################################

library(GenomicRanges)
library(readxl)
library(ggplot2)
library(SingleCellExperiment)
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)
library(ArchR)
library(chromVAR)

source("snATACseq/R/functions_TF.R")

# 1. read in data

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

paths <- "snRNAseq/processed_data"

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")

for(i in 1:length(peaks)){
    
    names(peaks[[i]]) <- paste0(names(peaks)[i], ".", names(peaks[[i]]))
}

sce <- readRDS(paste0(paths,
    "/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"))


# 2. remove all promoters regions

peaks <- lapply(peaks, function(x) {
    hits <- findOverlaps(x, promoters)
    ind <- setdiff(1:length(x), unique(hits@from))
    x[ind]
})

# 3. load disease genes
all_files <- list.files("scrna/processed_data/disease_enrichment_genes_sub")
all_diseases <- lapply(all_files, function(x) read.csv(paste0(paths,
    "/disease_enrichment_genes_sub/", x)))
all_diseases <- lapply(all_diseases, function(x) x[,2])
names(all_diseases) <- gsub(".csv", "", all_files)

all_diseases <- lapply(all_diseases, function(x)
    rowData(sce)$gene_ids[match(x, rowData(sce)$index)])

cell_types <- matrix(c("Astro", "Astro",
                       "L2-3_CUX2", "L2_3",
                       "L5-6_THEMIS", "L5_6",
                       "L5-6_TLE4", "L5_6",
                       "L4_RORB", "L4",
                       "VIP", "CGE_der",
                       "LAMP5_CA1", "CGE_der",
                       "ID2", "CGE_der",
                       "SST", "MGE_der",
                       "PV", "MGE_der",
                       "PV_SCUBE3", "MGE_der",
                       "Oligo", "Oligo",
                       "OPC", "Oligo",
                       "Micro", "Micro"
), ncol=2, byrow=T)

# 4. get peaks for disease, gene trend, cell type combinations

cre_diseases <- mapply(function(X,Y) obtain_peak_ids_dis(X,Y, peaks, cell_types),
       X=all_diseases, Y=names(all_diseases))
names_cre_diseases <- sapply(names(cre_diseases), function(x) 
  strsplit(x, "..", fixed=T)[[1]][2])
names_cre_diseases <- sapply(names_cre_diseases, function(x) 
  strsplit(x, "_", fixed=T)[[1]])
names_cre_diseases <- sapply(names_cre_diseases, function(x) 
  paste0(x[1:(length(x)-1)], collapse="_"))
new_names_cre_diseases <- cell_types[match(names_cre_diseases, 
  cell_types[,1]),2]

index <- split(1:length(new_names_cre_diseases), new_names_cre_diseases)

region_interest <- list()
for(i in 1:length(index)){
  
  tmp_tmp <-   cre_diseases[index[[i]]]
  #tmp_names <- sapply(names(tmp_tmp), function(x)
  #    strsplit(x, "..", fixed=T)[[1]][2])
  #index1 <- split(1:length(tmp_names), tmp_names)
  #tmp_tmp <- lapply(index1, function(x) unlist(tmp[x]))
  #names(tmp_tmp) <- paste0(names(index)[i], "..", names(tmp_tmp))
  
  region_interest[[names(index)[i]]] <- tmp_tmp
  
}

region_interest <- lapply(region_interest, function(x) unlist(x))

saveRDS(region_interest, file="snATACseq/processed_data/tf_diseases_new_sub_gr.RDS")

region_interest <- readRDS("snATACseq/processed_data/tf_diseases_new_sub_gr.RDS")

# 6. load ArchR project

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")
sc$Anno1 <- gsub(" ", "_", sc$Anno1)
sc$Anno1 <- gsub("/", "_", sc$Anno1)
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")
peaks <- peaks[-2] # remove vasculature

for(i in 1:length(peaks)){
    
    peaks[[i]]@elementMetadata <- peaks[[i]]@elementMetadata[, c(stages, "CRE",
          "peakType")]
    peaks[[i]]@elementMetadata$cellType <- names(peaks)[i]
}

# 7. enrichment for each cell type

enriched_diseases <- list()
for(i in 1:length(region_interest)){
    
    sc_sub <- sc[sc$Anno1==names(region_interest)[i],]
    
    sc_sub <- addPeakSet(sc_sub, peaks[[names(region_interest)[i]]], force=TRUE)
    sc_sub <- addPeakMatrix(sc_sub, force=TRUE)
    sc_sub <- addMotifAnnotations(ArchRProj = sc_sub, motifSet = "JASPAR2018", 
                          name = "Motif", force=TRUE)
    set.seed(10)
    bgdPeaks <- SummarizedExperiment::assay(getBgdPeaks(sc_sub, method="ArchR"))
    matches <- getMatches(sc_sub, NULL)
    
    enriched_diseases[[i]] <- lapply(region_interest[[i]], function(x)
      compute_enrich_traj(x, bgdPeaks, matches, 
                      names(region_interest)[i]))
}

names(enriched_diseases) <- names(region_interest)

saveRDS(enriched_diseases, file="snATACseq/processed_data/tf_diseases_sub_new.rds")
enriched_diseases <- readRDS("snATACseq/processed_data/tf_diseases_new.rds")

for(i in 1:length(enriched_diseases)){
    
  for(j in 1:length(enriched_diseases[[i]])){
        
      enriched_diseases[[i]][[j]]$trajectory <- names(enriched_diseases[[i]])[j]
      enriched_diseases[[i]][[j]]$pvalue <- 
        10^(-enriched_diseases[[i]][[j]]$mlog10p)
  }
    
  enriched_diseases[[i]] <- do.call(rbind, enriched_diseases[[i]])
  enriched_diseases[[i]]$adjusted_pvalue <- p.adjust(
    enriched_diseases[[i]]$pvalue, "fdr")
}

test <- lapply(enriched_diseases, function(x) x)
test <- do.call(rbind, test)

test$general_trend <- sapply(test$trajectory, function(x) 
    strsplit(x, "..", fixed=TRUE)[[1]][1])
test$pathway <- sapply(test$trajectory, function(x) 
    strsplit(x, "..", fixed=TRUE)[[1]][2])
test$pathway <- sapply(test$pathway, function(x) {
    tmp <- strsplit(x, "_", fixed=TRUE)[[1]]
    tmp[length(tmp)]
    })
test$trajectory <- sapply(test$trajectory, function(x) 
    strsplit(x, "..", fixed=TRUE)[[1]][2])
test$trajectory <- sapply(test$trajectory, function(x) {
    tmp <- strsplit(x, "_", fixed=TRUE)[[1]]
    paste0(tmp[1:(length(tmp)-1)], collapse="_")
  })
rownames(test) <- NULL


data.table::fwrite(test, 
  file="snATACseq/processed_data/tf_diseases_new.tsv", 
    row.names=FALSE, quote=FALSE, sep="\t")





