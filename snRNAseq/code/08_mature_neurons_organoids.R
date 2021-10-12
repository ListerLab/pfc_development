library(scater)
library(scran)
library(SingleCellExperiment)
library(batchelor)
library(FNN)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)

path1 <- "snRNAseq/processed_data/orgpred/"
all_samples <- c("Luc9228.RDS", "org2290.RDS", "org2432.RDS")

source("snATACseq/R/functions_nmf.R")

# Compare mature to immature neuronal organoid cells

## Comparison in individual datasets

all_markers_PN <- list()
sce <- list()

for(i in all_samples){
    
    sce[[i]] <- readRDS(paste0(path1, i))
    sce[[i]] <- logNormCounts(sce[[i]])

    
    ind_PN <- sce[[i]]$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
        "L5/6_THEMIS_TLE4", "PN_dev")
    ind_n_fetal <- sce$predstage != "Fetal"

    PN_cat <- rep(NA, ncol(sce[[i]]))
    PN_cat[ind_PN & ind_n_fetal] <- "mature"
    PN_cat[ind_PN & !ind_n_fetal] <- "immature"
    markers <- findMarkers(sce[[i]], groups=PN_cat, row.data=rowData(sce[[i]]))
    all_markers_PN[[i]] <- markers$mature
    
}

all_markers_PN_up <- lapply(all_markers_PN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC>0])
table_up_PN <- sort(table(unlist(all_markers_PN_up)))

all_markers_PN_down <- lapply(all_markers_PN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC<0])
table_down_PN <- sort(table(unlist(all_markers_PN_down)))

## Comparison in joined datasets

all_sce <- lapply(all_samples, function(x) readRDS(paste0(path1, x)))
all_sce <- multiBatchNorm(all_sce)

all_sce <- do.call(cbind, all_sce)

ind_PN <- all_sce$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
                                  "L5/6_THEMIS_TLE4", "PN_dev")
ind_n_fetal <- all_sce$predstage != "Fetal"

PN_cat <- rep(NA, ncol(all_sce))
PN_cat[ind_PN & ind_n_fetal] <- "mature"
PN_cat[ind_PN & !ind_n_fetal] <- "immature"
markers_PN <- findMarkers(all_sce, groups=PN_cat, row.data=rowData(all_sce), 
    block=all_sce$Sample)
markers_PN <- markers_PN$mature

ind_PN_down <- match(names(table_down_PN), markers_PN$Symbol)
markers_PN$down_ind <- NA
markers_PN$down_ind[ind_PN_down] <- table_down_PN
ind_PN_up <- match(names(table_up_PN), markers_PN$Symbol)
markers_PN$up_ind <- NA
markers_PN$up_ind[ind_PN_up] <- table_up_PN

write.table(markers_PN, file="snRNAseq/processed_data/Organoids_PN_markers.csv",
            quote=F, sep=",")


# Run permutation DEG

seeds <-  sample(1:1000000, 100)

nu_mat <- sum(PN_cat=="mature", na.rm=TRUE)
nu_immat <- sum(PN_cat=="immature", na.rm=TRUE)

all_sam_markers <- list()

for(i in 1:length(seeds)){
    
    set.seed(seeds[i])
    
    ind <- which(!is.na(PN_cat))
    ind_mat_sam <- sample(ind, nu_mat)
    ind_immat_sam <- setdiff(ind, ind_mat_sam)
    PN_sam <- PN_cat
    PN_sam[ind_mat_sam] <- "mature"
    PN_sam[ind_immat_sam] <- "immature"
    markers_tmp <- findMarkers(all_sce, groups=PN_sam, row.data=rowData(all_sce), 
        block=all_sce$Sample)$mature
    all_sam_markers[[i]] <-markers_tmp
}

lapply(all_sam_markers, function(x) sum(x$FDR<0.05 & x$summary.logFC>0))
all_gos_up_sam <- lapply(all_sam_markers, function(x) {
    tmp <- x$ID[x$FDR<0.05 & x$summary.logFC>0]
    go_term_analysis(tmp, markers_PN$ID)@result })
lapply(all_gos_up_sam, function(x)
    x$Description[x$p.adjust<0.05])
