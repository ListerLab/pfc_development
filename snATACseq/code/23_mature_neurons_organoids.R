library(scater)
library(scran)
library(SingleCellExperiment)
library(batchelor)
library(FNN)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)

path1 <- "paths_to_organoids"
all_samples <- c("Luc9228.RDS", "org2290.RDS", "org2432.RDS")

source("R/functions_nmf.R")

# Compare mature to immature neuronal organoid cells

## Comparison in individual datasets

all_markers_PN <- list()
all_markers_IN <- list()

for(i in all_samples){
    
    sce <- readRDS(paste0(path1, i))
    sce <- logNormCounts(sce)

    
    ind_PN <- sce$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
        "L5/6_THEMIS_TLE4", "PN_dev")
    ind_n_fetal <- sce$predstage != "Fetal"
    ind_IN <- sce$predcelltype %in% c("CGE_dev", "ID2", "MGE_dev", "PV",
        "SST", "VIP")
    
    PN_cat <- rep(NA, ncol(sce))
    PN_cat[ind_PN & ind_n_fetal] <- "mature"
    PN_cat[ind_PN & !ind_n_fetal] <- "immature"
    markers <- findMarkers(sce, groups=PN_cat, row.data=rowData(sce))
    all_markers_PN[[i]] <- markers$mature
    
    IN_cat <- rep(NA, ncol(sce))
    IN_cat[ind_IN & ind_n_fetal] <- "mature"
    IN_cat[ind_IN & !ind_n_fetal] <- "immature"
    markers <- findMarkers(sce, groups=IN_cat, row.data=rowData(sce))
    all_markers_IN[[i]] <- markers$mature
}

all_markers_PN_up <- lapply(all_markers_PN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC>0])
table_up_PN <- sort(table(unlist(all_markers_PN_up)))

all_markers_PN_down <- lapply(all_markers_PN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC<0])
table_down_PN <- sort(table(unlist(all_markers_PN_down)))

all_markers_IN_up <- lapply(all_markers_IN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC>0])
table_up_IN <- sort(table(unlist(all_markers_IN_up)))

all_markers_IN_down <- lapply(all_markers_IN, function(x)
    x$Symbol[x$FDR<0.05 & x$summary.logFC<0])
table_down_IN <- sort(table(unlist(all_markers_IN_down)))

## Comparison in joined datasets

all_sce <- lapply(all_samples, function(x) readRDS(paste0(path1, x)))
all_sce <- multiBatchNorm(all_sce)

all_sce <- do.call(cbind, all_sce)

ind_PN <- all_sce$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
                                  "L5/6_THEMIS_TLE4", "PN_dev")
ind_n_fetal <- all_sce$predstage != "Fetal"
ind_IN <- all_sce$predcelltype %in% c("CGE_dev", "ID2", "MGE_dev", "PV",
                                  "SST", "VIP")

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

write.table(markers_PN, file="snATACseq/processed_data/Organoids_PN_markers.csv",
            quote=F, sep=",")

IN_cat <- rep(NA, ncol(all_sce))
IN_cat[ind_IN & ind_n_fetal] <- "mature"
IN_cat[ind_IN & !ind_n_fetal] <- "immature"
markers_IN <- findMarkers(all_sce, groups=IN_cat, row.data=rowData(all_sce), 
    block=all_sce$Sample)
markers_IN <- markers_IN$mature

ind_IN_down <- match(names(table_down_IN), markers_IN$Symbol)
markers_IN$down_ind <- NA
markers_IN$down_ind[ind_IN_down] <- table_down_IN
ind_IN_up <- match(names(table_up_IN), markers_IN$Symbol)
markers_IN$up_ind <- NA
markers_IN$up_ind[ind_IN_up] <- table_up_IN

write.table(markers_IN, file="snATACseq/processed_data/Organoids_IN_markers.csv",
    quote=F, sep=",")



