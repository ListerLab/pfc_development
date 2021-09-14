# Make network

Sys.setenv(RETICULATE_MINICONDA_PATH="/home/sfreytag/working_data_03/programmes/razor/miniconda3")
library(scater)
library(qusage)
library(parallel)
library(reticulate)

use_condaenv("reticulate")
np <- import("numpy")
sklearn <- import("sklearn.linear_model")

motif_positions <- readRDS("snATACseq/clustering_annotation/Annotations/Motif-Positions-In-Peaks.rds")

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
peaks <- peaks[-2]
for(i in 1:length(peaks)){
    
    ind <- is.na(peaks[[i]]$CRE)
    peaks[[i]] <- peaks[[i]][!ind]
}

peaks_tmp <- peaks[[1]]
motif_tmp <- motif_positions[[1]]
motif_name <- names(motif_positions)[1]

make_network_motif <- function(peaks_tmp, motif_tmp, motif_name){
    
    ind <- findOverlaps(peaks_tmp, motif_tmp)
    cres <- peaks_tmp$CRE[ind@from]
    cres <- sapply(cres, function(x) strsplit(as.character(x), 
        "|", fixed=T)[[1]])
    cres <- unlist(cres)
    names(cres) <- NULL
    data.frame(CRE=unique(cres), motif=motif_name)
    
    
}

sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

all_networks <- list()

for(ii in 1:length(peaks)){

    network <- lapply(1:length(motif_positions), function(i)
        make_network_motif(peaks[[ii]], motif_positions[[i]], 
        names(motif_positions)[i]))
    network <- do.call(rbind, network)
    network$motif <- sapply(network$motif, function(x) 
        strsplit(x, "_", fixed=T)[[1]][1])
    network$CRE <- rowData(sce)$index[match(network$CRE, rowData(sce)$gene_ids)]
    network$motif <- sapply(network$motif, function(x) 
        strsplit(x, ".var", fixed=T)[[1]][1]) 
    all_networks[[ii]] <- network
    
    
}

names(all_networks) <- names(peaks)

saveRDS(all_networks, file="snATACseq/processed_data/All_networks.Rds")

cell_types <- c(L4="L4_RORB", L5_6="L5/6_TL4", L5_6="L5/6_THEMIS",
                L2_3="L2/3_CUX2")

run_linear_model <- function(target_gene, all_networks, sce) {
    
    predictors <- all_networks[all_networks[,1] %in% target_gene, 2]
    predictors <- setdiff(predictors, target_gene)
    tmp <- logcounts(sce[rowData(sce)$index %in% c(target_gene, predictors),])
    tmp <- t(tmp)
    x <- as.matrix(tmp[,-match(target_gene, colnames(tmp)), drop=FALSE])
    y <- as.vector(tmp[,target_gene,])
    predictors <- colnames(x)
    model_br = sklearn$BayesianRidge()
    model_br$fit(x, y)
    coef_mean = model_br$coef_
    coef_var = diag(model_br$sigma_)
    pvalue <- pnorm(abs(coef_mean), mean=0, sd=sqrt(coef_var), lower.tail=FALSE)*2
    
    if(sum(pvalue<0.05)>0){
        
        res <- data.frame(tfs=predictors[pvalue<0.05], target=target_gene,
                          strength=coef_mean[pvalue<0.05], pvalue=pvalue[pvalue<0.05])
        return(res)
        
    } else {
        
        return(NULL)
    }
}


for(i in 1:length(cell_types)){

    all_networks <- readRDS("snATACseq/processed_data/All_networks.Rds")
    all_networks <- all_networks[[names(cell_types)[i]]]
    tf_cell_type <- readRDS("snATACseq/processed_data/tf_cell_type_rna.RDS")
    tf_cell_type <- lapply(tf_cell_type, function(x) names(x)[!x==0])
    tf_cell_type <- tf_cell_type[[names(cell_types)[i]]]
    gc()
    
    all_networks[,2] <- toupper(all_networks[,2])
    tmp_tf <- unique(all_networks[,2])
    tmp_tf <- sapply(tmp_tf, function(x) strsplit(x, "..", fixed=T)[[1]])
    tmp_tf$`EWSR1.FLI1` <- c("EWSR1", "FLI1")

    ind <- mclapply(tmp_tf, function(x) match(x, 
        tf_cell_type), mc.cores=2)
    tf_not<- names(ind)[sapply(ind, function(x) any(is.na(x)))]
    all_networks <- all_networks[!all_networks[,2] %in% tf_not,]
    tmp <- mclapply(1:dim(all_networks)[1], function(x) {
        a <- strsplit(all_networks[x,2], "..", fixed=T)[[1]]
        cbind(all_networks[x,1], a)
    }, mc.cores=10)
    tmp <- do.call(rbind, tmp)
    all_networks <- tmp

    rm(tmp)
    gc()

    go_terms <- read.gmt("annotation/hsapiens.GO:BP.name.gmt")
    gos <- c(Ensheatment="GO:0007272")
    ind <- which(names(go_terms) %in% gos)
    go_ids <- go_terms[ind][[1]]

    sce <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")
    sce <- sce[, sce$major_clust==cell_types[i]]
    sce <- logNormCounts(sce)

    ind <- all_networks[,1] %in% c(unique(all_networks[,2]), go_ids)
    all_networks <- all_networks[ind,]

    target_genes <- intersect(unique(all_networks[,1]), go_ids)

    res_infancy <- mclapply(target_genes, function(x) 
     run_linear_model(x, all_networks, sce[,sce$stage_ids=="Infancy"]),
     mc.cores=10)
    res_infancy <- res_infancy[!sapply(res_infancy, function(x) is.null(x))]
    res_infancy <- do.call(rbind, res_infancy)
    saveRDS(res_infancy, file=paste0("Infancy_Network_Ensheatment_",
    cell_types[i], ".RDS"))
    res_adolescence <- mclapply(target_genes, function(x) 
        run_linear_model(x, all_networks, sce[,sce$stage_ids=="Adolescence"]),
        mc.cores=10)
    res_adolescence <- res_adolescence[!sapply(res_adolescence, 
        function(x) is.null(x))]
    res_adolescence <- do.call(rbind, res_adolescence)
    saveRDS(res_adolescence, file=paste0("Adolescence_Network_Ensheatment_",
    cell_types[i], ".RDS"))
}


