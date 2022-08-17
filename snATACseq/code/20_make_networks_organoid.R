# Make network

library(scater)
library(parallel)
library(reticulate)
library(readxl)

np <- import("numpy")
sklearn <- import("sklearn.linear_model")

cell_types <- c(L4="L4_RORB", L5_6="L5/6_TLE4", L5_6="L5/6_THEMIS",
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
    pvalue <- pnorm(abs(coef_mean), mean=0, sd=sqrt(coef_var), 
      lower.tail=FALSE)*2
    
    if(sum(pvalue<0.05)>0){
        
        res <- data.frame(tfs=predictors[pvalue<0.05], target=target_gene,
            strength=coef_mean[pvalue<0.05], pvalue=pvalue[pvalue<0.05])
        return(res)
        
    } else {
        
        return(NULL)
    }
}


motifs <- read_xlsx("annotation/motif_annotations.xlsx", sheet=2)
motifs <- motifs[motifs$Database=="Jaspar2018",]
motifs <- sapply(motifs$Motif, function(x) strsplit(x, "_")[[1]][1])
motifs <- toupper(motifs)
motifs <- sapply(motifs, function(x) strsplit(x, "+", fixed=TRUE)[[1]])
motifs <- unlist(motifs)

markers <- read.csv("snRNAseq/processed_data/Organoids_Brain_markers_PN.csv")
markers_up <- rownames(markers)[markers$FDR<0.05 & markers$summary.logFC>0]
markers_down <- rownames(markers)[markers$FDR<0.05 & markers$summary.logFC<0]

tf_up <- markers_up[markers_up %in% motifs]
tf_down <- markers_down[markers_down %in% motifs]

for(i in 1:length(cell_types)){

    all_networks <- readRDS("snATACseq//processed_data/all_networks.rds")
    all_networks <- all_networks[[names(cell_types)[i]]]
    tf_cell_type <- readRDS("snATACseq//processed_data/tf_not_cell_type_rna.RDS")
    tf_cell_type <- tf_cell_type[[cell_types[i]]]
    gc()
    
    all_networks[,2] <- toupper(all_networks[,2])
    tmp_tf <- unique(all_networks[,2])
    tmp_tf <- sapply(tmp_tf, function(x) strsplit(x, "..", fixed=T)[[1]])
    tmp_tf$`EWSR1.FLI1` <- c("EWSR1", "FLI1")

    ind <- mclapply(tmp_tf, function(x) match(x, 
        tf_cell_type), mc.cores=2)
    keep <- names(ind)[sapply(ind, function(x) all(is.na(x)))]
    all_networks <- all_networks[all_networks[,2] %in% keep,]
    tmp <- mclapply(1:dim(all_networks)[1], function(x) {
        a <- strsplit(all_networks[x,2], "..", fixed=T)[[1]]
        cbind(all_networks[x,1], a)
    }, mc.cores=10)
    tmp <- do.call(rbind, tmp)
    all_networks <- tmp

    rm(tmp)
    gc()

    organoids_gene <- read.csv("snRNAseq/processed_data/Organoids_Brain_markers_PN.csv")
    markers_brain <- rownames(organoids_gene)[organoids_gene$summary.logFC<0 &
      organoids_gene$FDR<0.05]
    markers_organoids <- rownames(organoids_gene)[organoids_gene$summary.logFC>0 &
      organoids_gene$FDR<0.05]

    sce <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")
    sce <- sce[, sce$major_clust==cell_types[i]]
    sce <- logNormCounts(sce)

    ind_brain <- all_networks[,1] %in% c(unique(all_networks[,2]), markers_brain)
    all_networks_brain <- all_networks[ind_brain,]
    all_networks_brain <- all_networks_brain[all_networks_brain[,2] %in% tf_down,]
    
    ind_organoids <- all_networks[,1] %in% c(unique(all_networks[,2]), 
      markers_organoids)
    all_networks_organoids <- all_networks[ind_organoids,]
    all_networks_organoids <- all_networks_brain[all_networks_organoids[,2] %in% tf_up,]
    
    target_genes_brain <- intersect(unique(all_networks[,1]), 
      unique(all_networks_brain[,1]))
      target_genes_organoids <- intersect(unique(all_networks[,1]), 
      unique(all_networks_organoids[,1]))

    res_neonatal_organoids <- mclapply(target_genes_organoids, function(x) 
     run_linear_model(x, all_networks_organoids, sce[,sce$stage_ids=="Neonatal"]),
     mc.cores=10)
    res_neonatal_brain <- mclapply(target_genes_brain, function(x) 
     run_linear_model(x, all_networks_brain, sce[,sce$stage_ids=="Neonatal"]),
     mc.cores=10)
    res_neonatal_organoids <- res_neonatal_organoids[
      !sapply(res_neonatal_organoids, function(x) is.null(x))]
    res_neonatal_organoids<- do.call(rbind, res_neonatal_organoids)
    res_neonatal_brain <- res_neonatal_brain[
      !sapply(res_neonatal_brain, function(x) is.null(x))]
    res_neonatal_brains<- do.call(rbind, res_neonatal_brain)
    saveRDS(res_neonatal_organoids, file=paste0("snATACseq//processed_data/",
    "neonatal_network_organoids_", gsub("/", "_", cell_types[i]), "_degs.rds"))
    saveRDS(res_neonatal_brains, file=paste0("snATACseq//processed_data/",
    "neonatal_network_brain_", gsub("/", "_", cell_types[i]), "_degs.rds"))
}


