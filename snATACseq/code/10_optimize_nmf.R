###############################################################################
#                                                                             #
#                                                                             #
# Optimize NMF                                                                #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(reticulate)
library(parallel)

np <- import("numpy")
sp <- import("scipy")
sklearn <- import("sklearn.decomposition")

source("snATACseq/R/functions_nmf.R")

# 1. read in data

frags_in_peaks <- readRDS("snATACseq/processed_data/frags_in_peaks.RDS")
pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")
colInfo.select <- readRDS(file="snATACseq/processed_data/colInfo.select.RDS")

# 2. prepare data

pearson_corr_sig$Gene <- as.character(pearson_corr_sig$Gene)
genes_sam <- unique(pearson_corr_sig$Gene)
peaks_sam <- split(rownames(pearson_corr_sig), pearson_corr_sig$Gene)
peaks_sam <- lapply(peaks_sam, function(x) sapply(x, function(y) 
    paste0(strsplit(y, ".", fixed=T)[[1]][1:2], collapse=".")))
peaks_sam <- lapply(peaks_sam, function(x) unique(x))

sam_atac <- lapply(1:length(genes_sam), function(i)
    frags_in_peaks[[genes_sam[i]]][peaks_sam[[i]],, drop=F])
sam_atac <- do.call(rbind, sam_atac)

sam_atac <- sam_atac[!duplicated(rownames(sam_atac)),]

# 3. run optimization step 1

run_NMF_stats <- function(rank, sam_atac, max_iter=200L, seed=0L){
    
    
    model = sklearn$NMF(n_components=rank, init='random', 
       random_state=seed, verbose="True", max_iter=max_iter)
    W = model$fit_transform(sam_atac)
    H = model$components_

    normH = t(apply(H, 1, function(x) x/sum(x)))
    normW = apply(W, 2, function(x) x/sum(x))
    entropy<- cal_entropy(normH)[1]
    sparse_H <- cal_sparseness(normH)[1]
    sparse_W <- cal_sparseness(normW)[1]
    RSS_MSS <- cal_rss_mse(W,H,sam_atac)
    
    print("run finished")
    
    return(c(`rank`= rank, `entropy`=entropy, `sparse_H`=sparse_H,
             `sparse_W`=sparse_W, `rss`=RSS_MSS[1], `mss`=RSS_MSS[2]))
}

ranks <- seq(30L, 75L, 5L)

res <- list()

for(i in 1:length(ranks)){
    
    res[[i]] <- run_NMF_stats(ranks[i], sam_atac) 
    
}

# 4. run optimization step 2

ranks <- seq(11L, 17L, 1L)
set.seed(10)
seeds <- sample(1L:10000L, 10L)

res <- list()

for(i in 1:length(ranks)){
    
    
    res[[i]] <- lapply(seeds, function(x)
        run_NMF_stats(ranks[i], sam_atac, seed=x)) 
    
}


res1 <- lapply(res, function(x) do.call(rbind, x))
res1 <- lapply(res1, function(x) colMeans(x))

do.call(rbind, res1)
