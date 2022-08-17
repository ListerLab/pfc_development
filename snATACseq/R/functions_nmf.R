cal_sparseness <- function(X){
    
    absSum = sum(abs(X))
    n = dim(X)[1]*dim(X)[2]
    squareSum = sum(X^2)
    numerator = sqrt(n) - (absSum / sqrt(squareSum))
    denominator = sqrt(n) - 1
    sparseness = numerator / denominator
    return(sparseness)
}


cal_rss_mse <- function (W, H, V){
    
    residualSquare = ((W%*%H) - V)^2
    rss = sum(residualSquare)
    mse = mean(residualSquare)
    out = c(rss, mse)
    return(out)
    
}

cal_entropy <- function(normH){
    
    k = dim(normH)[1]
    n = dim(normH)[2]
    e_list = list()
    
    for(i in 1:k){
        
        h_ij_list = list()
        
        for(j in 1:n){
            
            h_ij = normH[i,j]
            if (h_ij != 0){
                tmp = h_ij * log2(h_ij)
            } else{
                tmp = 0
            }
            h_ij_list[[j]] <- tmp
        }
        
        h_ij_sum = -sum(unlist(h_ij_list))
        e_list[[i]] <- (h_ij_sum)
    }
    
    e_list <- unlist(e_list)
    normInfoGain = 1 - sum(e_list) / (n * log2(k))
    return(c('normInfoGain'=normInfoGain, 'e_list'=e_list))
}


def_cell_class<- function(normH, contribute=0.3){
    
    index = 1:dim(normH)[2]
    colmax = apply(normH, 2, function(x) max(x))
    colsum = colSums(normH)
    contributes = colmax / colsum
    class0 = apply(normH, 2, function(x) which.max(x))
    unclass = which(contributes < contribute)
    unclassN = length(unclass)
    unclassPct = unclassN/dim(normH)[2]
    class1 = class0
    class1[unclass] = 100
    return(list('class0'=class0, 'class1'=class1, 'contributes'=contributes, 
                'unclass'=unclass, 'pct'=unclassPct))
}

cal_featureScore_kim <- function(W){
    
    k = dim(W)[2]
    m = dim(W)[1]
    fs_list = list()
    
    for(i in 1:m){
        
        rowsum <- sum(W[i,])
        p_iq_x_list = list()
        for(q in 1:k){
            p_iq = W[i,q] / rowsum
            if (p_iq != 0){
                tmp = p_iq * log2(p_iq)
            } else{
                tmp = 0
            }
            p_iq_x_list[[q]] <- tmp
        }
        fs = 1 + 1/log2(k) * sum(unlist(p_iq_x_list))
        fs_list[[i]] <- fs
    }
    
    fs_list <- unlist(fs_list)
    med = median(fs_list)
    mad = median(abs(fs_list - med))
    fs_cutoff = med + 3 * mad
    selt_fs_idx = which(fs_list >= fs_cutoff)
    return(list('fs'=fs_list, 'med'=med, 'mad'=mad, 
                'fs_cutoff'= fs_cutoff, 'selt_fs_idx'= selt_fs_idx))
}


def_region_class <- function(normW){
    
    rowmax <- apply(normW, 1, function(x) max(x))
    rowsum <- rowSums(normW)
    contributes <- rowmax/rowsum
    class0 <- apply(normW, 1, function(x) which.max(x))
    med_contribute <- median(normW)
    selt_med_idx <- which(contributes > med_contribute)
    return(list('class0'=class0, 'contributes'=contributes, 
                'med_contribute'=med_contribute, 'selt_med_idx'=selt_med_idx))
}

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

go_term_analysis <- function(genes, genes_background){
    
    
    genes <-  AnnotationDbi::select(org.Hs.eg.db, keys=genes, columns="ENTREZID", 
                                    keytype="ENSEMBL")[,2]
    genes <- genes[!is.na(genes)]
    genes_background <- AnnotationDbi::select(org.Hs.eg.db, keys=genes_background, 
      columns="ENTREZID", keytype="ENSEMBL")[,2]
    genes_background <- genes_background[!is.na(genes_background)]
    
    go_bp <- enrichGO(gene= genes,
                      universe = genes_background,
                      OrgDb = org.Hs.eg.db,
                      ont  = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05,
                      readable = TRUE)
    return(go_bp)
    
}
