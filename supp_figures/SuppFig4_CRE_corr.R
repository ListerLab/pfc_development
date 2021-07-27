library(GenomicRanges)
library(corrplot)
library(readxl)
library(jaccard)
library(viridis)
library(svglite)
library(ggplot2)

cis <- readRDS("../processed_data/cell_type_atac_peaks_filtered_anno_dge_stage_motif_CRE_gr_no_long.Rds")
cis <- lapply(cis, function(x) x[!is.na(x$CRE)])
cis <- as(cis, "GRangesList")

all_cis <- unlist(cis)
all_cis <- reduce(all_cis)

atac_samples <- read_xlsx("../annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sample_by_stages <- split(atac_samples$Sample, atac_samples$Stage)

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

overlap_index <- lapply(cis, function(x) 
    findOverlaps(all_cis, x))
overlap_index <- lapply(overlap_index, function(x)
    data.frame(V1=x@to, V2=x@from)
)

collect_access <- mapply(function(X, Y) {
    aa<- X[Y$V1]@elementMetadata[,stages]
    aa <- as.matrix(aa)
    ind <- split(1:nrow(aa),Y$V2)
    ab <- t(sapply(ind, function(xx) colSums(aa[xx, , drop=FALSE])))
    ab <- as.data.frame(ab)
    ab$peaks <- names(ind)
    return(ab)}, X=cis, Y=overlap_index, SIMPLIFY = FALSE)

for(i in 1:length(collect_access)){
    
    ncols <- ncol(collect_access[[i]])
    colnames(collect_access[[i]])[-ncols] <- paste0(
        names(collect_access)[i], "_",
        colnames(collect_access[[i]])[-ncols])
}

all_access <- merge(x=collect_access[[1]], y=collect_access[[2]],
                    by="peaks",  all=TRUE)
all_access[which(is.na(all_access), arr.ind = T)] <- 0

for(i in 3:length(collect_access)){
    
    all_access <- merge(x=all_access, y=collect_access[[i]],
                        by="peaks",all=TRUE)
    all_access[which(is.na(all_access), arr.ind = T)] <- 0
    
}

rownames(all_access) <- all_access$peaks
all_access <- all_access[,-1]

aa <- ncol(all_access)

mat <- sapply(1:aa, function(rr) sapply(1:aa, function(cc)
    jaccard(all_access[,rr], all_access[,cc])))
mat[which(is.na(mat), arr.ind = TRUE)] <- 0
colnames(mat) <- rownames(mat) <- colnames(all_access)
diag(mat) <- 1

corrplot(mat, method="color",  rect.col = "white", tl.col = "black",
         col=c(magma(256), magma(256)), is.corr = FALSE, cl.lim=c(0,1),
         diag=FALSE,tl.cex=0.8, mar=c(0,0,1,0), order="hclust")
```

## Overlap collect fpkms

```{r, fig.height=10, fig.width=10}
stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

overlap_index <- lapply(cis, function(x) 
    findOverlaps(all_cis, x))
overlap_index <- lapply(overlap_index, function(x)
    data.frame(V1=x@to, V2=x@from)
)

collect_fpkm_fun <- function(X, Y, sample_by_stages) {
    aa<- X[Y$V1]@elementMetadata
    aa<- aa[,grepl("_RL", colnames(aa))]
    aa <- as.matrix(aa)
    ind <- split(1:nrow(aa),Y$V2)
    ab <- t(sapply(ind, function(xx) colMeans(aa[xx, , drop=FALSE])))
    ac <- lapply(sample_by_stages, function(xx) {
        ind1 <- grepl(paste0(xx, collapse="|"), colnames(ab))
        if(sum(ind1)==0) {
            rep(0, length=nrow(ab))
        } else {
            rowMeans(ab[, ind1, drop=FALSE])
        }
    })
    ac <- do.call(cbind ,ac)
    ac <- as.data.frame(ac)
    ac$peaks <- names(ind)
    return(ac)
}

collect_fpkm<-   mapply(function(X, Y) 
    collect_fpkm_fun(X, Y, sample_by_stages), X=cis, Y=overlap_index, 
    SIMPLIFY = FALSE)


for(i in 1:length(collect_fpkm)){
    
    ncols <- ncol(collect_fpkm[[i]])
    colnames(collect_fpkm[[i]])[-ncols] <- paste0(
        names(collect_fpkm)[i], "_",
        colnames(collect_fpkm[[i]])[-ncols])
}

all_fpkm <- merge(x=collect_fpkm[[1]], y=collect_fpkm[[2]],
                  by="peaks",  all=TRUE)
all_fpkm[which(is.na(all_fpkm), arr.ind = T)] <- 0

for(i in 3:length(collect_fpkm)){
    
    all_fpkm <- merge(x=all_fpkm, y=collect_fpkm[[i]],
                      by="peaks",all=TRUE)
    all_fpkm[which(is.na(all_fpkm), arr.ind = T)] <- 0
    
}

rownames(all_fpkm) <- all_fpkm$peaks
all_fpkm <- all_fpkm[,-1]

aa <- ncol(all_fpkm)

mat <- sapply(1:aa, function(rr) sapply(1:aa, function(cc)
    cor(all_fpkm[,rr], all_fpkm[,cc], method="spearman")))
mat[which(is.na(mat), arr.ind = TRUE)] <- 0
colnames(mat) <- rownames(mat) <- colnames(all_access)
diag(mat) <- 1

corrplot(mat, method="color",  rect.col = "white", tl.col = "black",
         is.corr = FALSE, cl.lim=c(0,1),
         diag=FALSE, tl.cex=0.8, mar=c(0,0,1,0), order="hclust")

svg("supp_figures/SuppFig4_CREcorr.svg", height=9, width=9)
corrplot(mat, method="color",  rect.col = "white", tl.col = "black",
         is.corr = FALSE, cl.lim=c(0,1),
         diag=FALSE, tl.cex=0.8, mar=c(0,0,1,0), order="hclust")
dev.off()