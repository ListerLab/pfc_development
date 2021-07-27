###############################################################################
#                                                                             #
#                                                                             #
# Apply non-negative matrix factorization                                     #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(NMF)
library(reticulate)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

np <- import("numpy")
sp <- import("scipy")
sklearn <- import("sklearn.decomposition")

source("snATACseq/R/functions_nmf.R")

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

promoters <- promoters(txdb, upstream=2000, downstream=200)

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")

peaks <- lapply(peaks, function(x) {
    hits <- findOverlaps(x, promoters)
    ind <- setdiff(1:length(x), unique(hits@from))
    x[ind]
})
peaks <- unlist(as(peaks, "GRangesList"))

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

# 3. apply nmf

model = sklearn$NMF(n_components=43L, init='random', 
  random_state=0L, verbose="True", max_iter=2000L)
W = model$fit_transform(sam_atac)
H = model$components_

normH = t(apply(H, 1, function(x) x/sum(x)))
normW = apply(W, 2, function(x) x/sum(x))

# 4. find distinct CREs and genes

classes_new <- classes_new [!sapply(classes_new , function(x) is.na(x[1]))]
classes_new <- sapply(classes_new, function(x) getmode(x))

H_class <- def_cell_class(normH)

rna_genes_pseudo_bulk <- readRDS("snATACseq/processed_data/rna_pseudo_bulk.RDS")
sam_rna <- do.call(rbind, 
      rna_genes_pseudo_bulk[match(rownames(sam_atac_new_mat), 
                                  names(rna_genes_pseudo_bulk))])

saveRDS(H_class, "snATACseq/processed_data/nmf_sam_H_class.RDS")

peaks_nmf <- rownames(sam_atac)[o_ind]
sam_atac_new <- lapply(1:length(genes_sam), function(i){
    peaks_tmp <- intersect(peaks_sam[[i]], peaks_nmf)
    if(length(peaks_tmp)==0){
        return(NA)
    } else {
        return(colSums(frags_in_peaks[[genes_sam[i]]][peaks_tmp,, drop=F]))
    }
})
names(sam_atac_new) <- genes_sam
sam_atac_new <- sam_atac_new[!sapply(sam_atac_new, function(x) is.na(x[1]))]
sam_atac_new_mat <- do.call(rbind, sam_atac_new)

classes_new <- lapply(1:length(genes_sam), function(i){
    peaks_tmp <- which(peaks_nmf %in% peaks_sam[[i]])
    if(length(peaks_tmp)==0){
        return(NA)
    } else {
        return(class0[peaks_tmp])
    }
})

saveRDS(list(classes_new=classes_new, sam_atac_new=sam_atac_new_mat,
             rna_seq=sam_rna),
        file="snATACseq/processed_data/nmf_sam_atac_rna.RDS")

saveRDS(list(class0=class0, peaks_nmf=peaks_nmf),
        file="snATACseq/processed_data/peaks_nmf.RDS")

peaks_nmf <- readRDS("snATACseq/processed_data/peaks_nmf.RDS")
peaks_nmf <- split(peaks_nmf$peaks_nmf, peaks_nmf$class0)
peaks_nmf <- lapply(peaks_nmf, function(x)
    intersect(x, names(peaks)))
saveRDS(peaks_nmf, file="snATACseq/processed_data/peaks_nmf_distinct.RDS")

o_featureScore_kim <- cal_featureScore_kim(normW)
o_region_class <- def_region_class(normW)

o_ind <- intersect(o_featureScore_kim[["selt_fs_idx"]],
                   o_region_class[["selt_med_idx"]])
class0 = o_region_class[['class0']][o_ind]

peaks_max_mod <- split(rownames(sam_atac), o_region_class[['class0']])
peaks_max_mod <- lapply(peaks_max_mod, function(x)
    intersect(x, names(peaks)))
saveRDS(peaks_max_mod,
        file="snATACseq/processed_data/peaks_max_mod.RDS")


