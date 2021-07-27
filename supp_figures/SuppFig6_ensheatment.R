library(qusage)
library(rtracklayer)
library(SingleCellExperiment)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(reshape2)
library(cba)
library(corrplot)
library(cowplot)
library(viridis)
library(ggplot2)
library(dplyr)
library(jaccard)
library(cowplot)

# 1. read in data

cis <- readRDS(
    "processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
cis <- lapply(cis, function(x) x[!is.na(x$CRE)]) # only CRE elements
cis <- as(cis, "GRangesList")

cis <- cis[!names(cis)=="Vas"]

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sample_by_stages <- split(atac_samples$Sample, atac_samples$Stage)

txdb <- makeTxDbFromGFF("annotation/genes.gtf", format = "gtf")
seqlevelsStyle(txdb) <- "UCSC"

paths <- "../shiny_app_celltypes/"
go_terms <- read.gmt("annotation/hsapiens.GO:BP.name.gmt")
sce <- readRDS(paste0(paths, 
   "2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))

gene_trends <- read.csv(paste0(paths, "../scRNAseq/gene_cluster_ids.csv"))
gene_trends$gene_name <- rowData(sce)$gene_ids[match(gene_trends$gene_name, 
                                                     rownames(sce))]

# 2. remove promoters

promoters <- promoters(txdb, upstream=2000, downstream=200)

cis <- lapply(cis, function(x) {
    hits <- findOverlaps(x, promoters)
    ind <- setdiff(1:length(x), unique(hits@from))
    x[ind]
})


gene_trends <- read.csv(paste0(paths, "../scRNAseq/gene_cluster_ids.csv"))
gene_trends$gene_name <- rowData(sce)$gene_ids[match(gene_trends$gene_name, 
                                                     rownames(sce))]

# 4. get regions per trajectory in ensheatment pathway

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
                       "OPC", "OPC",
                       "Micro", "Micro"
), ncol=2, byrow=T)

gene_trends <- lapply(names(cis), function(x) 
    gene_trends[gene_trends$major_clust %in% cell_types[cell_types[,2]==x,1],])
names(gene_trends) <- names(cis)

trajectories_cell_types <- lapply(gene_trends, function(x) x$gene_name)

ind <- which(names(go_terms) %in% "GO:0007272")
go_id <- go_terms[[ind]]

go_id <- rowData(sce)$gene_ids[match(go_id, rowData(sce)$index)]
go_id <- go_id[!is.na(go_id)]

subset_peaks <- function(peaks_tmp, trajectories_cell_type, go_id){
    
    genes_tmp <- intersect(go_id, trajectories_cell_type)
    
    ind <- unique(unlist(sapply(genes_tmp, function(x) grep(x, peaks_tmp$CRE))))
    peaks_tmp <- peaks_tmp[ind]
    peaks_tmp
}

cis_go_id <- mapply(function(X,Y) subset_peaks(X, Y, go_id),  X=cis, 
                    Y=trajectories_cell_types)


cis_go_id_all <- unlist(as(cis_go_id, "GRangesList"), recursive = TRUE)
cis_go_id_all <- reduce(cis_go_id_all)

cis_go_id <- cis_go_id[!lengths(cis_go_id)<=10]

# 4. find fpkm for each stage 

collect_fpkm_fun <- function(X, sample_by_stages) {
    aa<- X@elementMetadata
    aa<- aa[,grepl("_RL", colnames(aa))]
    ab <- as.matrix(aa)
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
    rownames(ac) <- names(X)
    return(ac)
}

collect_fpkm <- lapply(cis_go_id, function(x) 
    collect_fpkm_fun(x, sample_by_stages))

collect_peaks <- lapply(cis_go_id, function(x) 
    x@elementMetadata[,stages])

# 5. plot corrplots

mat <- list()

for(i in 1:length(collect_fpkm)){
    
    aa <- ncol(collect_fpkm[[i]]) 
    mat[[i]] <- sapply(1:aa, function(rr) sapply(1:aa, function(cc)
        cor(collect_fpkm[[i]][,rr], collect_fpkm[[i]][,cc], 
            method="spearman")))
    colnames(mat[[i]]) <- rownames(mat[[i]]) <- stages
    mat[[i]][which(is.na(mat[[i]]), arr.ind=T)] <- 0
    
    aa <- ncol(collect_peaks[[i]])
    tmp <- sapply(1:aa, function(rr) sapply(1:aa, function(cc)
        jaccard(collect_peaks[[i]][,rr], collect_peaks[[i]][,cc])))
    diag(tmp) <- 0
    colnames(tmp) <- rownames(tmp) <- stages
    
    col2 <- colorRampPalette(c("#67001F", "#D6604D",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))
    mat_tmp <- melt(mat[[i]])
    colnames(mat_tmp) <- c("X", "Y", "Correlation")
    tmp_tmp <- melt(tmp)
    colnames(tmp_tmp) <- c("X", "Y", "Jaccard")
    
    g1 <- ggplot(mat_tmp, aes(x=X, y=Y, fill=Correlation)) + geom_tile() +
        theme_classic() + scale_fill_gradientn(colours=col2(20), limits=c(-0.5,1)) +
        theme(axis.title=element_blank(), 
              panel.background=element_rect(fill="white", colour="white"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ggtitle(names(collect_fpkm)[i])
    
    gg <- plot_grid(plotlist=list(g1,g2))
    ggsave(g1, file=paste0("supp_figures/SuppFig6_ensheatment/", 
        names(collect_fpkm)[i], ".svg"), height=4.7, width=10)
    
}