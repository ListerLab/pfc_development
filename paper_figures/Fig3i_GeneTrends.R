# ATACseq trends

library(rtracklayer)
library(readxl)
library(reshape2)
library(dplyr)
library(parallel)
library(SingleCellExperiment)
library(ggplot2)
library(directlabels)

atac_meta <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_meta <- atac_meta[atac_meta$Used=="Yes",]

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm=T))}

paths <- "snRNAseq/processed_data/"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
sce <- readRDS(paste0(paths, 
    "2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_stage_CRE_gr.Rds")
cis <- lapply(peaks, function(x) x[!is.na(x$CRE)])
cis <- as(cis, "GRangesList")

cis <- cis[!names(cis)=="Vas"]

collect_fpkm_fun <- function(X) {
    aa <- X@elementMetadata
    aa <- aa[,grepl("_RL", colnames(aa))]
    aa <- as.matrix(aa)
    aa <- t(apply(aa, 1, function(x) range01(x)))
    aa <- as.data.frame(aa)
    aa$peaks <- names(X)
    aa <- melt(aa, id="peaks")
    return(aa)
}

gene_trends$gene_name <- rowData(sce)$gene_ids[match(gene_trends$gene_name, 
    rownames(sce))]
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

gene_trends <- lapply(names(cis),  function(x)
    gene_trends[gene_trends$major_clust %in%
                    cell_types[cell_types[,2]==x,1],])
names(gene_trends) <- names(cis)

for(i in 1:length(cis)){
    
    ind <- sapply(cis[[i]]$CRE, function(x) strsplit(x, "|", fixed=T)[[1]])
    ind <- sapply(ind, function(x) any(x %in% gene_trends[[i]]$gene_name))
    cis[[i]] <- cis[[i]][ind]
}


collect_fpkm <- lapply(cis, function(x) 
    collect_fpkm_fun(x))

for(i in 1:length(collect_fpkm)){
    
    collect_fpkm[[i]]$cell_type <- names(collect_fpkm)[i]
    collect_fpkm[[i]]$variable <- gsub(paste0(names(collect_fpkm)[i], "_"), "",
          collect_fpkm[[i]]$variable)
    collect_fpkm[[i]]$stages <-atac_meta$Stage[
        match(collect_fpkm[[i]]$variable, atac_meta$Sample)]
}

#get associated genes and gene trends

attach_gene_trends <- function(tmp_fpkm, cis_tmp, gene_trends_tmp){
    
    cre_genes <- cis_tmp$CRE[match(tmp_fpkm$peaks, 
            names(cis_tmp))]
    all_trends <- lapply(cre_genes, function(a)
        gene_trends_tmp[grep(a, gene_trends_tmp[,3]),4])
    ind <- lapply(1:length(all_trends), function(x)
        rep(x, length(all_trends[[x]])))
    ind <- unlist(ind)
    tmp_fpkm <- tmp_fpkm[ind,]
    tmp_fpkm$all_trends <- unlist(all_trends)
    return(tmp_fpkm)
    
}

collect_fpkm_new <- mclapply(1:length(cis), function(i) 
    attach_gene_trends(collect_fpkm[[i]], 
        cis[[i]], gene_trends[[i]]), mc.cores=10)


collect_fpkm_new <- do.call(rbind, collect_fpkm_new)
collect_fpkm_new <- collect_fpkm_new[!is.null(unlist
    (collect_fpkm_new$all_trends)),]

new_trend <- matrix(c(0, "up",
                      1, "down",
                      2, "down",
                      3, "up",
                      4, "down",
                      5, "inter up",
                      6, "up",
                      7, "up",
                      8, "inter down",
                      9, "down",
                      10, "down",
                      11, "inter up",
                      12, "inter down",
                      13, "inter up"), byrow=TRUE, ncol=2)

collect_fpkm_new$trend_new <- new_trend[
    as.numeric(collect_fpkm_new$all_trends)+1,2]

saveRDS(collect_fpkm_new, file="snATACseq/processed_data/collect_fpkm_new.RDS")
collect_fpkm_new <- readRDS("snATACseq/processed_data/collect_fpkm_new.RDS")

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

collect_fpkm_new$stages <- factor(collect_fpkm_new$stages, levels=stages)

cols <- c(down="#02577F", `inter down`="#D9D9D7",
          `inter up`="#FFD433", up="#F3704A")

collect_fpkm_new$stages_num <- as.numeric(collect_fpkm_new$stages)

dlabs <- collect_fpkm_new %>%
    group_by(all_trends) %>% 
    arrange(stages_num) %>% 
    filter(stages_num %in% c(first(stages_num), last(stages_num)))

g1 <- ggplot(collect_fpkm_new, aes(x=stages_num, y=value, 
    group=all_trends, col=trend_new)) +
    stat_summary(fun.y=mean, geom="line", size = 1) +
    theme_classic() + 
    stat_summary(data = dlabs, aes(label = all_trends, 
        hjust = ifelse(stages_num == first(stages_num), 1.2, -.2)), 
        fun = mean, geom = "text", color = "black")  +
    theme(axis.title.x = element_blank()) + ylab("Normalized fpkm") +
    scale_color_manual(values=cols) + scale_x_continuous(labels=stages,
        breaks=1:6,expand = expansion(mult = .1))
ggsave(g1, file="paper_figures/Fig3i_GeneTrends.svg", height=4, width=4)
    



    