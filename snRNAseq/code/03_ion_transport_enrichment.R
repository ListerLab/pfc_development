library(readxl)
library(scater)
library(ggplot2)
library(viridis)
library(fgsea)
library(ComplexHeatmap)
library(ggsci)

source("snRNAseq/R/functions_disease_enrichment.R")

paths <- "snRNAseq/processed_data/"

sce <- readRDS(paste0(paths,
   "2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"))

neuro_trans <- read_xlsx("annotation/Ion_transporters.xlsx")
neuro_trans$DEG <- neuro_trans$`Gene symbol` %in% gene_trends$gene_name

write.table(neuro_trans, "annotation/Ion_transporters_updated.tsv",
          row.names=F, quote=F, col.names=T, sep="\t")
neuro_trans <- neuro_trans[!is.na(neuro_trans$`Channel Type`),]
neuro_trans$`Channel Ion` <- paste0(neuro_trans$`Channel Type`, "..",
  neuro_trans$Ion)

neuro_trans$Ca <- neuro_trans$Ion %in% c("Ca")
neuro_trans$K <- neuro_trans$Ion %in% c("K")
neuro_trans$Na <- neuro_trans$Ion %in% c("Na")
neuro_trans$NT <- neuro_trans$Ion %in% c("NT")
neuro_trans$comb <-  neuro_trans$Ion %in% c("Na K", "Ca Na K")
neuro_trans <- as.data.frame(neuro_trans)

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
gene_trends <- split(gene_trends, gene_trends$major_clust)

cell_types <- names(input)

universe <- find_back_genes(sce, 0.05, input, cell_types)

universe_cell_type <- lapply(gene_trends, function(x) x$gene_name)

types_neuros <- unique(neuro_trans$`Channel Ion`)

enrichment <- function(genes_neuros, genes, universe){
    
    genes <- split(genes$gene_name, genes$major_trend)
    res <- lapply(genes, function(x) fora(genes_neuros, x, universe))
    for(i in 1:length(res)){
        res[[i]]$up_down <- names(res)[i]
    }
    res <- do.call(rbind, res)
    return(res)
    
}

genes_neuros <- lapply(types_neuros, function(x) 
    neuro_trans$`Gene symbol`[neuro_trans$`Channel Ion`==x])
names(genes_neuros) <- types_neuros
genes_neuros <- lapply(genes_neuros, function(x)
    intersect(x, unlist(universe)))

res_enrich <- lapply(1:length(gene_trends), function(x)
    enrichment(genes_neuros, gene_trends[[x]], universe_cell_type[[x]]))
names(res_enrich) <- names(gene_trends)
for(i in 1:length(res_enrich)){
    
    res_enrich[[i]]$cell_type <- names(res_enrich)[i]
    res_enrich[[i]]$padj <- p.adjust(res_enrich[[i]]$pval, method="fdr")
    
}
res_enrich <- do.call(rbind, res_enrich)
res_enrich$padj[res_enrich$padj>0.05] <- NA

saveRDS(res_enrich, file="snRNAseq/processed_data/ion_transport_enrich.rds")



