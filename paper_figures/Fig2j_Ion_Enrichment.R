library(readxl)
library(scater)
library(ggplot2)
library(viridis)
library(fgsea)
library(ComplexHeatmap)
library(ggsci)

paths <- "snRNAseq/processed_data/"

sce <- readRDS(paste0(paths,
    "/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))

neuro_trans <- read_xlsx("annotation/Ion_transporters.xlsx")
neuro_trans$Ca <- neuro_trans$Ion %in% c("Ca")
neuro_trans$K <- neuro_trans$Ion %in% c("K")
neuro_trans$Na <- neuro_trans$Ion %in% c("Na")
neuro_trans$NT <- neuro_trans$Ion %in% c("NT")
neuro_trans$comb <-  neuro_trans$Ion %in% c("Na K", "Ca Na K")
neuro_trans <- as.data.frame(neuro_trans)

down_up <- c(`0`="up" , `1`="down", `2`="down", `3`="up", `4`="down", `5`="interup", 
             `6`="up", `7`="up", `8`="interdown", `9`="down", `10`="down", `11`="interup", 
             `12`="interdown", `13`="interup")

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
gene_trends$up_down <- down_up[match(gene_trends$gene_trend,names(down_up))]
gene_trends <- split(gene_trends, gene_trends$major_clust)

universe_cell_type <- lapply(gene_trends, function(x) x$gene_name)

types_neuros <- c("Ca", "Na", "K", "NT")

universe <- rowData(sce)$index

enrichment <- function(genes_neuros, genes, universe){
    
    genes <- split(genes$gene_name, genes$up_down)
    res <- lapply(genes, function(x) fora(genes_neuros, x, universe))
    for(i in 1:length(res)){
        res[[i]]$up_down <- names(res)[i]
    }
    res <- do.call(rbind, res)
    return(res)
    
}

genes_neuros <- lapply(types_neuros, function(x) 
    neuro_trans$`Gene symbol`[neuro_trans[,x]])
names(genes_neuros) <- types_neuros
genes_neuros <- lapply(genes_neuros, function(x)
    intersect(x, universe))
    
res_enrich <- lapply(1:length(gene_trends), function(x)
    enrichment(genes_neuros, gene_trends[[x]], universe_cell_type[[x]]))
names(res_enrich) <- names(gene_trends)
for(i in 1:length(res_enrich)){
    
    res_enrich[[i]]$cell_type <- names(res_enrich)[i]
    res_enrich[[i]]$padj <- p.adjust(res_enrich[[i]]$pval, method="fdr")
    
}
res_enrich <- do.call(rbind, res_enrich)
res_enrich$padj[res_enrich$padj>0.05] <- NA

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

res_enrich$cell_type_up_down <- paste0(res_enrich$cell_type, "_", res_enrich$up_down)
res_enrich$cell_type_up_down <- factor(res_enrich$cell_type_up_down, levels=rev(c("L2-3_CUX2_interup",
    "L4_RORB_up", "L5-6_TLE4_up", "L5-6_TLE4_interup", "ID2_up", "SST_up", "PV_SCUBE3_up", "PV_SCUBE3_interup")))
 
g1 <- ggplot(res_enrich[!is.na(res_enrich$padj),], 
       aes(x=pathway, y=cell_type_up_down, fill=padj)) + 
    geom_point(aes(size=size), pch=21, colour="black") +
    theme_light() +  scale_fill_material("grey")

ggsave(g1, file="paper_figures/Fig2j_Ion_Enrichment.svg", height=7, width=4)

lt <- res_enrich[!is.na(res_enrich$padj)&res_enrich$pathway=="K"&
            res_enrich$up_down=="up",]$overlapGene
names(lt) <- res_enrich[!is.na(res_enrich$padj)&res_enrich$pathway=="K"&
          res_enrich$up_down=="up",]$cell_type
m = make_comb_mat(list_to_matrix(lt))
UpSet(m)
