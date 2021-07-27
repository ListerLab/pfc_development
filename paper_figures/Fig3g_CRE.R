# CRE plot

library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ComplexHeatmap)
library(gridExtra)
library(viridis)
library(org.Hs.eg.db)
library(cba)
library(clusterProfiler)
library(SingleCellExperiment)
library(AnnotationDb)

source("snATACseq/R/functions_nmf.R")

plot_list_old <- readRDS("snATACseq/processed_data/plot_atac_rna_corr.RDS")
info_type <- plot_list_old$info$type
info_stages_rna = plot_list_old$info$stages_rna

H_class <- readRDS("snATACseq/processed_data/nmf_sam_H_class.RDS")
tab_cluster_H<- table(info_type, H_class$class0)
tab_stage_H <- table(info_stages_rna, H_class$class0)

# Pie charts

cluster_of_interest <- c(9, 40, 32, 18, 41, 16, 23, 33, 24)

plot_pie <- function(x){
    
    dat_H <- data.frame(cluster=rownames(tab_cluster_H),
                        value=tab_cluster_H[,x])
    dat_stage <- data.frame(stage=rownames(tab_stage_H),
                        value=tab_stage_H[,x])
    
    col_stages_rna <- viridis(6, direction=1)
    names(col_stages_rna) <- c("Fetal", "Neonatal", "Infancy", "Childhood", 
            "Adolescence", "Adult")
    col_types <- c(Astro= '#ffc857', `IN dev`='#c6a785', 
       `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
        Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN dev`='#e2cdd0',
       `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')
    
    g1 <- ggplot(dat_stage, aes(x="", y=value, fill=stage))  + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=col_stages_rna) + theme_void() +
        guides(fill=FALSE)
    ggsave(g1, file=paste0("paper_figures/pie_charts/", x, "_pie_chart_stage.svg"),
           height=2, width=2)
    g2 <- ggplot(dat_H, aes(x="", y=value, fill=cluster))  + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=col_types) + theme_void() +
        guides(fill=FALSE)
    ggsave(g2, file=paste0("paper_figures/pie_charts/", x, "_pie_chart_cluster.svg"),
           height=2, width=2)
    
}

lapply(cluster_of_interest, function(x) plot_pie(x))

# Heatmap

plot_list <- readRDS("snATACseq/processed_data/nmf_sam_atac_rna.RDS")
plot_list_atac <- plot_list$sam_atac_new
plot_list_rna <- plot_list$rna_seq

ind_classes <- split(1:nrow(plot_list_rna), plot_list$classes_new)
plot_list_rna_classes <- t(sapply(ind_classes, function(x)
    colSums(plot_list_rna[x,, drop=FALSE])))
rownames(plot_list_rna_classes) <- names(ind_classes)
hc <- hclust(dist(plot_list_rna_classes))
order_row <- order.optimal(dist(plot_list_rna_classes), hc$merge)
classes_order <- rownames(plot_list_rna_classes)[order_row$order]
plot_list$classes_new <- factor(plot_list$classes_new, levels=classes_order)

my_hclust_row <- order(plot_list$classes_new)
dat_rna <- t(apply(plot_list_rna, 1, function(x) scale(x)))
dat_atac <- t(apply(plot_list_atac,1, function(x) scale(x)))

hc <- hclust(dist(t(dat_rna[my_hclust_row,])))
order_col <- order.optimal(dist(t(dat_rna[my_hclust_row,])), hc$merge)
my_hclust_cell <- order_col$order


dat_rna <- dat_rna[my_hclust_row, my_hclust_cell]
dat_atac <- dat_atac[my_hclust_row, my_hclust_cell]

col_stages_rna <- viridis(6, direction=1)
names(col_stages_rna) <- c("Fetal", "Neonatal", "Infancy", "Childhood", 
                           "Adolescence", "Adult")
col_types <- c(Astro= '#ffc857', `IN dev`='#c6a785', 
      `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
       Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN dev`='#e2cdd0',
      `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')

ha_left = HeatmapAnnotation(cell = info_type[my_hclust_cell],
                            stages = factor(info_stages_rna[my_hclust_cell],
                                            level=names(col_stages_rna)),
                            annotation_name_side = "left",
                            col = list(stages = col_stages_rna,
                                       cell=col_types))
ha_right = HeatmapAnnotation(cell = info_type[my_hclust_cell],
                             stages = factor(info_stages_rna[my_hclust_cell],
                                             level=names(col_stages_rna)),
                             annotation_name_side = "right",
                             col = list(stages = col_stages_rna,
                                        cell=col_types))

cluster_col<- c('#155145', '#d3d3d3', '#FE1886', '#BDC90F', '#4903FB', '#26FFA8',
    '#A184C0', '#910716', '#06E802', '#0798F1', '#68B167', '#D433FC',
    '#D96F4E', '#FDF853', '#602683', '#72D6F5', '#577A05', '#FDB2AB',
    '#3274A3', '#A8EE91', '#78FA2B', '#7254EB', '#FC3314', '#1C2DCA',
    '#970FCC')
names(cluster_col) <- levels(plot_list$classes_new)

png("paper_figures/Fig3g_CRE.png",  width=10, height=12, unit="cm",
    res=1200)
f1 = circlize::colorRamp2(seq(-4, 4, length = 3), c("blue", "#EEEEEE", "red"))
rowAnnotation(cluster = as.factor(plot_list$classes_new[my_hclust_row]),
              show_legend=FALSE, 
              col=list(cluster=cluster_col)) +
Heatmap(dat_rna,  top_annotation = ha_left, col=f1,
        cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = F,
        show_column_names = F, name = "scRNAseq", show_heatmap_legend = T) +
Heatmap(dat_atac, top_annotation = ha_right, col=f1,
        cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = F,
        show_column_names = F, name="scATACseq", show_heatmap_legend = T)
dev.off()

# Pathway Analysis

sce <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")

all_genes <- lapply(ind_classes, function(x) rownames(plot_list_atac)[x])
all_genes <- all_genes[lengths(all_genes)>10]

res <- lapply(all_genes, function(x) go_term_analysis(x, rowData(sce)$gene_ids))

for (i in 1:length(res)){
    
    res[[i]] <- res[[i]]@result
    res[[i]] <- as.data.frame(res[[i]])
    res[[i]]$cluster <- names(res)[i]
    
}

res <- do.call(rbind, res)

write.table(res, file="snATACseq/processed_data/Go_terms_nmf_cluster.tsv", 
            sep="\t", quote=F)
