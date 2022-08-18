library(readxl)
library(scater)
library(ggplot2)
library(viridis)
library(fgsea)
library(ComplexHeatmap)
library(ggsci)

res_enrich <- readRDS("snRNAseq/processed_data/Ion_transport.RDS")

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

res_enrich$cell_type_up_down <- paste0(res_enrich$cell_type, "_", res_enrich$up_down)
res_enrich$cell_type_up_down <- factor(res_enrich$cell_type_up_down, levels=rev(c("L2-3_CUX2_interup",
    "L4_RORB_up", "L5-6_TLE4_up", "L5-6_TLE4_interup", "ID2_up", "SST_up", "PV_SCUBE3_up", "PV_SCUBE3_interup")))
 
g1 <- ggplot(res_enrich[!is.na(res_enrich$padj),], 
       aes(x=pathway, y=cell_type_up_down, fill=padj)) + 
    geom_point(aes(size=size), pch=21, colour="black") +
    theme_light() +  scale_fill_material("grey")

ggsave(g1, file="paper_figures/Fig2_Ion_Enrichment.svg", height=7, width=4)

lt <- res_enrich[!is.na(res_enrich$padj)&res_enrich$pathway=="K"&
            res_enrich$up_down=="up",]$overlapGene
names(lt) <- res_enrich[!is.na(res_enrich$padj)&res_enrich$pathway=="K"&
          res_enrich$up_down=="up",]$cell_type
m = make_comb_mat(list_to_matrix(lt))
UpSet(m)
