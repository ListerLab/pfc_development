library(readxl)
library(scater)
library(ggplot2)
library(viridis)
library(fgsea)
library(ComplexHeatmap)
library(ggsci)

res_enrich <- readRDS("snRNAseq/processed_data/ion_transport_enrich.rds")
res_enrich <- res_enrich[!is.na(res_enrich$padj),]

res_enrich$cell_type_up_down <- paste0(res_enrich$cell_type, "_", 
  res_enrich$up_down)
res_enrich$cell_type_up_down <- factor(res_enrich$cell_type_up_down, 
    levels=rev(c("L2-3_CUX2_interup", "L5-6_THEMIS_interup",
   "L5-6_TLE4_interup", "ID2_up", "PV_SCUBE3_interup")))
 
g1 <- ggplot(res_enrich[!is.na(res_enrich$padj),], 
       aes(x=pathway, y=cell_type_up_down, fill=padj)) + 
    geom_point(aes(size=size), pch=21, colour="black") +
    theme_light() +  scale_fill_material("grey") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.4, hjust=0.4))
  

ggsave(g1, file="main_figures/ion_enrichment_scrna.svg", height=7, width=4)


