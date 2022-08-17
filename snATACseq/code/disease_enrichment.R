###############################################################################
#                                                                             #
#                                                                             #
# Disease enrichment plot                                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(data.table)
library(disgenet2r)
library(scater)
library(fgsea)
library(ggplot2)
library(viridis)
library(cowplot)
library(cba)
library(dplyr)
library(scran)
library(reshape2)
disgenet_api_key <- "b1287e0476b5bc0593e130d55ce1479cc7b8e84b"

# 1. read in data

sce <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")

source("snRNAseq/R/functions_disease_enrichment.R")

gene_trends <- read.csv("snRNAseq/processed_data/gene_cluster_ids.csv")

# 2. pick genes

diseases <- fread("annotation/disgenet_all_diseases.tsv", header=T, 
  data.table=F)

use_diseases_sub_name <- c("Alcoholic Intoxication, Chronic", 
  "Alzheimer's Disease", "Amnesia", "Attention deficit hyperactivity disorder", 
  "Autism Spectrum Disorders", "Bipolar Disorder", "Dementia", 
  "Learning Disabilities", "Huntington Disease", "Impaired cognition", 
  "Major Depressive Disorder","Manic", "Mild cognitive disorder", "Narcolepsy", 
  "Nicotine Dependence", "Schizoaffective Disorder", "Schizophrenia", 
  "Substance Dependence", "Generalized seizures", "Anemia",
  "Thrombophilia", "Hemophilia A")

use_diseases_sub <- diseases$diseaseId[diseases$diseaseName %in% 
        use_diseases_sub_name]
genes_diseases_sub <- disease2gene(disease = use_diseases_sub, database = "ALL", 
      verbose = TRUE, api_key = disgenet_api_key)
genes_diseases_sub <- split(genes_diseases_sub@qresult$gene_symbol, 
                        genes_diseases_sub@qresult$disease_name)

blue_colours <- c("#C6DBEF", "#6BAED6",  "#2171B5", "#08306B")

gene_trends$names <- paste0(gene_trends$major_trend, "..", 
                            gene_trends$major_clust) 
genes_interest <- split(gene_trends$gene_name, gene_trends$names)
genes_interest <- genes_interest[lengths(genes_interest)>=10]

genes_dev <- split(gene_trends$gene_name, gene_trends$major_clust)

res <- disease_enrich_celltype_trend(genes_diseases_sub, genes_interest, 
  genes_dev=genes_dev, sce=sce, input=input, exprs_thresh=0.05, 
  fdr_all=FALSE, back_cell=TRUE, same_range = FALSE)

res_plot <- res$res
rm_celltype_trend <- sapply(split(res_plot$padj, res_plot$celltype_trend), 
    function(x) all(is.na(x)))
res_plot <- res_plot[!res_plot$celltype_trend %in% 
    names(rm_celltype_trend)[rm_celltype_trend],]

res_plot$celltype <- recode(res_plot$celltype, "LAMP5_CA1"="LAMP5_NOS1",
"LAMP5_CCP2"="LAMP5_CCK",
"LAMP5_NMBR"="LAMP5_NDNF")

cell_types <- rev(c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4",
  "VIP", "ID2", "LAMP5_CA1", "SST", "PV", "PV_SCUBE3", "Astro",
  "OPC", "Oligo"))
general_trends <- c("up", "interup", "interdown", "down")
levels_trend <- paste0(rep(general_trends, each=length(cell_types)),
  "..", rep(cell_types, length(general_trends)))
res_plot$celltype_trend <- paste0(res_plot$trend, "..", res_plot$celltype)
res_plot$celltype_trend <- factor(res_plot$celltype_trend, levels=
    levels_trend)

order_xy <- find_hierarchical_order_1(res_plot, y=FALSE)
#order_becca <- c("Attention deficit hyperactivity disorder", "Amnesia", 
#  "Autism Spectrum Disorders", "Developmental regression",
#  "Generalized seizures", "Learning Disabilities", "Profound Mental Retardation",
#  "Alcoholic Intoxication, Chronic",
#   "Bipolar Disorder", "Major Depressive Disorder", "Manic", "Narcolepsy", 
#  "Nicotine Dependence", "Schizoaffective Disorder", "Schizophrenia", 
#  "Substance Dependence", "Alzheimer's Disease", "Dementia", "Huntington Disease", 
#  "Impaired cognition", "Mild cognitive disorder", "Anemia",
#  "Thrombophilia", "Hemophilia A")
res_plot$pathway <- factor(res_plot$pathway, levels=order_xy$x)


g1 <- ggplot(res_plot, aes(x=pathway, y=celltype_trend, fill=padj)) + 
  geom_tile() + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours = 
  rev(grDevices::colorRampPalette(colors=blue_colours)(20)), 
  na.value = "transparent") + 
  ggtitle("Disease enrichment") +
  background_grid()

ggsave(g1, file="main_figures/disease_enrichment_final.png",
       height=7, width=5.5)

for(i in 1:length(res$overlap)){
    
   ind <- res$res$celltype_trend==names(res$overlap)[i]
   ind1 <- !is.na(res$res$padj[ind])
   
   tmp <- res$overlap[[i]][ind1]
   
   if(!length(tmp)==0){
   
      for(ii in 1:length(tmp)){
        
        write.csv(tmp[[ii]],
          paste0("snRNAseq/processed_data/disease_enrichment_genes/", 
          names(res$overlap)[i], "_", names(tmp)[ii], ".csv"), 
          quote=F)
      }
   }
    
}


gene_trends$names <- paste0(gene_trends$gene_trend, "..", 
                            gene_trends$major_clust) 
genes_interest <- split(gene_trends$gene_name, gene_trends$names)
genes_interest <- genes_interest[lengths(genes_interest)>=10]

genes_dev <- split(gene_trends$gene_name, gene_trends$major_clust)

res <- disease_enrich_celltype_trend(genes_diseases_sub, genes_interest, 
  genes_dev=genes_dev, sce=sce, input=input, exprs_thresh=0.05, 
  fdr_all=FALSE, back_cell=TRUE, same_range = FALSE)

res_plot <- res$res
rm_celltype_trend <- sapply(split(res_plot$padj, res_plot$celltype_trend), 
    function(x) all(is.na(x)))
res_plot <- res_plot[!res_plot$celltype_trend %in% 
    names(rm_celltype_trend)[rm_celltype_trend],]


res_plot$celltype <- recode(res_plot$celltype, "LAMP5_CA1"="LAMP5_NOS1",
"LAMP5_CCP2"="LAMP5_CCK",
"LAMP5_NMBR"="LAMP5_NDNF")

cell_types <- rev(c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4",
  "VIP", "ID2", "LAMP5_NOS1", "SST", "PV", "PV_SCUBE3", "Astro",
  "OPC", "Oligo"))
general_trends <- paste0("G", 1:14)
levels_trend <- paste0(rep(general_trends, each=length(cell_types)),
  "..", rep(cell_types, length(general_trends)))
res_plot$celltype_trend <- paste0(res_plot$trend, "..", res_plot$celltype)
res_plot$celltype_trend <- factor(res_plot$celltype_trend, levels=
    levels_trend)

order_xy <- find_hierarchical_order_1(res_plot, y=FALSE)
#order_becca <- c("Attention deficit hyperactivity disorder", "Amnesia", 
#  "Autism Spectrum Disorders", "Developmental regression",
#  "Generalized seizures", "Learning Disabilities", "Profound Mental Retardation",
#  "Alcoholic Intoxication, Chronic",
#   "Bipolar Disorder", "Major Depressive Disorder", "Manic", "Narcolepsy", 
#  "Nicotine Dependence", "Schizoaffective Disorder", "Schizophrenia", 
#  "Substance Dependence", "Alzheimer's Disease", "Dementia", "Huntington Disease", 
#  "Impaired cognition", "Mild cognitive disorder", "Anemia",
#  "Thrombophilia", "Hemophilia A")
res_plot$pathway <- factor(res_plot$pathway, levels=order_xy$x)


g1 <- ggplot(res_plot, aes(x=pathway, y=celltype_trend, fill=padj)) + 
  geom_tile() + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8),
  axis.text.y = element_text(size=8)) +
  scale_fill_gradientn(colours = 
  rev(grDevices::colorRampPalette(colors=blue_colours)(20)), 
  na.value = "transparent") + 
  ggtitle("Disease enrichment") +
  background_grid()

ggsave(g1, file="main_figures/disease_enrichment_subtype_final.svg",
       height=7, width=5.5)

for(i in 1:length(res$overlap)){
    
   ind <- res$res$celltype_trend==names(res$overlap)[i]
   ind1 <- !is.na(res$res$padj[ind])
   
   tmp <- res$overlap[[i]][ind1]
   
   if(!length(tmp)==0){
   
      for(ii in 1:length(tmp)){
        
        write.csv(tmp[[ii]],
          paste0("snRNAseq/processed_data/disease_enrichment_genes_sub/", 
          names(res$overlap)[i], "_", names(tmp)[ii], ".csv"), 
          quote=F)
      }
   }
    
}

