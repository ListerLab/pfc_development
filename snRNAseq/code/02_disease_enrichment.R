###############################################################################
#                                                                             #
#                                                                             #
# Disease enrichment analysis for different stages                            #
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

blue_colours <- c("#C6DBEF", "#6BAED6",  "#2171B5", "#08306B")

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
  "Developmental regression", "Profound Mental Retardation", 
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

## 3. cell type and major trend - background cell type genes

res <- disease_enrich_celltype_trend(genes_diseases_sub, genes_interest, 
  genes_dev=genes_dev, sce=sce, input=input, exprs_thresh=0.05, 
  fdr_all=FALSE, back_cell=TRUE, same_range = TRUE)

res_plot <- res$res
rm_celltype_trend <- sapply(split(res_plot$padj, res_plot$celltype_trend), 
    function(x) all(is.na(x)))
res_plot <- res_plot[!res_plot$celltype_trend %in% 
    names(rm_celltype_trend)[rm_celltype_trend],]

cell_types <- c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4",
  "VIP", "ID2", "LAMP5_CA1", "SST", "PV", "PV_SCUBE3", "Astro",
  "OPC", "Oligo")
general_trends <- c("up", "interup", "interdown", "down")
levels_trend <- paste0(rep(cell_types, each=length(general_trends)),
  "..", rep(general_trends, length(cell_types)))
res_plot$celltype_trend <- factor(res_plot$celltype_trend, levels=
    levels_trend)

order_xy <- find_hierarchical_order_1(res_plot, y=FALSE)
res_plot$pathway <- factor(res_plot$pathway, levels=order_xy$x)


g10 <- ggplot(res_plot, aes(x=pathway, y=celltype_trend, fill=padj)) + 
  geom_tile() + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours = 
  rev(grDevices::colorRampPalette(colors=blue_colours)(20)), 
  na.value = "transparent") + 
  ggtitle("Disease enrichment") +
  background_grid()

ggsave(g10, file="main_figures/disease_enrichment_final.svg",
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

# 4. cell-specific DEGs adult stage

deg <- read.table("snRNAseq/processed_data/degs_major_traj/results_file.txt",
  sep="\t")
colnames(deg) <- gsub("L5.6", "L5_6", colnames(deg) )
colnames(deg)  <- gsub("L2.3", "L2_3", colnames(deg) )

all_combos <- colnames(deg)[grepl("Res.", colnames(deg))]
all_combos <- gsub("Res.", "", all_combos)
all_combos <- sapply(all_combos, function(x) strsplit(x, ".", fixed=TRUE)[[1]])
all_combos <- t(all_combos)

de_list <- lapply(rownames(all_combos), function(x) 
  deg[,grepl(paste0(x, "$"), colnames(deg))])
for(i in 1:length(de_list)) {
  
  colnames(de_list[[i]]) <- c("Coef", "t", "p.value", "Res")
  
}

celltype_markers <- combineMarkers(de_list, all_combos, pval.field = "p.value",
  effect.field = "Coef", pval.type = "some")
genes_interest <- lapply(celltype_markers, function(x) 
  rownames(x)[x$FDR<0.05 & x$summary.Coef>5])
names(genes_interest)[c(5,6,8)] <- c("L2-3_CUX2", "L5-6_TLE4", "L5-6_THEMIS")
genes_back <- rownames(deg)

res <- disease_enrich_celltype_specific(genes_diseases_sub, genes_interest, 
    genes_back, fdr_all=FALSE) 

res$res$celltype <- recode(res$res$celltype, "LAMP5_CA1"="LAMP5_NOS1",
"LAMP5_CCP2"="LAMP5_CCK",
"LAMP5_NMBR"="LAMP5_NDNF")

order_xy <- find_hierarchical_order(res$res)
res$res$pathway <- factor(res$res$pathway, levels=order_xy$x)
res$res$celltype <- factor(res$res$celltype, levels=order_xy$y)

g1 <- ggplot(res$res, aes(x=pathway, y=celltype, fill=padj)) + 
  geom_tile() + theme_classic() +
  theme(axis.title = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_gradientn(colours = 
  rev(grDevices::colorRampPalette(colors=blue_colours)(20)), 
  na.value = "transparent") + 
  ggtitle("Cell type specific expression") +
  background_grid()


g1 

ggsave(g1, file="supp_figures/disease_enrichment_celltype_specific_adult.png")



