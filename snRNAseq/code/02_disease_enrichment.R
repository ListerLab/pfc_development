###############################################################################
#                                                                             #
#                                                                             #
# Disease enrichment analysis for different stages                            #
#                                                                             #
#                                                                             #    
###############################################################################

library(data.table)
library(disgenet2r)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(viridis)
library(scater)
library(cba)
library(dplyr)
disgenet_api_key <- ""

sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

diseases <- fread("annotation/all_diseases.tsv", header=T, data.table=F)

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")
cell_types <- unique(sce$major_clust)

make_disease_enrich_mat <- function(genes_diseases_sub, genes_stage_clust,
        stages, cell_types) {

    all_res <- list()

    for(j in 1:length(cell_types)){
    
        print(cell_types[j])
    
        all_res[[cell_types[j]]] <- list()
        all_genes_stages <- genes_stage_clust[grepl(cell_types[j], 
                names(genes_stage_clust))]
        all_genes <- unique(unlist(all_genes_stages))
      
    
        for(i in 1:length(stages)){
        
            print(stages[i])
            tmp_genes <- all_genes_stages[grepl(paste0("_", stages[i], "$"), 
                names(all_genes_stages))] 
            tmp_genes <- unlist(tmp_genes)
            if(length(tmp_genes)<10) {
                
                all_res[[cell_types[j]]][[paste0("Cluster_", stages[i])]] <- NULL
                
            } else {
        
                over <- fora(genes_diseases_sub, tmp_genes, all_genes, 
                     minSize = 1, maxSize = Inf)
        
                all_res[[cell_types[j]]][[paste0("Cluster_", stages[i])]] <- over
            }
        }
        
    }

    all_res2 <- all_res

    for(j in 1:length(cell_types)){
    
        for(i in 1:length(stages)){
        
            if(is.null(all_res2[[cell_types[j]]][[paste0("Cluster_", stages[i])]])) next
        
            all_res2[[cell_types[j]]][[paste0("Cluster_", stages[i])]]$Cell_type <- cell_types[j]
            all_res2[[cell_types[j]]][[paste0("Cluster_", stages[i])]]$Stage <- stages[i]
        
        }
    
        all_res2[[cell_types[j]]] <- all_res2[[cell_types[j]]][!
            sapply(all_res2[[cell_types[j]]], function(x) is.null(x))]
    }

    all_res2 <- lapply(all_res2, function(x) Reduce(rbind, x))
    all_res2 <- Reduce(rbind, all_res2)
    all_res2$Cell_type_stages <- paste0(all_res2$Cell_type, "_", all_res2$Stage)
    sample_order <- paste0(rep(cell_types, each=length(stages)), "_", 
        rep(rev(stages), length(cell_types)))
    all_res2$Cell_type_stages <- factor(all_res2$Cell_type_stages, 
                                    levels=sample_order)
    miss <- setdiff(sample_order, unique(all_res2$Cell_type_stages))
    all_diseases <- unique(all_res2$pathway)
    dat <- as.data.frame(matrix(NA, nrow=length(miss)*length(all_diseases), 
                            ncol=ncol(all_res2)))
    colnames(dat) <- colnames(all_res2)
    dat$Cell_type_stages <- rep(miss, each=length(all_diseases))
    dat$pathway <- rep(all_diseases, length(miss))

    all_res2 <- rbind(all_res2, dat)
    
    return(all_res2)
    
}

# Gene-disease associations
use_diseases_sub_name <- diseases$diseaseName[grepl("F03", diseases$diseaseClass) & 
                                           diseases$diseaseType=="disease"]
use_diseases_sub_name <- c("Alcoholic Intoxication, Chronic", 
    "Alzheimer's Disease", "Amnesia", "Attention deficit hyperactivity disorder", 
    "Autism Spectrum Disorders", "Bipolar Disorder", "Dementia", 
    "Developmental regression", 
    "Huntington Disease", "Impaired cognition", "Major Depressive Disorder",
    "Manic", "Mild cognitive disorder", "Narcolepsy", "Nicotine Dependence",
    "Schizoaffective Disorder", "Schizophrenia", "Substance Dependence")

use_diseases_sub <- diseases$diseaseId[diseases$diseaseName %in% 
        use_diseases_sub_name & diseases$diseaseType=="disease"]
genes_diseases_sub <- disease2gene(disease = use_diseases_sub, database = "ALL", 
      verbose = TRUE, api_key = disgenet_api_key)
genes_diseases_sub <- split(genes_diseases_sub@qresult$gene_symbol, 
                        genes_diseases_sub@qresult$disease_name)

# Genes expressed in gene trends

gene_trends <- read.csv("snRNAseq/processed_data/gene_cluster_ids.csv")
gene_trends$names <- paste0(gene_trends$major_clust, "_", 
        gene_trends$gene_trend) 
gene_stage_clust <- split(gene_trends$gene_name, gene_trends$names)

gene_stage_clust <- lapply(gene_stage_clust, function(x) 
    rowData(sce)$index[match(x, rownames(sce))])

cell_types <- unique(gene_trends$major_clust)
stages <- unique(gene_trends$gene_trend)
all_res2 <- make_disease_enrich_mat(genes_diseases_sub, gene_stage_clust,
                stages, cell_types)
all_res2$padj <- p.adjust(all_res2$pval, method="fdr") 
all_res3 <- apply(all_res2, 2, function(x) as.character(x))
write.table(all_res3[,-6], file="snRNAseq/processed_data/Disease_enrichment.tsv", sep="\t", 
    quote=F)
all_res2$padj[all_res2$padj>0.05] <- NA
all_res2$padj[all_res2$overlap<5] <- NA

system("mkdir snRNAseq/processed_data/disease_enrichment_genes")

split_row_order <- as.character(split_row_order)
for(i in split_row_order){
    
    tmp <- all_res2[all_res2$Cell_type_stages==i,]
    tmp <- tmp[!is.na(tmp$padj),]
    
    for(ii in 1:nrow(tmp)){
        
        write.csv(tmp$overlapGenes[ii],
            paste0("snRNAseq/processed_data/disease_enrichment_genes/", i, "_", 
                   tmp$pathway[ii], ".csv"), quote=F)
    }
    
}

