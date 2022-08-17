library(scater)
library(scran)
library(SingleCellExperiment)
library(batchelor)
library(FNN)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)

path1 <- "snRNAseq/processed_data/orgpred/"
all_samples <- c("Luc9228.RDS", "org2290.RDS", "org2432.RDS")

sce <- readRDS("snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")
umap <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")
reducedDim(sce, "UMAP") <- reducedDim(umap, "UMAP")

sce <- logNormCounts(sce)

rm(umap)
gc()

source("scatac/R/functions_nmf.R")

walk_along <- function (neigh){
  
  final_vec <- vector(length = 0)
  for(i in 1:nrow(neigh)){
      
      if(i==1){
        
        final_vec <- neigh[i,1]
        
      } else {
        
        tmp <- neigh[i,which(!neigh[i,] %in% final_vec)[1]]
        final_vec <- c(final_vec, tmp)
        
      }
  }
  
  return(final_vec)
}

# Compare mature neuronal organoid cells to similar maturity brain calls
## Comparison in individual datasets

all_markers_PN <- list()
sce_organoids <- list()
all_neighbors <- list()

for(i in all_samples){
    
    sce_organoids[[i]] <- readRDS(paste0(path1, i))
    sce_organoids[[i]] <- scater::logNormCounts(sce_organoids[[i]])
    umap <- readRDS(paste0("snRNAseq/processed_data/orgpred/outumaps/",
        gsub(".RDS", "", i), "_umap.RDS"))
    rownames(umap) <- umap$barcode
    reducedDim(sce_organoids[[i]], "UMAP") <- umap[,-3]

    ind_PN <- sce_organoids[[i]]$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
        "L5/6_THEMIS_TLE4", "PN_dev")
    ind_mature <- sce_organoids[[i]]$predstage != "Fetal"
    
    neighbors <- get.knnx(reducedDim(sce[, sce$cell_type=="PN"], "UMAP"),
      reducedDim(sce_organoids[[i]], "UMAP")[ind_PN&ind_mature,], k=50, 
      algorithm="kd_tree")
    neighbors <- walk_along(neighbors$nn.index)
    all_neighbors[[i]] <- colnames(sce[, sce$cell_type=="PN"])[neighbors] 
    rownames(sce_organoids[[i]]) <- rowData(sce_organoids[[i]])$Symbol
    
    all_genes <- intersect(rownames(sce), rownames(sce_organoids[[i]]))
    tmp_org <- sce_organoids[[i]][all_genes,ind_PN&ind_mature]
    reducedDim(tmp_org) <- NULL
    colData(tmp_org) <- DataFrame(barcode=tmp_org$Barcode,
      Condition="Organoid")
    tmp_sce <- sce[all_genes, neighbors]
    reducedDim(tmp_sce, "UMAP") <- NULL
    reducedDim(tmp_sce, "PCA") <- NULL
    reducedDim(tmp_sce, "UMAT") <- NULL
    colData(tmp_sce) <- DataFrame(barcode=colnames(tmp_sce),
      Condition="Brain")
    
    new_sce <- cbind(tmp_org, tmp_sce)
    
    markers <- findMarkers(new_sce, groups=new_sce$Condition)
    all_markers_PN[[i]] <- markers$Organoid
    
}

all_markers_PN_up <- lapply(all_markers_PN, function(x)
    rownames(x)[x$FDR<0.05 & x$summary.logFC>0])
table_up_PN <- sort(table(unlist(all_markers_PN_up)))
head(table_up_PN )
all_markers_PN_down <- lapply(all_markers_PN, function(x)
    rownames(x)[x$FDR<0.05 & x$summary.logFC<0])
table_down_PN <- sort(table(unlist(all_markers_PN_down)))
head(table_down_PN )

## Comparison in joined datasets

all_sce <- lapply(all_samples, function(x) readRDS(paste0(path1, x)))
all_sce <- multiBatchNorm(all_sce)

all_sce <- do.call(cbind, all_sce)
rownames(all_sce) <- rowData(all_sce)$Symbol

ind_PN <- all_sce$predcelltype %in% c("L2/3_CUX2", "L4_RORB",
                                  "L5/6_THEMIS_TLE4", "PN_dev")
ind_mature <- all_sce$predstage != "Fetal"


all_genes <- intersect(rownames(sce), rownames(all_sce))
tmp_org <- all_sce[all_genes,ind_PN&ind_mature]
reducedDim(tmp_org) <- NULL
colData(tmp_org) <- DataFrame(barcode=tmp_org$Barcode,
      Condition="Organoid")
tmp_sce <- sce[, sce$cell_type=="PN"]
tmp_sce <- tmp_sce[all_genes,  unique(unlist(all_neighbors))]
reducedDim(tmp_sce, "UMAP") <- NULL
reducedDim(tmp_sce, "PCA") <- NULL
reducedDim(tmp_sce, "UMAT") <- NULL
colData(tmp_sce) <- DataFrame(barcode=colnames(tmp_sce),
      Condition="Brain")

new_sce <- cbind(tmp_org, tmp_sce)
    
markers <- findMarkers(new_sce, groups=new_sce$Condition)

write.table(markers$Organoid, 
      file="snRNAseq/processed_data/Organoids_Brain_markers_PN.csv",
            quote=F, sep=",")

markers_up <- rownames(markers$Organoid)[markers$Organoid$FDR<0.05 & 
        markers$Organoid$summary.logFC>0]
markers_down <- rownames(markers$Organoid)[markers$Organoid$FDR<0.05 & 
        markers$Organoid$summary.logFC<0]

go_term_analysis <- function(genes, genes_background){
    
    
    genes <-  AnnotationDbi::select(org.Hs.eg.db, keys=genes, columns="ENTREZID", 
                                    keytype="SYMBOL")[,2]
    genes <- genes[!is.na(genes)]
    genes_background <- AnnotationDbi::select(org.Hs.eg.db, keys=genes_background, 
      columns="ENTREZID", keytype="SYMBOL")[,2]
    genes_background <- genes_background[!is.na(genes_background)]
    
    go_bp <- enrichGO(gene= genes,
                      universe = genes_background,
                      OrgDb = org.Hs.eg.db,
                      ont  = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff  = 0.05,
                      readable = TRUE)
    return(go_bp)
    
}


go_pn_up <- go_term_analysis(markers_up, rownames(markers$Organoid))
go_pn_down <- go_term_analysis(markers_down, rownames(markers$Organoid))

go_pn_up_tab <- go_pn_up@result
go_pn_down_tab <- go_pn_down@result

go_pn_up_tab$geneID <- unlist(go_pn_up_tab$geneID)
go_pn_down_tab$geneID <- unlist(go_pn_down_tab$geneID)

write.table(go_pn_up_tab, file="snRNAseq/processed_data/GO_brain_organoid_up.tsv",
    sep="\t", quote=FALSE, col.names = TRUE)
write.table(go_pn_down_tab, file="snRNAseq/processed_data/GO_brain_organoid_down.tsv",
    sep="\t", quote=FALSE, col.names = TRUE)

ind <- match(unlist(all_neighbors), colnames(sce))
umap_org <- reducedDim(sce, "UMAP")
umap_organoid <- do.call(rbind, lapply(sce_organoids, function(x) 
  reducedDim(x, "UMAP")))
umap_organoid <- umap_organoid[ind_mature&ind_PN,]
umap_org <- as.data.frame(umap_org)  
colnames(umap_org) <- c("UMAP1", "UMAP2")
umap_organoid <- as.data.frame(umap_organoid)

g1 <- ggplot(umap_org[,1:2], aes(x=UMAP1, y=UMAP2)) +
    geom_point(col="black", size=1.5) +
    geom_point(col="white", size=0.8) +
    theme_classic() +
    geom_point(data=umap_org[ind,1:2], 
        alpha=0.5, size=1, col="cornflowerblue") +
    geom_point(data=umap_organoid[,1:2], alpha=0.5, size=0.5, col="#990099")


ggsave(g1, file="main_figures/ogranoids_brain_umap_PN.png", width=7, height=7)

saveRDS(umap_org, file="snRNAseq/processed_data/umap_orginial.RDS")
saveRDS(umap_organoid, file="snRNAseq/processed_data/umap_organoid.RDS")
