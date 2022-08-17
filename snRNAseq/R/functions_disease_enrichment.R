# Functions disease enrichment

input <- list(Astro=c(which(sce$`Astro_GFAP_dev-traj`=="dev-traj"), 
    which(sce$`Astro_SLC1A2_dev-traj`=="dev-traj")), 
    Oligo=c(which(sce$`Oligo_dev-traj`=="dev-traj"), which(
    sce$`OPC_MBP_dev-traj`=="dev-traj")),
    OPC=c(which(sce$`OPC_MBP_dev-traj`=="dev-traj"),
    which(sce$`OPC_dev-traj`=="dev-traj")),
    Micro=c(which(sce$`Micro_dev-traj`=="dev-traj")),
    L4_RORB=c(which(sce$`L4_RORB_LRRK1_dev-traj`=="dev-traj"),
    which(sce$`L4_RORB_MET_dev-traj`=="dev-traj"),
    which(sce$`L4_RORB_MME_dev-traj`=="dev-traj")),
    `L5-6_TLE4`=c(which(sce$`L5/6_TLE4_HTR2C_dev-traj`=="dev-traj"),
     which(sce$`L5/6_TLE4_SCUBE1_dev-traj`=="dev-traj"),
     which(sce$`L5/6_TLE4_SORCS1_dev-traj`=="dev-traj")),
    `L5-6_THEMIS`=c(which(sce$`L5/6_THEMIS_CNR1_dev-traj`=="dev-traj"), 
      which(sce$`L5/6_THEMIS_NTNG2_dev-traj`=="dev-traj")),
    `L2-3_CUX2`=c(which(sce$`L2_CUX2_LAMP5_dev-traj`=="dev-traj"),
    which(sce$`L3_CUX2_PRSS12_dev-traj`=="dev-traj")),
    `LAMP5_CA1`=c(which(sce$`LAMP5_LCP2_dev-traj`=="dev-traj"),
    which(sce$`LAMP5_NMBR_dev-traj`=="dev-traj"),
    which(sce$`LAMP5_CA1_dev-traj`=="dev-traj")),
    ID2= c(which(sce$`ID2_CSMD1_dev-traj`=="dev-traj"),
    which(sce$`CCK_MC4R_dev-traj`=="dev-traj"),
    which(sce$`CCK_RELN_dev-traj`=="dev-traj"),
    which(sce$`CCK_SYT6_dev-traj`=="dev-traj")),
    VIP= c(which(sce$`VIP_ABI3BP_dev-traj`=="dev-traj"),
    which(sce$`VIP_ADAMTSL1_dev-traj`=="dev-traj"),
    which(sce$`VIP_CHRM2_dev-traj`=="dev-traj"),
    which(sce$`VIP_CRH_dev-traj`=="dev-traj"),
    which(sce$`VIP_DPP6_dev-traj`=='dev-traj'),
    which(sce$`VIP_HS3ST3A1_dev-traj`=='dev-traj'),
    which(sce$`VIP_KIRREL3_dev-traj`=="dev-traj")),
    `PV_SCUBE3`=c(which(sce$`PV_SCUBE3_dev-traj`=="dev-traj")), 
    `PV`=c(which(sce$`PV_SST_dev-traj`=="dev-traj"),
     which(sce$`PV_SULF1_dev-traj`=="dev-traj"),
     which(sce$`PV_WFDC2_dev-traj`=="dev-traj")),
    `SST`=c(which(sce$`SST_ADGRG6_dev-traj`=="dev-traj"),
     which(sce$`SST_B3GAT2_dev-traj`=="dev-traj"),
     which(sce$`SST_BRINP3_dev-traj`=="dev-traj"),
     which(sce$`SST_CALB1_dev-traj`=="dev-traj"),
     which(sce$`SST_NPY_dev-traj`=="dev-traj"),
     which(sce$`SST_STK32A_dev-traj`=="dev-traj"),
     which(sce$`SST_TH_dev-traj`=="dev-traj")),
    Vas=c(which(sce$`Vas_CLDN5_dev-traj`=="dev-traj"),
        which(sce$`Vas_PDGFRB_dev-traj`=="dev-traj"),
        which(sce$`Vas_TBX18_dev-traj`=="dev-traj")))

# find background genes

find_back_genes <- function(sce, exprs_thresh, input, cell_types){
  
  input1 <- input[cell_types]
  no_exprs <- lapply(input1, function(x) scater::nexprs(sce[,x], byrow=TRUE))
  tmp_name <- mapply(x=no_exprs, y=lengths(input1), function(x,y) 
    names(x)[x>=y*exprs_thresh], SIMPLIFY = FALSE)
  names(tmp_name) <- cell_types
  
  return(tmp_name)
  
}

same_range <- function(genes_back, genes_interest, sce, input, cell_types){
  
  input <- input[cell_types]
  all_ranges <- mapply(X=input, Y=genes_interest, function(X,Y) {
    tmp <- counts(sce[rownames(sce) %in% Y, X])
    tmp <- rowMeans(tmp)
    c(min(tmp), max(tmp))
  }, SIMPLIFY = FALSE)
  means_back <- mapply(X=input, Y=genes_back, function(X,Y){
    tmp <- counts(sce[rownames(sce) %in% Y, X])
    tmp <- rowMeans(tmp)
    tmp
  })
  gene_back <- mapply(X=means_back, Y=all_ranges, function(X,Y){
    names(X)[X>=Y[1] & X<=Y[2]]
  })
  return(gene_back)
  
}

# disease enrichment in cell-type versus expressed in all

disease_enrich_all <- function(genes_diseases_sub, genes_interest, 
    sce, input=input, exprs_thresh, fdr_all=TRUE) {
  
  # prepare background
  
  genes_back <- nexprs(sce, byrow=TRUE)
  genes_back <- names(genes_back)[genes_back>dim(sce)[2]*exprs_thresh]
  genes_back <- rep(list(genes_back), length(genes_interest))
  
  # do enrichment
  
  res <- mapply(int_genes=genes_interest, back_genes=genes_back, 
      function(int_genes, back_genes) fora(genes_diseases_sub, int_genes, 
      back_genes, minSize = 1, maxSize = Inf),  SIMPLIFY=FALSE)
  names(res) <- names(genes_interest)
  
  # make matrix and split out overlap genes
  
  overlap_genes <-lapply(res, function(x) x$overlapGenes)
  
  for(i in 1:length(res)[[1]]){
    
    res[[i]] <- res[[i]][, -"overlapGenes"]
    res[[i]] <- as.data.frame(res[[i]]) 
    res[[i]]$celltype <- names(res)[i]
    
  } 
  
  mat_res <- do.call(rbind, res)
  
  if(fdr_all){
  
    mat_res$padj <- p.adjust(mat_res$pval, method="fdr")
  
  }
  mat_res$padj[mat_res$padj>0.05] <- NA 
  
  return(list(res=mat_res, overlap=overlap_genes))
  
}
  

# disease enrichment for devDEGs versus expressed in cell type

disease_enrich_celltype <- function(genes_diseases_sub, genes_interest, 
    sce, input=input, exprs_thresh, fdr_all=TRUE, same_range=FALSE) {
  
  # prepare background
  
  cell_types <- names(genes_interest)
  genes_back <- find_back_genes(sce=sce, exprs_thresh=exprs_thresh, 
      input=input, cell_types=cell_types)
  
  if(same_range){
    
    genes_back <- same_range(genes_back, genes_interest, sce, input, cell_types)
  }
  
  # do enrichment
  
  res <- mapply(int_genes=genes_interest, back_genes=genes_back, 
      function(int_genes, back_genes) fora(genes_diseases_sub, int_genes, 
      back_genes, minSize = 1, maxSize = Inf),  SIMPLIFY=FALSE)
  names(res) <- cell_types
  
  # make matrix and split out overlap genes
  
  overlap_genes <- lapply(res, function(x) x$overlapGenes)
  for(i in 1:length(overlap_genes)){
    
    names(overlap_genes[[i]]) <- 
      res[[i]]$pathway
    
  }
  
  for(i in 1:length(res)[[1]]){
    
    res[[i]] <- res[[i]][, -"overlapGenes"]
    res[[i]] <- as.data.frame(res[[i]]) 
    res[[i]]$celltype <- names(res)[i]
    
  } 
  
  mat_res <- do.call(rbind, res)
  
  if(fdr_all){
  
    mat_res$padj <- p.adjust(mat_res$pval, method="fdr")
  
  }
  mat_res$padj[mat_res$padj>0.05] <- NA 
  
  return(list(res=mat_res, overlap=overlap_genes))
  
}

# disease enrichment of devDEGs in trend versus background

disease_enrich_celltype_trend <- function(genes_diseases_sub, genes_interest, 
  genes_dev, sce, input=input, exprs_thresh, fdr_all=TRUE, back_cell=FALSE,
  same_range=FALSE) {
  
  # prepare background
  
  cell_types <- sapply(names(genes_interest), function(x) 
    strsplit(x, "..", fixed=TRUE)[[1]][2])
  
  if(back_cell) {
  
    genes_back <- find_back_genes(sce=sce, exprs_thresh=exprs_thresh, 
      input=input, cell_types=names(genes_dev))
  
  } else {
    
    genes_back <- genes_dev[cell_types]
    
  }
  
  if(same_range){
    
    genes_back <- same_range(genes_back, genes_interest, sce, input, cell_types)
  }
  
  
  # do enrichment
  
  res <- mapply(int_genes=genes_interest, back_genes=genes_back, 
   function(int_genes, back_genes) fora(genes_diseases_sub, int_genes, 
   back_genes, minSize = 1, maxSize = Inf),  SIMPLIFY=FALSE)
  names(res) <- names(genes_interest)
  
  # make matrix and split out overlap genes
  
  overlap_genes <- lapply(res, function(x) x$overlapGenes)
  
  for(i in 1:length(overlap_genes)){
    
    names(overlap_genes[[i]]) <- 
      res[[i]]$pathway
    
  }
  
  
  for(i in 1:length(res)[[1]]){
    
    res[[i]] <- res[[i]][, -"overlapGenes"]
    res[[i]] <- as.data.frame(res[[i]]) 
    res[[i]]$celltype <- strsplit(names(res)[i], "..", fixed=TRUE)[[1]][2]
    res[[i]]$trend <- strsplit(names(res)[i], "..", fixed=TRUE)[[1]][1]
    res[[i]]$celltype_trend <- names(res)[i]
    
  } 
  
  mat_res <- do.call(rbind, res)
  
  if(fdr_all){
    
    mat_res$padj <- p.adjust(mat_res$pval, method="fdr")
    
  } else {
    
    index_split <- split(1:nrow(mat_res), mat_res$celltype)
    for(i in index_split){
      mat_res$padj[i] <- p.adjust(mat_res$pval[i], method="fdr")
    }
    
  }
  
  
  mat_res$padj[mat_res$padj>0.05] <- NA 
  
  return(list(res=mat_res, overlap=overlap_genes))
  
}
  
disease_enrich_celltype_specific <- function(genes_diseases_sub, genes_interest, 
    genes_back, fdr_all=TRUE) {
  
  # prepare background
  
   if(length(genes_back) != length(genes_interest)){
  
    genes_back <- rep(list(genes_back), length(genes_interest))
  
  }
  
  # do enrichment
  
  res <- mapply(int_genes=genes_interest, back_genes=genes_back, 
      function(int_genes, back_genes) fora(genes_diseases_sub, int_genes, 
      back_genes, minSize = 1, maxSize = Inf),  SIMPLIFY=FALSE)
  names(res) <- names(genes_interest)
  
  # make matrix and split out overlap genes
  
  overlap_genes <-lapply(res, function(x) x$overlapGenes)
  
  for(i in 1:length(res)[[1]]){
    
    res[[i]] <- res[[i]][, -"overlapGenes"]
    res[[i]] <- as.data.frame(res[[i]]) 
    res[[i]]$celltype <- names(res)[i]
    
  } 
  
  mat_res <- do.call(rbind, res)
  
  if(fdr_all){
  
    mat_res$padj <- p.adjust(mat_res$pval, method="fdr")
  
  }
  mat_res$padj[mat_res$padj>0.05] <- NA 
  
  return(list(res=mat_res, overlap=overlap_genes))
  
}
   
# figure out hierarchical order

find_hierarchical_order <- function(tmp, x=TRUE, y=TRUE){
  
  tmp <- acast(tmp, celltype~pathway, value.var="padj")
  tmp[is.na(tmp)] <- 1
  
  if(!x){
    
    tmp <- tmp[x,]
    
  } 
  
  if(!y) {
    
    tmp <- tmp[,y]
    
  }
  
  if(y){
    
    hc <- hclust(dist(tmp))
    order_row <- order.optimal(dist(tmp), hc$merge)
    y <- names(order_row$order)[order_row$order]
    tmp <- tmp[y,]
    
  }
  
  if(x){
    
    hc <- hclust(dist(t(tmp)))
    order_row <- order.optimal(dist(t(tmp)), hc$merge)
    x <- names(order_row$order)[order_row$order]
    tmp <- tmp[,x]
    
  }
  
  return(list(x=x, y=y))
  
  
}


find_hierarchical_order_1 <- function(tmp, x=TRUE, y=TRUE){
  
  tmp <- acast(tmp, celltype_trend~pathway, value.var="padj")
  tmp[is.na(tmp)] <- 1

  if(y){
    
    hc <- hclust(dist(tmp))
    order_row <- order.optimal(dist(tmp), hc$merge)
    y <- names(order_row$order)[order_row$order]
    tmp <- tmp[y,]
    
  }
  
  if(x){
    
    hc <- hclust(dist(t(tmp)))
    order_row <- order.optimal(dist(t(tmp)), hc$merge)
    x <- names(order_row$order)[order_row$order]
    tmp <- tmp[,x]
    
  }
  
  return(list(x=x, y=y))
  
  
}
