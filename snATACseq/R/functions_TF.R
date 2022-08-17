compute_enrich <- function(tmp_interest, bgdPeaks, matches){
    
    idx <- match(tmp_interest@elementMetadata$new_name, rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
        as.vector(bgdPeaks[idx,])))
    return(res)
    
}


compute_enrich_traj <- function(tmp_interest, bgdPeaks, matches, traj){
    
    names(tmp_interest) <- sapply(names(tmp_interest), function(x)
        strsplit(x, "id_", fixed=TRUE)[[1]][2])
    names(tmp_interest) <- paste0("id_", names(tmp_interest))
    idx <- match(names(tmp_interest), rownames(matches))
    
    res <- ArchR:::.computeEnrichment(matches, idx, c(idx, 
          as.vector(bgdPeaks[idx,])))
    return(res)
    
}

obtain_peak_ids <- function(peak_tmp, genes_cluster_trends, cluster_tmp, 
                            cell_types){
    
    cluster_tmp <- cell_types[cell_types[,2] %in% cluster_tmp, 1]
    ind <- sapply(cluster_tmp, function(x) 
        grepl(paste0(x, "."), names(genes_cluster_trends)))
    if(class(ind)[[1]]=="matrix"){
        ind <- apply(ind, 1, function(x) any(x))
    }
    tmp <- genes_cluster_trends[ind]
    ind <- lapply(tmp, function(x) unlist(sapply(x, function(y) 
        grep(y, peak_tmp$CRE))))
    ind <- lapply(ind, function(x) x[!is.na(x)])
    split_tmp <- lapply(ind, function(x) peak_tmp[x])
    split_tmp
}

obtain_peak_ids_dis <- function(disease_tmp, disease_name_tmp, peaks, 
                            cell_types){
    
    ind <- which(str_count(disease_name_tmp, cell_types[,1])==1)
    tmp <- peaks[[cell_types[ind[1],2]]]
    ind <- unlist(sapply(disease_tmp, function(y) 
        grep(y, tmp$CRE)))
    tmp[ind]
}

find_name <- function(i, tf_dat, index){
    
   tmp <- tf_dat[i,]
    tmp <- matrix(as.character(tmp), nrow=length(index[[i]]), 
                  ncol=dim(tf_dat)[2] ,byrow = TRUE)
    tmp <- as.data.frame(tmp)
    tmp$family_name <- motif_cluster$Name[match(
        motif_anno$Cluster_ID[index[[i]]], motif_cluster$Cluster_ID)]
    return(tmp)
}


motif_cluster <- read_xlsx("annotation/motif_annotations.xlsx", sheet=1)
motif_anno <- read_xlsx("annotation/motif_annotations.xlsx", sheet=2)
motif_anno <- motif_anno[motif_anno$Database=="Jaspar2018",]
motif_anno$Motif_short <- sapply(motif_anno$Motif, function(x)
    strsplit(x, "_")[[1]][1])
motif_anno$Motif_short <- toupper(motif_anno$Motif_short)

attach_family_name <- function(tf_dat, motif_cluster){

  index <- lapply(tf_dat$feature_short, function(x) 
      which(motif_anno$Motif_short %in% x))
  tf_dat <- tf_dat[!lengths(index)==0,]
  index <- index[!lengths(index)==0]
  
  tf_dat_new <- 
    lapply(1:length(index), function(x) find_name(x, tf_dat, index))
  tf_dat_new <- do.call(rbind, tf_dat_new)

  colnames(tf_dat_new)[1:(ncol(tf_dat_new)-1)] <- colnames(tf_dat)     
  
  return(tf_dat_new)
}

tf_not <- readRDS("snATACseq/processed_data/tf_not_rna.RDS")

remove_non_expressed_tf <- function(tf_dat, tf_not) {
  
  ind <- !(tf_dat$feature_short %in% tf_not)
  tf_dat <- tf_dat[ind,]
  
  return(tf_dat)
  
}

tf_not_cell_type <- readRDS("snATACseq/processed_data/tf_not_cell_type_rna.RDS")

cell_types <- matrix(c("Astro", "Astro",
                       "L2-3_CUX2", "L2_3",
                       "L5-6_THEMIS", "L5_6",
                       "L5-6_TLE4", "L5_6",
                       "L4_RORB", "L4",
                       "VIP", "CGE_der",
                       "LAMP5_CA1", "CGE_der",
                       "ID2", "CGE_der",
                       "SST", "MGE_der",
                       "PV", "MGE_der",
                       "PV_SCUBE3", "MGE_der",
                       "Oligo", "Oligo",
                       "OPC", "Oligo",
                       "Micro", "Micro"
), ncol=2, byrow=T)  

names_tf_not_cell_type_new <- cell_types[match(names(tf_not_cell_type), 
    cell_types[,1]),2]
index1 <- split(1:length(names_tf_not_cell_type_new), names_tf_not_cell_type_new)
tf_not_cell_type <- lapply(index1, function(x) Reduce(intersect, 
    tf_not_cell_type[x]))
  
remove_non_expressed_tf_cell_type <- function(tf_dat, tf_not_cell_type) {
  
  index <- split(1:nrow(tf_dat), tf_dat$trajectory)
  tmp <- lapply(index, function(x) tf_dat[x,])
  tf_not_cell_type <- tf_not_cell_type[names(tmp)]
  tf_dat_new <- mapply(function(X,Y) X[!X$feature_short %in% Y,], 
    X=tmp, Y=tf_not_cell_type, SIMPLIFY = FALSE)
  tf_dat_new <- do.call(rbind, tf_dat_new)
  
  return(tf_dat_new)
  
}

find_motif_order <- function(tmp) {

  a <- acast(tmp,  family_name ~ cell_type, mean, fill=1,
           drop=T, value.var="p.summary")
  hc <- hclust(dist(a))
  order_row <- order.optimal(dist(a), hc$merge)
  levels_motifs <- names(order_row$order)[order_row$order]
  
  return(levels_motifs)
  
}
