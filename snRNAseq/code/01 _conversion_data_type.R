###############################################################################
#                                                                             #
#                                                                             #
#  Make SingleCellExperiment from anndata                                     #
#                                                                             #
#                                                                             #    
###############################################################################

library(SingleCellExperiment)
library(rhdf5)
library(Matrix)

link <- "paths_to_h5ad_folder/"
file_name <- "2020-12-18_whole-tissue_post-restaged-GABA-clustering.h5ad"
h5ls(paste0(link, file_name))

xx <- h5read(paste0(link, file_name), "X/data")
xx <- c(xx)
indices <- h5read(paste0(link, file_name), "X/indices")
indices <- c(indices)
indptr <- h5read(paste0(link, file_name), "X/indptr")
indptr <- c(indptr)
obs__categories <-  h5read(paste0(link, file_name), "obs/__categories")
obs <-  h5read(paste0(link, file_name), "obs")
var <- h5read(paste0(link, file_name), "var")
var__categories <- h5read(paste0(link, file_name), "var/__categories")

sce <- sparseMatrix(i = indices+1 , p=indptr, x=xx)
rownames(sce) <- c(var$`_index`)
colnames(sce) <- c(obs$`_index`)

rm("xx", "indices", "indptr")
gc()

var <- var[-which(names(var)=="__categories")]
names_var <- names(var)
names_var[which(names_var=="_index")] <- "index"

func <- paste0("rowData <- DataFrame(",
        paste0(paste0(names_var, "= c(var[[", 1:length(names_var),"]])",
               collapse=","), ")"))

eval(parse(text=func))

obs <- obs[-which(names(obs)=="__categories")]
names_obs <- names(obs)
names_obs[which(names_obs=="_index")] <- "index"
names_obs <- gsub(" ", "_", names_obs)
names_obs <- gsub("*", "", names_obs, fixed=TRUE)
names_obs <- gsub("-", "_", names_obs, fixed=TRUE)


func <- paste0("colData <- DataFrame(",
        paste0(paste0(names_obs, "= c(obs[[", 1:length(names_obs),"]])",
                 collapse=","), ", check.names=F)"))

eval(parse(text=func))
colnames(colData) <- sapply(colnames(colData), function(x)
        strsplit(x, ".", fixed=T)[[1]][1])


names_obs <- names(obs__categories)
names_obs <- gsub(" ", "_", names_obs)
names_obs <- gsub("*", "", names_obs, fixed=TRUE)
names_obs <- gsub("-", "_", names_obs, fixed=TRUE)

change <- which(sapply(obs__categories, class)!="array")

if(length(change)>0){
    
    for(i in change){
        
        obs__categories[[i]] <- unlist(obs__categories[[i]])
    }
    
}


#convert the categorical variables to factors


funcs <- lapply(1:length(names_obs), function(xx)
    paste0("colData$", names_obs[xx], "<- as.factor(c(obs__categories[[",
           xx, "]][as.numeric(as.character(colData$", names_obs[xx], "))+1]))"))

for(i in funcs){ 
    eval(parse(text=i))
}

sc <- SingleCellExperiment(assays=SimpleList(counts=sce), 
   rowData=rowData, colData=colData)

sc <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

red_UMAP <-  h5read(paste0(link, file_name), "obsm")

reducedDim(sc, "PCA") <- t(as.matrix(red_UMAP$X_pca))
reducedDim(sc, "UMAP") <- t(as.matrix(red_UMAP$X_umap))
reducedDim(sc, "UMAT") <- t(as.matrix(red_UMAP$X_umat))

rowData(sc)$highly_variable<-rowData(sc)$highly_variable>1

outlink <- 'snRNAseq/processed_data'

saveRDS(sc, file=paste0(outlink, gsub(".h5ad",".RDS", file_name)))

