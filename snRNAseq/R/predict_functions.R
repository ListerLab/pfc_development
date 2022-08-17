getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

predict_cells <- function(commap, test=commap$dataset=="query", 
                          train=commap$dataset=="ref", col, knn,
                          method="mode"){
    
    neighbours <- get.knnx(commap[train,1:2], commap[test,1:2], 
                           k=knn, algorithm="kd_tree")
    if(method=="mode"){
        predict <- apply(neighbours$nn.index, 1, function(x)
            getmode(commap[x, col]))
    } 
    if(method=="mean"){
        if(mode(commap[, col])!="numeric") stop("Not possible!")
        predict <- apply(neighbours$nn.index, 1, function(x)
            mean(commap[x, col]))
    }
    if(method=="median"){
        if(mode(commap[, col])!="numeric") stop("Not possible!")
        predict <- apply(neighbours$nn.index, 1, function(x)
            median(commap[x, col]))
    } 
    if(!method%in% c("median", "mode", "mean", "weighted_mode")){
        stop("Choose method: mode, mean, median, weighted_mode")
    }
    return(predict)
}


predict_all <- function(commap, test=commap$dataset=="query", 
         train=commap$dataset=="ref", knn, age_method="mode",
         method_all="mode") {
    
    age_predict <- as.numeric(predict_cells(commap, test=test, 
          train=train, "Age", knn=knn, method=age_method))
    if(method_all=="mode"){
        predict <- getmode(age_predict)
    } 
    if(method_all=="mean"){
        predict <- mean(age_predict)
    }
    if(method_all=="median"){
        predict <- median(age_predict)
    } 
    if(!method_all%in% c("median", "mode", "mean")){
        stop("Choose method: mode, mean, median")
    }
    
    return(predict)
    
}

predict_celltypes <- function(commap, test=commap$dataset=="query", 
                        train=commap$dataset=="ref", knn, age_method="mode",
                        method_celltypes="mode") {
    
    age_predict <- as.numeric(predict_cells(commap, test=test, 
         train=train, "Age", knn=knn, method=age_method))
    celltype_predict <- predict_cells(commap, test=test, 
        train=train, "Anno", knn=knn, method="mode")
    age_predict <- split(age_predict, celltype_predict)
    if(method_celltypes=="mode"){
        predict <- sapply(age_predict, function(x) getmode(x))
    } 
    if(method_celltypes=="mean"){
        predict <- sapply(age_predict, function(x) mean(x))
    }
    if(method_celltypes=="median"){
        predict <- sapply(age_predict, function(x) median(x))
    } 
    if(!method_celltypes%in% c("median", "mode", "mean")){
        stop("Choose method: mode, mean, median")
    }
    
    return(predict)
    
}

plot_histogram <- function(commap, test=commap$dataset=="query", 
                           train=commap$dataset=="ref", knn,
                           method="mode"){
    
    predict_age <-as.numeric(predict_cells(commap, test=test, 
          train=train, col="Age", knn=knn, method=method))
    predict_type <- predict_cells(commap, test=test, 
          train=train, col="Anno", knn=knn, method="mode")
    predict_stage <- predict_cells(commap, test=test, 
          train=train, "Stage", knn=knn, method="mode")
    
    predict_df <- data.frame(Age=predict_age, Type=predict_type, 
            Stage=predict_stage)
    
    ages_dat <- data.frame(c("ga22","ga38","60d","1yr","10yr","20yr"),
         c(-0.59206997, -0.06903606,  0.29393111,  1.35923668,  3.59432467,
         4.28690546))
    ages_dat <- ages_dat[order(as.numeric(ages_dat[,2])),]
    
    ggplot(predict_df, aes(x=Age, fill=Type)) + geom_histogram(col="black",
        bins=30) + 
        theme_classic() + scale_fill_manual(values=col) + scale_x_continuous(
            breaks=as.numeric(ages_dat[,2]), labels=ages_dat[,1])
}
