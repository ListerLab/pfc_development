################################################################################
#                                                                              #
#                                                                              #
# Propoortion testing                                                          #
#                                                                              #
#                                                                              #    
################################################################################

library(scater)
library(scran)
library(edgeR)
library(ggplot2)
library(cowplot)
library(reshape2)
library(cba)

find_hierarchical_order <- function(tmp, x=TRUE, y=TRUE){
  
  tmp <- acast(tmp, Cell_Type~Contrasts, value.var="FDR")
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

sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

abundances <- table(sce$major_clust, sce$batch) 
abundances <- unclass(abundances) 
head(abundances)

extra.info <- colData(sce)[match(colnames(abundances), sce$batch),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab
y.ab$samples$stage_ids <- factor(y.ab$samples$stage_ids, levels =
  c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))

keep <- filterByExpr(y.ab, group=y.ab$samples$stage_ids)
y.ab <- y.ab[keep,]
summary(keep)

y.ab$samples$Lot2 <- y.ab$samples$Library.Prep.Lot==2

design <- model.matrix(~0+stage_ids+Sex+chem+Lot2, y.ab$samples)
colnames(design) <- gsub("stage_ids", "", colnames(design))
y.contrasts <- makeContrasts(`Neonatal vs Fetal`=Neonatal-Fetal, 
    `Infancy vs Fetal`=Infancy-Fetal, `Childhood vs Fetal`=Childhood-Fetal,
    `Adolescence vs Fetal`=Adolescence-Fetal, `Adult vs Fetal`=Adult-Fetal, 
    `Infancy vs Neonatal`=Infancy-Neonatal, `Childhood vs Neonatal`=Childhood-Neonatal,
    `Adolescence vs Neonatal`=Adolescence-Neonatal, `Adult vs Neonatal`=Adult-Neonatal,
    `Childhood vs Infancy`=Childhood-Infancy, `Adolescence vs Infancy`=Adolescence-Infancy,
    `Adult vs Infancy`=Adult-Infancy, `Adolescence vs Childhood`=Adolescence-Childhood,
    `Adult vs Childhood`=Adult-Childhood, `Adult vs Adolescence`=Adult-Adolescence,
    levels=design)
y.contrasts <- makeContrasts(`Neonatal vs Fetal`=Neonatal-Fetal, 
    `Infancy vs Neonatal`=Infancy-Neonatal, 
    `Childhood vs Infancy`=Childhood-Infancy, 
    `Adolescence vs Childhood`=Adolescence-Childhood,
  `Adult vs Adolescence`=Adult-Adolescence,
    levels=design)

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)
fit <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit$var.prior)

decideTest <- lapply(colnames(y.contrasts), function(x){
  
  qlf <- glmQLFTest(fit, contrast=y.contrasts[,x])
  summary(decideTests(qlf))
})

decideTest <- do.call(cbind, decideTest)
decideTest

all_tests <- list()

for(i in 1:length(colnames(y.contrasts))){

  all_tests[[i]] <- 
    topTags(glmQLFTest(fit, contrast=y.contrasts[,colnames(y.contrasts)[i]]))
   all_tests[[i]] <- as.data.frame(all_tests[[i]])
   all_tests[[i]]$Contrasts <- colnames(y.contrasts)[i]
   all_tests[[i]]$Cell_Type <- rownames(all_tests[[i]])
    
}

all_tests <- do.call(rbind, all_tests)
all_tests
all_tests$FDR[all_tests$FDR>0.05] <- NA

order_xy <- find_hierarchical_order(all_tests)
all_tests$Cell_Type <- factor(all_tests$Cell_Type, levels=order_xy$y)
all_tests$Contrasts <- factor(all_tests$Contrasts, levels=colnames(y.contrasts))
all_tests$FDR <- all_tests$FDR*sign(all_tests$logFC)
all_tests$logFC[is.na(all_tests$FDR)] <- NA 

ggplot(all_tests, aes(x=Cell_Type, y=Contrasts, fill=logFC)) + geom_tile() + theme_classic() +
  theme(axis.title = element_blank(), 
  axis.text.x = element_text(angle = 45, hjust=1, size=8)) +
  scale_fill_gradientn(colours =colorspace::diverge_hcl(7) , 
  na.value = "transparent") + 
  ggtitle("Abundance testing") +
  background_grid()
