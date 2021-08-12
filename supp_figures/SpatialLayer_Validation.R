library(spatialLIBD)
library(ggplot2)
library(dplyr)

# 1. Load spatial data
sce <- fetch_data(type = 'sce')

# 2. Make pictures for PCDH15

dat <- data.frame(Layer=sce$layer_guess_reordered_short,
    PCDH15=sce@assays@data$logcounts[which(rowData(sce)$gene_name=="PCDH15"),],
    Sample=sce$sample_name)

dat <- dat %>% group_by(Sample, Layer) %>% summarize(PCDH15=mean(PCDH15)) 
dat <- dat %>% filter(!is.na(Layer))

g1 <- ggplot(dat, aes(x=Layer, y=PCDH15)) + geom_boxplot() + 
    geom_point() + theme_classic() 
ggsave(g1, file=paste0("supp_figures/PCDH15_spatial.svg"),
       height=4, width=7)

# 3. Make pictures for HTR2A

dat <- data.frame(Layer=sce$layer_guess_reordered_short,
    HTR2A=sce@assays@data$logcounts[which(rowData(sce)$gene_name=="HTR2A"),],
                  Sample=sce$sample_name)

dat <- dat %>% group_by(Sample, Layer) %>% summarize(HTR2A=mean(HTR2A)) 
dat <- dat %>% filter(!is.na(Layer))

g1 <- ggplot(dat, aes(x=Layer, y=HTR2A)) + geom_boxplot() + 
    geom_point() + theme_classic() 
ggsave(g1, file=paste0("supp_figures/HTR2A_spatial.svg"),
       height=4, width=7)

# 4. Diseases genes

paths <- "snRNAseq/processed_data/disease_enrichment_genes/"
combs <- c("L5-6_TLE4_13_Schizophrenia", "L5-6_TLE4_6_Schizophrenia",
           "L5-6_THEMIS_6_Schizophrenia", 
           "L2-3_CUX2_13_Autism Spectrum Disorders")


for (i in 1:length(combs)){

    disease_genes <- read.csv(paste0(paths, combs[i], ".csv"), header=TRUE)
    disease_genes <- disease_genes[,2]
    
    
    tmp <- colMeans(as.matrix(sce@assays@data$logcounts[rowData(sce)$gene_name %in% 
                                  disease_genes,]))
    
    dat <- data.frame(Layer=sce$layer_guess_reordered_short,
    Genes=tmp, Sample=sce$sample_name)
    dat <- dat %>% group_by(Sample, Layer) %>% summarize(Genes=mean(Genes)) 
    dat <- dat %>% filter(!is.na(Layer))

    g1 <- ggplot(dat, aes(x=Layer, y=Genes)) + geom_boxplot() + 
    geom_point() + theme_classic() + ggtitle(combs[i]) 
    ggsave(g1, file=paste0("supp_figures/", combs[i], "_spatial.svg"),
           height=4, width=7)
    
}

