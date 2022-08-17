library(ArchR)
library(SummarizedExperiment)
library(ggplot.multistats)
library(readxl)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
library(seqLogo)
library(JASPAR2018)

addArchRGenome("hg19")
addArchRThreads(threads = 10)

atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("clustering_annotation")

tf_diseases <-readRDS("processed_data/tf_diseases_gr.RDS")
diseases_motif <- read.table("processed_data/tf_diseases.tsv", sep="\t")
diseases_motif <- diseases_motif[diseases_motif$padjust<0.05,]
id <- diseases_motif$cell_type[1]
motif_interest <- diseases_motif$feature[1]

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm=T))}

collect_fpkm_fun <- function(X) {
    aa <- X@elementMetadata
    aa <- aa[,grepl("_RL", colnames(aa))]
    aa <- as.matrix(aa)
    aa <- t(apply(aa, 1, function(x) range01(x)))
    aa <- as.data.frame(aa)
    aa$peaks <- names(X)
    aa <- melt(aa, id="peaks")
    return(aa)
}

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")
cell_types_conversions <- matrix(c("Astro", "Astro",
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
    "OPC", "OPC",
    "Micro", "Micro"
), ncol=2, byrow=T)

fpkm_stage_plot <- function(id, motif_interest, path, col) {
    
    gr_interest <- tf_diseases[[id]]
    i <- sapply(cell_types_conversion[,1], function(x) grepl(x, id))
    i <- cell_types_conversion[i,2]
    sc_sub <- sc[sc$Anno1 %in% i[1],]

    motifPositions <- getPositions(sc_sub)
    mm <- motifPositions[[motif_interest]]
    mm <- gr_interest[findOverlaps(gr_interest, mm)@from]
    
    mm_df <- data.frame(chr=seqnames(unique(mm)), start=start(unique(mm)),
        end=end(unique(mm)), CRE=unique(mm)$CRE, cell_type=strsplit(names(mm)[1],
        ".", fixed=T)[[1]][1])
    write.table(mm_df, file=paste0(id, "_", "regions_motifs", 
        motif_interest, ".txt"), 
        quote=FALSE, sep="\t")

    collect_fpkm <- collect_fpkm_fun(mm)

    collect_fpkm$variable <- gsub(paste0(i, "_"), "",
                collect_fpkm$variable)
    collect_fpkm$stages <-atac_samples$Stage[
        match(collect_fpkm$variable, atac_samples$Sample)]
    collect_fpkm$stages <- factor(collect_fpkm$stages,
        levels=stages)

    collect_fpkm_mean <- collect_fpkm %>% 
        dplyr::group_by(stages) %>% 
        dplyr::summarise(value = mean(value))
    collect_fpkm_mean$peaks <- 1

    collect_fpkm<- collect_fpkm %>% 
        dplyr::group_by(stages, peaks) %>% 
        dplyr::summarise(value = mean(value))
    collect_fpkm_mean$peaks <- 1
    
    print(collect_fpkm_mean)
    

        g1 <- ggplot(collect_fpkm, aes(x=stages, y=value,group=peaks)) + 
            geom_line(size = 0.2, alpha=0.2) +
            geom_line(data = collect_fpkm_mean, alpha = 1, size = 3) +
            theme_classic() 

        ggsave(g1, file=paste0(path, "/Fig5_", id, "_", motif_interest, ".pdf")
    
}

fpkm_stage_plot("Oligo_4_Autism Spectrum Disorders", "SOX13_406",
    "paper_figures/", col='#255f85')
fpkm_stage_plot("L4_RORB_6_Bipolar Disorder", "FOS..JUND_434",
 "paper_figures/", col='#a86a72')

opts <- list()
opts[["species"]] <- 9606
hs_motifs <- TFBSTools::getMatrixSet(JASPAR2018, opts)

svg("paper_figures/Fig5_SOX13.svg", height=2, width=4)
ind <- grep("SOX13$", sapply(hs_motifs, function(x) x@name))
seqLogo(apply(hs_motifs[[ind]]@profileMatrix, 
    2, function(x) x/sum(x)), ic.scale=FALSE, xaxis=FALSE, yaxis=FALSE)
dev.off()

svg("paper+figures/Fig5_FOS..JUND.svg", height=2, width=4)
ind <- grep("FOS..JUND$", sapply(hs_motifs, function(x) x@name))
seqLogo(apply(hs_motifs[[ind]]@profileMatrix, 
    2, function(x) x/sum(x)), ic.scale=FALSE, xaxis=FALSE, yaxis=FALSE)
dev.off()


paths <-  ""
sce <- readRDS(paste0(paths, 
     "2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))
all_file <- list.files("region_files")

for(i in 1:length(all_file)){

df <- read.table(paste0("region_files/", all_file[i]), sep="\t")
df$position <- paste0(df$chr, ":", df$start, "-", df$end)

df$CRE_new <- sapply(df$CRE, function(x) 
    strsplit(x, "|", fixed=T)[[1]])
df$CRE_new <- sapply(df$CRE_new, function(x) 
    rowData(sce)$index[match(x, rowData(sce)$gene_ids)])
df$CRE_new <- sapply(df$CRE_new, function(x)
    paste0(x, collapse="|"))
write.table(df, file=paste0("region_files/", all_file[i]), sep="\t",
      quote=F)

}

