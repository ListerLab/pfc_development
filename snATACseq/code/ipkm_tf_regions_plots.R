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

atac_samples <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]

sc <- loadArchRProject("snATACseq/clustering_final")

tf_diseases <-readRDS("snATACseq/processed_data/tf_diseases_gr.RDS")
tf_diseases <- unlist(tf_diseases, recursive = FALSE, use.names =TRUE)
diseases_motif <- fread("snATACseq/processed_data/tf_diseases_new.tsv", 
  sep="\t", data.table = FALSE)
diseases_motif <- diseases_motif[diseases_motif$adjusted_pvalue<0.05,]
diseases_motif$cell_type <- paste0(diseases_motif$trajectory, "..",
      diseases_motif$general_trend, "_", diseases_motif$pathway)
id <- diseases_motif$cell_type[1]
motif_interest <- diseases_motif$feature[1]

sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")
cell_types_conversions <- matrix(c("Astro", "Astro",
    "L2-3_CUX2", "L2/3",
    "L5-6_THEMIS", "L5/6",
    "L5-6_TLE4", "L5/6",
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

make_region_files <- function(id, motif_interest) {
    
    gr_interest <- tf_diseases[[grep(paste0(id, "$"), names(tf_diseases))]]
    i <- sapply(cell_types_conversions[,1], function(x) grepl(x, id))
    i <- cell_types_conversions[i,2]
    sc_sub <- sc[sc$Anno1 %in% i[1],]

    motifPositions <- getPositions(sc_sub)
    mm <- motifPositions[[motif_interest]]
    mm1 <- gr_interest[findOverlaps(gr_interest, mm)@from]
    mm <- mm[findOverlaps(gr_interest, mm)@to]
    
    
    mm_df <- data.frame(chr=seqnames(mm1), start=start(mm1),
        end=end(mm1), CRE=mm1$CRE, cell_type=strsplit(names(mm1)[1],
        ".", fixed=T)[[1]][1], start_motif=start(mm), end_motid=end(mm))
    mm_df <- mm_df[!duplicated(mm_df),]
    id <- gsub(",", "", id)
    id <- gsub("'", "", id)
    id <- gsub(" ", "_", id)
    mm_df$Symbol <- sapply(mm_df$CRE, function(x) {
      all_ens <- strsplit(x, "|", fixed=TRUE)[[1]]
      all_symbol <- rowData(sce)$index[match(all_ens, rowData(sce)$gene_ids)]
      all_symbol <- paste0(all_symbol, collapse = "|")
      all_symbol}
    )
    
    write.table(mm_df, file=paste0("snATACseq/processed_data/region_files/",
        id, "_", "regions_motifs", motif_interest, ".txt"), quote=FALSE, 
        sep="\t")
    
}

apply(diseases_motif, 1, function(x) make_region_files(x[16], x[1]))

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


fpkm_stage_plot <- function(id, motif_interest, col){
  
    gr_interest <- tf_diseases[[grep(paste0(id, "$"), names(tf_diseases))]]
    i <- sapply(cell_types_conversions[,1], function(x) grepl(x, id))
    i <- cell_types_conversions[i,2]
    sc_sub <- sc[sc$Anno1 %in% i[1],]

    motifPositions <- getPositions(sc_sub)
    mm <- motifPositions[[motif_interest]]
    mm <- gr_interest[findOverlaps(gr_interest, mm)@from]

    collect_fpkm <- collect_fpkm_fun(mm)

    collect_fpkm$variable <- gsub(paste0(gsub("/", "_", i), "_"), "",
                collect_fpkm$variable)
    collect_fpkm$stages <-atac_samples$Stage[
        match(collect_fpkm$variable, atac_samples$Sample)]
    collect_fpkm$stages <- factor(collect_fpkm$stages,
        levels=stages)

    collect_fpkm_mean <- collect_fpkm %>% 
        dplyr::group_by(stages) %>% 
        dplyr::summarise(value = mean(value))
    collect_fpkm_mean$peaks <- 1
    collect_fpkm_mean <- as.data.frame(collect_fpkm_mean)
    
    print(collect_fpkm_mean)
    
    id <- gsub(",", "", id)
    id <- gsub("'", "", id)
    id <- gsub(" ", "_", id)

    g1 <- ggplot(collect_fpkm_mean, aes(x=stages, y=value, group=peaks)) + 
        geom_line(size=3, color=col) +
        theme_classic() + scale_y_continuous(limits=c(0,1))

    ggsave(g1, file=paste0("snATACseq/region_figures/", id, "_", 
          motif_interest, ".pdf"))
    
}

fpkm_stage_plot("L5-6_THEMIS..interup_Bipolar Disorder",
  "SREBF1_73", col="#a86a72")
fpkm_stage_plot("L5-6_THEMIS..interup_Bipolar Disorder",
  "SREBF2_74", col="#a86a72")
fpkm_stage_plot("L4_RORB..up_Alzheimer's Disease",
  "FOSL1_40", col="#8b3843")
fpkm_stage_plot("L4_RORB..up_Schizophrenia",
  "FOS..JUN_418", col="#8b3843")
fpkm_stage_plot("L4_RORB..up_Schizophrenia",
  "FOSL1..JUNB_430", col="#8b3843")
fpkm_stage_plot("L2-3_CUX2..up_Schizophrenia",
  "BACH2_371", col="#6e0614")
fpkm_stage_plot("L5-6_TLE4..up_Bipolar Disorder",
  "FOS..JUND_434", col="#c59ba1")


opts <- list()
opts[["species"]] <- 9606
hs_motifs <- TFBSTools::getMatrixSet(JASPAR2018, opts)


svg("main_figures/FOSL1.svg", height=2, width=4)
ind <- grep("FOSL1$", sapply(hs_motifs, function(x) x@name))
seqLogo(apply(hs_motifs[[ind]]@profileMatrix, 
    2, function(x) x/sum(x)), ic.scale=FALSE, xaxis=FALSE, yaxis=FALSE)
dev.off()

