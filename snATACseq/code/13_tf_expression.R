###############################################################################
#                                                                             #
#                                                                             #
# Find transcription factors not expressed                                    #
#                                                                             #
#                                                                             #    
###############################################################################
options(scipen=999)

library(JASPAR2018)
library(scater)

# 1. load TFs and scRNA-seq data

opts <- list()
opts[["species"]] <- 9606
hs_motifs <- TFBSTools::getMatrixSet(JASPAR2018, opts)

paths <- "snRNAseq/processed_data/"
sce <- readRDS(paste0(paths, 
  "/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS"))

# 2. get symbols to match

symbols <- sapply(hs_motifs, function(x) x@tags$symbol)
symbols[sapply(symbols, is.null)] <- NA
name_motif <- sapply(hs_motifs, function(x) x@name)
dat <- data.frame(motif=unlist(name_motif), symbol=unlist(symbols))
dat$symbol2 <- dat$motif
dat$symbol2 <- toupper(dat$symbol2)
dat$symbol2 <- gsub("(VAR.2)", "", dat$symbol2, fixed=T)
dat$symbol2 <- gsub("(VAR.3)", "", dat$symbol2, fixed=T)
dat$symbol2 <- sapply(dat$symbol2, function(x) strsplit(x, "::", fixed=T)[[1]])
dat$symbol2$`EWSR1-FLI1` <- c("EWSR1", "FLI1")

# 3. save TFs not expressed anywhere

dat$ind <- sapply(dat$symbol2, function(x) match(x, rowData(sce)$index))
tf_all <- names(dat$ind)[sapply(dat$ind, function(x) all(is.na(x)))]
saveRDS(tf_all, "snATACseq/processed_data/tf_not_rna.RDS")

# 4. identify TF expression in trajectories

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
     which(sce$`SST_TH_dev-traj`=="dev-traj")))

tf_cell_type <- lapply(input, function(x)
    rowSums(counts(sce[,x])))
tf_cell_type <- lapply(tf_cell_type, function(x)
    x[!x==0])

all_tfs_cell_type <- lapply(tf_cell_type, function(y)
  sapply(dat$symbol2, function(x) match(x, names(y))))
all_tfs_cell_type <- lapply(all_tfs_cell_type,
    function(x) names(x)[sapply(x, function(xx) all(is.na(xx)))])

# 5. save trajectory-specific TF expressions

saveRDS(all_tfs_cell_type, "snATACseq/processed_data/tf_not_cell_type_rna.RDS")
