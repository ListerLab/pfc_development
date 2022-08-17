################################################################################
#                                                                              #
#                                                                              #
# Quality control and UMAP plots for snATAC-seq data                           #
#                                                                              #
#                                                                              #    
################################################################################

library(ggplot2) 
library(SingleCellExperiment)
library(scuttle)
library(Polychrome)
library(readxl)
library(rcartocolor)

## 1. Quality plots prior TSS filtering

### 1a. Load data and relevel

cell_dat <- readRDS("snATACseq/processed_data/scatac_cell_data_prior_filtering.rds")
cell_dat<- as.data.frame(cell_dat)

relevel_dat <- unique(cell_dat[, c("Sample", "arcsin_ages")])
new_levels <- rev(relevel_dat$Sample[order(relevel_dat$arcsin_ages, 
          decreasing=FALSE)])

cell_dat$Sample <- factor(cell_dat$Sample, levels=new_levels)

my_colors = rev(c(carto_pal(12, "Pastel"), carto_pal(5, "Antique")))

### 1b. Calculate TSS enrichment outliers

tss_outliers <- list()

for(i in 1:length(new_levels)){
  
  sample_i <- new_levels[i]
  tss_enrich <- cell_dat$TSSEnrichment[cell_dat$Sample==sample_i]
  tss_outliers[[i]] <- isOutlier(tss_enrich, nmads=1,  type="lower")
  
}
names(tss_outliers) <- new_levels

lower_threshold <- data.frame(Sample=names(tss_outliers),
  Threshold=sapply(tss_outliers, function(x) attr(x, "thresholds")[1]))

lower_threshold$y <- c(0:(nrow(lower_threshold)-1))+0.525
lower_threshold$yend <- c(2:(nrow(lower_threshold)+1))-0.525

### 1c. Make TSS Enrichment boxplot prior filtering

g1 <- ggplot(cell_dat, aes(x=TSSEnrichment, y=Sample, fill=Sample)) + 
  geom_boxplot() +
  theme_classic() + scale_fill_manual(values=my_colors) + 
  geom_segment(data=lower_threshold, 
    aes(x = Threshold, xend = Threshold, y=y, yend=yend, linetype = "Cutoff")) +
  geom_segment(y=-650, yend = -650, x=40, xend=47, linetype = "dashed") +
  scale_linetype_manual("Cutoff", values=c("Cutoff"="dashed"))

ggsave(g1, file="supp_figures/tss_enrichment_boxplot.svg",
       height=4.88, width=4.2)

rm(list=ls())

## 2. UMAP plots after TSS filtering

### 2a. Read in data and relevel 

sce_atac <- readRDS("snATACseq/processed_data/umap_gene_score_all.rds")

relevel_dat <- unique(colData(sce_atac)[, c("Sample", "arcsin_ages")])
new_levels <- rev(relevel_dat$Sample[order(relevel_dat$arcsin_ages, 
          decreasing=FALSE)])
colData(sce_atac)$Sample <- as.character(colData(sce_atac)$Sample)
colData(sce_atac)$Sample <- factor(colData(sce_atac)$Sample, 
            levels=as.character(new_levels))

annotation_cluster <- read_xlsx("annotation/scatac_clustering_annotation.xlsx", 
                                sheet="first_clustering")
sce_atac$Anno <- annotation_cluster$Final_Annotation[match(sce_atac$Clusters,
              paste0("C", annotation_cluster$Cluster))]

umap <- reducedDim(sce_atac, "UMAP")
colnames(umap) <- c("UMAP1", "UMAP2")
umap <- cbind(umap, colData(sce_atac))
umap <- as.data.frame(umap)
set.seed(12)
umap <- umap[sample(1:nrow(umap)),]

col <- c(Astro= '#ffc857', `IN dev`='#c6d5c0', `L2/3`= '#6e0614', L4= '#8b3843', 
         `L5/6`='#a86a72', Micro='#484848',OPC='#92afc2', Oligo='#255f85', 
         `PN dev`='#e2cdd0', `Low quality`='lightgrey', `MGE der`='#b44622',
         `CGE der`='#1c5701', Vas='#a3a3a3')

### 2b. UMAP plot for annotation

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=Anno), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_manual(values=col) + 
  guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

ggsave(g1, file="supp_figures/umap_prior_filtering_anno.png", height=7, width=7)

### 2c. UMAP plot for doublet score

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=log10(DoubletScore+1)), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_viridis_c() 

ggsave(g1, file="supp_figures/umap_prior_filtering_doublet.png", height=7, 
       width=7)

### 2d. UMAP plot for promoter ratio

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=PromoterRatio), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_viridis_c() 

ggsave(g1, file="supp_figures/umap_prior_filtering_promoter.png", height=7, 
       width=7)

rm(list=ls())
gc()

## 3. Quality control and UMAP plot after final filtering

### 3a. Read in data and relevel

sce_atac <- readRDS("snATACseq/processed_data/umap_gene_score_rm.rds")

relevel_dat <- unique(colData(sce_atac)[, c("Sample", "arcsin_ages")])
new_levels <- rev(relevel_dat$Sample[order(relevel_dat$arcsin_ages, 
    decreasing=FALSE)])
colData(sce_atac)$Sample <- as.character(colData(sce_atac)$Sample)
colData(sce_atac)$Sample <- factor(colData(sce_atac)$Sample, 
                                   levels=as.character(new_levels))
colData(sce_atac)$Stage <- factor(colData(sce_atac)$Stage, levels=c(
  "Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult"))

annotation_cluster <- read_xlsx("annotation/scatac_clustering_annotation.xlsx", 
    sheet="second_clustering")
sce_atac$Anno <- annotation_cluster$Final_Annotation[match(sce_atac$Clusters,
   paste0("C", annotation_cluster$Cluster))]
sce_atac$Anno <- factor(sce_atac$Anno , levels=rev(c('L2/3', 'L4', 
    'L5/6', 'PN dev', 'MGE der', 'CGE der', 'IN dev', 
    'Astro', 'Oligo', 'OPC', 'Micro', 'Vas')))

umap <- reducedDim(sce_atac, "UMAP")
colnames(umap) <- c("UMAP1", "UMAP2")
umap <- cbind(umap, colData(sce_atac))
umap <- as.data.frame(umap)
set.seed(12)
umap <- umap[sample(1:nrow(umap)),]

col <- c(Astro= '#ffc857', `IN dev`='#c6d5c0', `L2/3`= '#6e0614', L4= '#8b3843', 
         `L5/6`='#a86a72', Micro='#484848',OPC='#92afc2', Oligo='#255f85', 
         `PN dev`='#e2cdd0', `Low quality`='lightgrey', `MGE der`='#b44622',
         `CGE der`='#1c5701', Vas='#a3a3a3')
my_colors = rev(c(carto_pal(12, "Pastel"), carto_pal(5, "Antique")))

### 3b. Make fragment boxplot after filtering

g1 <- ggplot(umap, aes(x=nFrags, y=Sample, fill=Sample)) + geom_boxplot() +
  theme_classic() + scale_fill_manual(values=my_colors) 

ggsave(g1, file="supp_figures/nfrags_boxplot.svg",
       height=4.88, width=4.2)

### 3c. Make promoter read boxplot after filtering

g1 <- ggplot(umap, aes(x=ReadsInPromoter, y=Sample, fill=Sample)) + 
  geom_boxplot() + theme_classic() + scale_fill_manual(values=my_colors) 

ggsave(g1, file="supp_figures/reads_promoters_boxplot.svg",
       height=4.88, width=4.2)

### 3d. Make promoter ratio boxplot after filtering

g1 <- ggplot(umap, aes(x=PromoterRatio, y=Sample, fill=Sample)) + 
  geom_boxplot() + theme_classic() + scale_fill_manual(values=my_colors) 

ggsave(g1, file="supp_figures/promoter_ratios_boxplot.svg",
       height=4.88, width=4.2)

### 3e. Make FRIP boxplot after filtering

g1 <- ggplot(umap, aes(x=FRIP, y=Sample, fill=Sample)) + 
  geom_boxplot() + theme_classic() + scale_fill_manual(values=my_colors) 

ggsave(g1, file="supp_figures/frip_boxplot.svg",
       height=4.88, width=4.2)

### 3f. Make read in peaks boxplot after filtering

g1 <- ggplot(umap, aes(x=ReadsInPeaks, y=Sample, fill=Sample)) + 
  geom_boxplot() + theme_classic() + scale_fill_manual(values=my_colors) 

ggsave(g1, file="supp_figures/read_peaks_boxplot.svg",
       height=4.88, width=4.2)

### 3g. Barplot for composition by stage

g1 <- ggplot(data = umap, mapping = aes(y = Stage, fill = Anno)) + 
  geom_bar(position = "fill", col="black", size=0.2) + 
  scale_fill_manual(values = col) + 
  scale_y_discrete(limits = rev(levels(umap$Stage))) + theme_classic() +
  guides(fill=FALSE)

ggsave(g1, file="main_figures/composition_stage_barplot.svg", 
       height=2.81*2, width=8.98*2, units="cm")

### 3h. Barplot for composition by sample

g1 <- ggplot(data = umap, mapping = aes(y = Sample, fill = Anno)) + 
  geom_bar(position = "fill", col="black", size=0.2) + 
  scale_fill_manual(values = col) + 
  theme_classic() +
  guides(fill=FALSE)

ggsave(g1, file="supp_figures/composition_sample_barplot.svg", 
       height=8.81*2, width=8.98*2, units="cm")

### 3i. UMAP of annotation

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=Anno), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_manual(values=col) + 
  guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

ggsave(g1, file="main_figures/umap_anno.png", height=7, width=7)

### 3j. UMAP of age

ages_dat <- data.frame(c("ga22","ga38","60d","1yr","10yr","20yr"),
  c(-0.59206997, -0.06903606,  0.29393111,  1.35923668,  3.59432467,
  4.28690546))
ages_dat <- ages_dat[order(as.numeric(ages_dat[,2])),]

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=arcsin_ages), size=0.285, alpha=0.5) +
  theme_classic() + scale_colour_gradientn(colours =
    grDevices::colorRampPalette(colors=c("#512568", "#443682", "#3D6B93",
    "#20988C", "#98CA43", "#F9E51B"))(20),
     breaks=as.numeric(ages_dat[,2]), labels=ages_dat[,1])

ggsave(g1, file="main_figures/umap_ages.png", height=5.3, width=4.9)

### 3k. UMAP by sample

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=Sample), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_manual(values=my_colors) + 
  guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

ggsave(g1, file="supp_figures/umap_sample.png", height=7, width=7)

### 3l. UMAP by clusters

p36 <- palette36.colors(36)
p36 <- rep(p36,3)
names(p36) <- NULL
p36 <- p36[1:length(unique(umap$Clusters))]
names(p36) <- unique(umap$Clusters)

g1 <- ggplot(umap, aes(x=UMAP1, y=UMAP2)) + 
  geom_point(col="black", size=0.8) +
  geom_point(col="white", size=0.3) +
  geom_point(aes(col=Clusters), size=0.285, alpha=0.5) +
  theme_classic() +  scale_color_manual(values=p36) + 
  guides(color= FALSE)

ggsave(g1, file="supp_figures/umap_clusters.png", height=7, width=7)

### 3m. Nuclei counts per sample

g1 <- ggplot(umap, aes(y=Sample, fill=Sample)) + 
  geom_bar() +theme_classic() + scale_fill_manual(values=my_colors) + 
  guides(fill=FALSE) + xlab("Number of Cells") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g1, file="supp_figures/nuclei_count_sample_barplot.svg", 
  height=7, width=7)

### 3n. Nuclei counts per stage

col_stage <- c(Fetal="#512568", Neonatal="#443682", Infancy="#3D6B93",
               Childhood="#20988C", Adolescence="#98CA43", Adult="#F9E51B")

g1 <- ggplot(umap, aes(x=Stage, fill=Stage)) + 
  geom_bar() +theme_classic() + scale_fill_manual(values=col_stage) + 
  guides(fill=FALSE) + ylab("Number of Cells") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g1, file="supp_figures/nuclei_count_stage_barplot.svg", 
       height=7, width=7)
 
### 3o. Annotation 

unique(umap$predictedGroup)

list_anno <- list(IN_dev = c("MGE_dev-1", "MGE_dev-2", "ID2_dev", "PV_dev", 
  "SST_CALB1_dev", "CGE_dev", "SST_ADGRG6_dev", "VIP_dev", "PV_SULF1_dev", 
  "PV_SCUBE3_dev","ID2_dev", "VIP_PCDH20_dev"),
PN_dev = c("L4_RORB_dev-fetal", "L5/6_TLE4_dev", "PN_dev", "L2_CUX2_LAMP5_dev",
  "L2/3_CUX2_dev-fetal", "L2/3_CUX2_dev-1", "L2/3_CUX2_dev-3", 
  "L5/6_THEMIS_dev-1", "L2/3_CUX2_dev-2", "L5/6_THEMIS_dev-2", 
  "L2/3_CUX2_dev-6", "L4_RORB_dev-2", "L2/3_CUX2_dev-4", "L2/3_CUX2_dev-5",
  "PV_dev"),
Micro = c("Micro_out", "Micro"),
Astro = c("Astro_SLC1A2", "Astro_GFAP", "Astro_dev-2", "Astro_dev-1", 
  "Astro_SLC1A2_dev", "Astro_dev-3"),
MGE_der = c("SST_STK32A", "SST_ADGRG6", "PV_WFDC2", "PV_SCUBE3", "SST_B3GAT2",
  "SST_CALB1", "PV_SST", "PV_SULF1", "SST_TH", "SST_BRINP3"),
CGE_der = c("VIP_ABI3BP", "CCK_MC4R", "LAMP5_NMBR", "LAMP5_LCP2", "VIP_HS3ST3A1",
  "VIP_CRH", "VIP_CHRM2", "LAMP5_CA1", "VIP-KIRREL3", "VIP_DPP6", "CCK_RELN",
  "ID2_CSMD1", "VIP_ADAMTSL1"),
Oligo = c("Oligo_mat", "Oligo-2", "Oligo-3", "Oligo-1", "Oligo-4", "Oligo-5"),
L2_3 = c("L3_CUX2_PRSS12", "L2_CUX2_LAMP5"),
L4 = c("L4_RORB_LRRK1", "L4_RORB_MET", "L4_RORB_MME"),
L5_6 = c("L5/6_TLE4_SCUBE1", "L5/6_TLE4_SORCS1", "L5/6_THEMIS_CNR1", 
         "L5/6_THEMIS_NTNG2", "L5/6_TLE4_HTR2C"),
Vas = c("Vas_TBX18", "Vas_PDGFRB", "Vas_CLDN5"),
OPC = c("OPC")) 
  
plot_lists <- list()

for(i in names(list_anno)){
  
  umap$tmp <- umap$predictedGroup %in% list_anno[[i]]
  
  plot_lists[[i]] <- ggplot(umap, aes_string(x="UMAP1", y="UMAP2", color="tmp")) +
    geom_point(alpha=0.5, size=0.2) +
    theme_classic() + scale_color_manual(values=c("black", "red")) + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
}

system(paste0("mkdir supp_figures/cluster_figures"))

for(i in names(list_anno)){
  ggsave(plot_lists[[i]], file=paste0("supp_figures/cluster_figures/", i, ".png"))
}

