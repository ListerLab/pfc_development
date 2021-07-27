# Supplementary Figures

library(SingleCellExperiment)
library(Polychrome)
library(rcartocolor)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Figure unfiltered UMAP plot

save_objs <- readRDS("../UMAP_GeneScore_all.rds")
atac_samples <- read_xlsx("annotation/scATACseq_neuronal_maturation.xlsx")
atac_samples <- atac_samples[atac_samples$Used=="Yes",]


df <- save_objs[[3]]
df$Sample <- save_objs[[2]]$Sample
df$Sex <- atac_samples$Sex[as.numeric(match(df$Sample, atac_samples$Sample))]
df$Stage <- save_objs[[2]]$Stage
df$Cluster <- save_objs[[2]]$Clusters
df$DoubletScores <- save_objs[[2]]$DoubletScore
df$TSSEnrichment <- save_objs[[2]]$TSSEnrichment
colnames(df)[1:2] <- c("UMAP1", "UMAP2")
cell_types <- c(C1="Low quality",C2="Low quality",C3="Low quality",
    C4="Low quality", C5="Low quality", C6="Low quality",
    C7="Low quality", C8="Low quality", C9="Astro", C10="Low quality",
    C11="Astro", C12="Astro", C13="Astro",C14="Astro", C15= "Astro",
    C16="Low quality",C17="Astro", C18="Astro",C19="Astro",C20="Low quality",
    C21="Astro",C22="Low quality",C23="Astro",C24="Micro",
    C25="Micro",C26="Micro", C27="Micro",C28="Micro",C29="Low quality",
    C30="L4", C31="L4",
    C32="L4", C33="L4",C34="L2/3", C35="L4", C36="L2/3",
    C37="PN dev", C38="Low quality", C39="Low quality", C40="L5/6", 
    C41= "L2/3",
    C42="MGE der", C43="MGE der", C44="MGE der", C45="CGE der",
    C46="CGE der", C47="L5/6", C48="Low quality",
    C49="L2/3", C50="L2/3", C51="L2/3", C52="Vas",
    C53="CGE der", C54="IN dev", C55="IN dev", C56="IN dev",
    C57="Low quality", C58="PN dev", C59="PN dev", C60="PN dev", 
    C61="PN dev", C62="L5/6", C63="L5/6", C64="L5/6",
    C65="Oligo",  C66="Oligo", C67="Oligo", C68="Oligo", C69= "Oligo",
    C70="Oligo",C71="OPC", C72="OPC", C73="OPC",
    C74="OPC", C75="Low quality", C76="OPC", C77="OPC", C78="OPC",
    C79="OPC", C80="OPC")
names(cell_types) <- paste0("C", 1:length(cell_types))
df$Anno <- cell_types[match(df$Cluster, names(cell_types))]
df <- df[sample(1:nrow(df), nrow(df)),]
set.seed(10)

### Annotation

col <- c(Astro= '#ffc857', `IN dev`='#c6a785', 
         `L2/3`= '#6e0614', L4= '#8b3843', `L5/6`='#a86a72',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN dev`='#e2cdd0',
         `MGE der`='#b44622',`CGE der`='#1c5701', Vas='#a3a3a3',
         `Low quality`='lightgrey')

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Anno), size=0.285, alpha=0.5) +
    theme_classic() +  scale_color_manual(values=col) + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

ggsave(g1, file="supp_figures/SuppFig1_Annotation.png", height=7, width=7)


g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Sex), size=0.285, alpha=0.5) +
    theme_classic() + 
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))

### Cluster

p36 <- palette36.colors(36)
p36 <- rep(p36,3)
names(p36) <- NULL
p36 <- p36[1:length(unique(df$Cluster))]
names(p36) <- unique(df$Cluster)

indexes <- split(1:nrow(df), df$Cluster)
center <- t(sapply(indexes, function(x) c(mean(df[x,1]), mean(df[x,2]))))
colnames(center) <- c("UMAP1", "UMAP2")
center <- as.data.frame(center)
center$label <- rownames(center)

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Cluster), size=0.285, alpha=0.5) +
    theme_classic() + scale_colour_manual(values=p36) +
    guides(color=FALSE) +
    ggrepel::geom_label_repel(data=center, aes(x=UMAP1, y=UMAP2, label=label),
         colour="black", label.size = NA, fill=NA, max.overlaps = 100) 

ggsave(g1, file="supp_figures/SuppFig1_Clusters.png", height=7, width=7)

### TSS enrichment

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=TSSEnrichment), size=0.285, alpha=0.5) +
    theme_classic() + scale_color_viridis_c() 

ggsave(g1, file="supp_figures/SuppFig1_TSSEnrich.png", height=7, width=7)

### Doublet Score

g1 <- ggplot(df[order(df$DoubletScores, decreasing=FALSE),], aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=DoubletScores), size=0.285, alpha=0.5) +
    theme_classic() + scale_color_viridis_c() 

ggsave(g1, file="supp_figures/SuppFig1_DoubletScore.png", height=7, width=7)

## Figure filtered UMAP plot Supp 2

save_objs <- readRDS("../UMAP_GeneScore_rm.rds")

df <- save_objs[[3]]
df$Sample <- save_objs[[2]]$Sample
df$Sex <- atac_samples$Sex[as.numeric(match(df$Sample, atac_samples$Sample))]
df$Stage <- save_objs[[2]]$Stage
df$Cluster <- save_objs[[2]]$Clusters
df$Group <- save_objs[[2]]$predictedGroup
colnames(df)[1:2] <- c("UMAP1", "UMAP2")
df$Age <- save_objs[[2]]$arcsin_ages
df$Age <- as.numeric(as.character(df$Age))
cell_types <- c(C1="Micro",C2="Micro",C3="Micro",
    C4="Oligo", C5="Oligo", C6="Oligo",
    C7="Oligo", C8="Oligo", C9="Oligo", C10="MGE der",
    C11="IN dev", C12="IN dev", C13="L5/6",C14="CGE der", 
    C15= "MGE der", C16="MGE der", C17="PN dev", C18="PN dev",
    C19="L2/3", C20="L2/3", C21="L2/3",C22="L4",
    C23="L4",C24="L4",C25="L4",C26="L4",
    C27="L2/3", C28="L2/3", C29="L2/3", C30="PN dev", C31="L5/6",
    C32="L5/6", C33="L5/6",C34="PN dev",C35="PN dev", C36="L5/6",
    C37="PN dev", C38="L5/6", C39="L5/6", C40="IN dev", 
    C41= "MGE der", C42="CGE der", C43="CGE der", C44="OPC", 
    C45="Vas", C46="Astro",C47="Astro", C48="Astro",
    C49="Astro", C50="Astro", C51="Astro", C52="Astro",
    C53="Astro", C54="Astro", C55="Astro", C56="Astro",
    C57="OPC", C58="OPC", C59="OPC", C60="OPC", 
    C61="OPC", C62="OPC", C63="OPC", C64="OPC",
    C65="OPC")

names(cell_types) <- paste0("C", 1:length(cell_types))
df$Anno <- cell_types[match(df$Cluster, names(cell_types))]
df$L5_6 <- grepl("L5/6", df$Group)
df$L4 <- grepl("L4", df$Group)
df$L2_3 <- grepl("L2/3|L2|L3", df$Group)
df$CGE_der <- grepl("VIP|ID2", df$Group)
df$MGE_der <- grepl("SST|PV", df$Group)
df$IN_dev <- grepl("ID2_dev|MGE_dev-2|MGE_dev-1|CGE_dev|VIP_dev", df$Group)
df$Oligo <- grepl("Oligo", df$Group)
df$OPC <- grepl("OPC", df$Group)
df$Astro <- grepl("Astro", df$Group)
df$Vas <- grepl("Vas", df$Group)
df$Micro <- grepl("Micro", df$Group)
df <- df[sample(1:nrow(df), nrow(df)),]
set.seed(10)
df$Sample <- as.character(df$Sample)

### UMAP clusters

p36 <- palette36.colors(36)
p36 <- rep(p36,3)
names(p36) <- NULL
p36 <- p36[1:length(unique(df$Cluster))]
names(p36) <- unique(df$Cluster)

indexes <- split(1:nrow(df), df$Cluster)
center <- t(sapply(indexes, function(x) c(mean(df[x,1]), mean(df[x,2]))))
colnames(center) <- c("UMAP1", "UMAP2")
center <- as.data.frame(center)
center$label <- rownames(center)

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Cluster), size=0.285, alpha=0.5) +
    theme_classic() + scale_colour_manual(values=p36) +
    guides(color=FALSE) +
    ggrepel::geom_label_repel(data=center, aes(x=UMAP1, y=UMAP2, label=label),
    colour="black", label.size = NA, fill=NA, max.overlaps = 100) 

ggsave(g1, file="supp_figures/SuppFig2_Clusters.png", height=7, width=7)

### UMAP plot samples

my_colors = c(carto_pal(12, "Pastel"), carto_pal(5, "Antique"))

g1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) + 
    geom_point(col="black", size=0.8) +
    geom_point(col="white", size=0.3) +
    geom_point(aes(col=Sample), size=0.285, alpha=0.5) +
    theme_classic() + scale_colour_manual(values=my_colors) +
    guides(color= guide_legend(override.aes = list(size = 3, alpha=1)))


ggsave(g1, file="supp_figures/SuppFig2_Samples.png", height=7, width=7)

### Barplot samples

num_cell_sample <- data.frame(num_cells=table(df$Sample), 
                             sample=names(table(df$Sample)))

relevel_dat <- unique(save_objs[[2]][, c("Sample", "arcsin_ages")])
new_levels <- relevel_dat$Sample[order(relevel_dat$arcsin_ages, decreasing=FALSE)]

my_colors = c(carto_pal(12, "Pastel"), carto_pal(5, "Antique"))

num_cell_sample$sample <- factor(num_cell_sample$sample, 
    level=as.character(new_levels))

g1 <- ggplot(num_cell_sample, aes(x=sample, y=num_cells.Freq, fill=sample)) + 
    geom_col() +theme_classic() + scale_fill_manual(values=my_colors) + 
    guides(fill=FALSE) + ylab("Number of Cells") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g1, file="supp_figures/SuppFig2_CellNum_Sample.svg", 
       height=5.51, width=5.13, units="cm")

### Barplot stages

num_cell_stage <- 
    data.frame(num_cells=table(df$Stage), stage=names(table(df$Stage)))
num_cell_stage$stage <- factor(num_cell_stage$stage, level=c("Fetal", "Neonatal",
      "Infancy", "Childhood", "Adolescence", "Adult"))

col_stage <- c(Fetal="#512568", Neonatal="#443682", Infancy="#3D6B93",
               Childhood="#20988C", Adolescence="#98CA43", Adult="#F9E51B")

g1 <- ggplot(num_cell_stage, aes(x=stage, y=num_cells.Freq, fill=stage)) + 
    geom_col() +theme_classic() + scale_fill_manual(values=col_stage) + 
    guides(fill=FALSE) + ylab("Number of Cells") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(g1, file="supp_figures/SuppFig2_CellNum_Stage.svg", 
       height=5.51, width=5.13, units="cm")

### Composition detailed

df1 <- as.data.frame(save_objs[[2]])[,c("Stage","Clusters", "Sample")]
df1$Anno <- cell_types[match(df1$Cluster, names(cell_types))]

relevel_dat <- unique(save_objs[[2]][, c("Sample", "arcsin_ages")])
new_levels <- relevel_dat$Sample[order(relevel_dat$arcsin_ages, decreasing=FALSE)]

df1$Stage <- factor(df1$Stage, levels = c("Fetal", "Neonatal", "Infancy", 
       "Childhood", "Adolescence", "Adult"))
df1$Anno <- factor(df1$Anno, levels=rev(c('L2/3', 'L4', 
       'L5/6', 'PN dev', 'MGE der', 'CGE der', 'IN dev', 
        'Astro', 'Oligo', 'OPC', 'Micro', 'Vas')))
df1$Sample <- factor(df1$Sample, levels = as.character(new_levels))
g1 <- ggplot(data = df1, mapping = aes(y = Sample, fill = Anno)) + 
    geom_bar(position = "fill", colour="black") + scale_fill_manual(values = col) + 
    scale_y_discrete(limits = rev(levels(df1$Sample))) + theme_classic()

ggsave(g1, file="supp_figures/SuppFig2_CompSamples.svg", height=8, width=6)

### Annotation of Clusters

cell_types <- c("PN_dev", "L5_6", "L4", "L2_3", "CGE_der", "MGE_der", "IN_dev",
                "Oligo", "OPC", "Astro", "Vas", "Micro")

plot_lists <- list()

for(i in cell_types){
    plot_lists[[i]] <- ggplot(df, aes_string(x="UMAP1", y="UMAP2", color=i)) +
        geom_point(alpha=0.5, size=0.2) +
        theme_classic() + scale_color_viridis_d() + 
        guides(color = guide_legend(override.aes = list(size = 3, alpha=1)))
}

system(paste0("mkdir ", paths, "cluster_figures"))

for(i in cell_types){
    ggsave(plot_lists[[i]], file=paste0(paths, "supp_figures/cluster_figures/", i, ".png"))
}

## Peak quality

### UpSet plot

peak_grl <- readRDS("processed_data/cell_type_atac_peaks_filtered_anno_dge_stage_motif_gr_no_long.Rds")
peak_grl <- as(peak_grl, "GRangesList")
union_gr <- unlist(peak_grl) %>% reduce()

get_hits <- function(x){
    hits <- overlapsAny(query = union_gr, peak_grl[[x]])
    return(hits)
}

hits_df <- lapply(1:length(peak_grl), get_hits) %>% do.call(cbind, .)
hits_df <- data.frame(hits_df + 0)
colnames(hits_df) <- names(peak_grl)

UpSetR::upset(data = hits_df, order.by = "freq",
    decreasing = TRUE, nsets = 10, nintersects = 20)
#save by hand

### ChromHMM Fetal

dat <- read.table("processed_data/ChromHMM_FetalPeaks.txt", sep="\t")
dat$Tissue <- rownames(dat)
dat <- melt(dat, id="Tissue")
dat$Tissue <- factor(dat$Tissue, levels=rev(c("Brain Angular Gyrus", 
    "Brain Anterior Caudate", "Brain Cingulate Gyrus", "Brain Germinal Matrix",
    "Brain Hippocampus Middle", "Brain Inferior Temporal Lobe", 
    "Brain Dorsolateral Prefrontal Cortex", "Brain Substantia Nigra",
    "Fetal Brain Male", "Fetal Brain Female", "Fetal Lung", "Fetal Thymus",
    "Esophagus", "Left Ventricle", "Liver", "Lung", "Ovary", "Pancreas", 
    "Placenta", "Thymus",  "Spleen")))
blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5", "#08306B")

g1 <- ggplot(dat, aes(y=Tissue, x=variable, fill=value)) + 
    geom_tile() + theme_classic() + 
    scale_fill_gradientn(colours =
    grDevices::colorRampPalette(colors=blue_colours)(20)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) 

ggsave(g1, file="supp_figures/SuppFig3_ChromHMM_Fetal.svg", height=4, width=3.5)

dat <- read.table("processed_data/ChromHMM_Peaks_CellTypes.txt", sep="\t")
dat <- melt(dat, id=c("tissue", "cell_type"))
dat$tissue <- factor(dat$tissue, levels=rev(c("Brain Angular Gyrus", 
    "Brain Anterior Caudate", "Brain Cingulate Gyrus", "Brain Germinal Matrix",
    "Brain Hippocampus Middle", "Brain Inferior Temporal Lobe", 
    "Brain Dorsolateral Prefrontal Cortex", "Brain Substantia Nigra",
    "Fetal Brain Male", "Fetal Brain Female", "Fetal Lung", "Fetal Thymus",
    "Esophagus", "Left Ventricle", "Liver", "Lung", "Ovary", "Pancreas", 
    "Placenta", "Thymus",  "Spleen")))


g1 <- ggplot(dat, aes(y=tissue, x=variable, fill=value)) + 
    geom_tile() + theme_classic() + 
    scale_fill_gradientn(colours =
                             grDevices::colorRampPalette(colors=blue_colours)(20)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) + facet_grid(cols=vars(cell_type))

ggsave(g1, file="supp_figures/SuppFig3_ChromHMM_peaks_celltypes.svg", height=4, width=10)

### Annotation genome

processed_data <- readRDS(file=
  "processed_data/cell_type_atac_peaks_filtered_anno_dge_stage_gr.Rds")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
chromInfo <- seqinfo(txdb)
chromRanges <- GRanges(seqnames(chromInfo), 
                       IRanges(start = 1, end = seqlengths(chromInfo)))

# Create a GRanges object containing the merged exon ranges
exonRanges <- exons(txdb)
exonRanges <- reduce(exonRanges)

# Get the widths of the exon and chromosome ranges
exonWidths <- width(exonRanges)
chromWidths <- width(chromRanges)

# Calculate the total width of the exon and chromosome ranges
exonTotal <- sum(exonWidths)
chromTotal <- sum(chromWidths)

# Calculate the percentage of exonic nucleotides
pctExon <- (exonTotal / chromTotal) * 100

transcripts <- transcripts(txdb)

# Exon regions
exons <- exons(txdb)
strand(exons) <- "*"
exons <- reduce(exons)

# Intron regions
introns <- intronsByTranscript(txdb)
introns <- unlist(introns)
strand(introns) <- "*"
introns <- reduce(introns)

# 5 prime UTR
utr5 <- fiveUTRsByTranscript(txdb)
utr5 <- unlist(utr5)
strand(utr5) <- "*"
utr5 <- reduce(utr5)

# 3 prime UTR
utr3 <- threeUTRsByTranscript(txdb)
utr3 <- unlist(utr3)
strand(utr3) <- "*"
utr3 <- reduce(utr3)

# CDS regions
cds <- cds(txdb)
strand(cds) <- "*"
cds <- reduce(cds)

# Gene regions
genes <- genes(txdb)
strand(genes) <- "*"
genes <- reduce(genes)

# Promoter regions
promoters <- promoters(txdb, upstream = 500, downstream = 100)
promoters <- reduce(promoters)
strand(promoters) <- "*"

# Intergenic regions (Credit to Vince Buffalo and his book for this one)
tx <- reduce(transcripts)
strand(tx) <- "*"
intergenic <- GenomicRanges::setdiff(chromRanges, tx)
intergenic <- reduce(intergenic)

#collect everything in one list
all_annos <- list(exons=exons, introns=introns, utr5=utr5, utr3=utr3, cds=cds,
                  promoters=promoters, intergenic=intergenic)

enrich <- function(peaks, all_annos, chromTotal){
    
    res <- sapply(all_annos, function(y)
        c(observed=sum(width(GenomicRanges::intersect(peaks, y)))/sum(width(peaks))*100,
          expected=sum(width(y))/chromTotal*100))
    
    res <- melt(res)
    colnames(res) <- c("Group", "Region", "Percentage")
    res$width <- 0.8
    res$width[res$Group=="observed"] = 0.6
    res <- res[rev(order(res$Group)),]
    res$alpha <- 1
    res$alpha[res$Group=="expected"] = 0.6
    
    return(res)
}

enrich_plot <- function(res, title){
    
    g1 <- ggplot(res, aes(x=Percentage, y=Region, 
                          group=Group, fill=Group)) + 
        geom_bar(stat="identity",
                 position = "identity",width =res$width, alpha=res$alpha) +
        scale_fill_manual(values=c(expected="grey", observed="#17becf" )) +
        theme_classic() + ggtitle(title)
    
    return(g1)
    
}

res <- lapply(processed_data, function(x) 
    enrich(x, all_annos, chromTotal))

res_all <- res[[1]]
res_all[,3] <- apply(sapply(res, function(x) x[,3]), 1, mean)

g1 <- enrich_plot(res_all,"All peaks")

ggsave(g1, file="supp_figures/SuppFig4_AnnotationGenome.svg", height=5, 
       width=5.5)

## CRE element quality

### Correlation distribution

pearson_corr <- readRDS("processed_data/pearson_corr.RDS")
pearson_corr_background <- readRDS("processed_data/pearson_corr_background.RDS")

nn <- 0.025*nrow(pearson_corr_background)
thresh <- sort(abs(pearson_corr_background$Correlation), decreasing = T)[nn+1]

set.seed(4)
dat <- data.frame(Correlation=c(sample(pearson_corr$Correlation, 
    nrow(pearson_corr_background)), pearson_corr_background$Correlation),
    Group=c(rep("Real", nrow(pearson_corr_background)), rep("Background", 
    nrow(pearson_corr_background))))

g1 <- ggplot(dat, aes(x=Correlation, fill=Group)) + geom_density(alpha=0.4) +
    geom_vline(xintercept=thresh, linetype="dashed") + theme_classic() +
    annotate("text", label = paste0("threshold=", round(thresh,2)), 
    x = thresh+0.30, y = 4)

ggsave(g1, file="supp_figures/SuppFig4_CorrelationTresh.svg", height=4, width=6)


##Number
col <- c(Astro= '#ffc857',
         `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
         Micro='#484848',Oligo='#255f85', 
         `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')

pearson_corr_sig <- readRDS("processed_data/pearson_corr_sig.RDS")
cell_types <- sapply(pearson_corr_sig$peak_name, function(x) 
    strsplit(x, ".", fixed=T)[[1]][1])

df <- as.data.frame(table(cell_types))
df$cell_types <- factor(df$cell_types, levels=c("L2_3", "L4", "L5_6",
    "MGE_der", "CGE_der", "Astro", "Oligo", "Micro", "Vas"))
ggplot(df, aes(y=Freq, x=cell_types, fill=cell_types)) + geom_col() +
    theme_classic() + scale_fill_manual(values=col)

### ChromHMM 

dat <- read.table("processed_data/ChromHMM_CRE.txt", sep="\t")
dat$Tissue <- rownames(dat)
dat <- melt(dat, id="Tissue")
dat$Tissue <- factor(dat$Tissue, levels=rev(c("Brain Angular Gyrus", 
  "Brain Anterior Caudate", "Brain Cingulate Gyrus", "Brain Germinal Matrix",
  "Brain Hippocampus Middle", "Brain Inferior Temporal Lobe", 
  "Brain Dorsolateral Prefrontal Cortex", "Brain Substantia Nigra",
  "Fetal Brain Male", "Fetal Brain Female", "Fetal Lung", "Fetal Thymus",
  "Esophagus", "Left Ventricle", "Liver", "Lung", "Ovary", "Pancreas", 
   "Placenta", "Thymus",  "Spleen")))


blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(dat, aes(y=Tissue, x=variable, fill=value)) + 
    geom_tile() + theme_classic() + 
    scale_fill_gradientn(colours =
   grDevices::colorRampPalette(colors=blue_colours)(20)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) 

ggsave(g1, file="supp_figures/SuppFig4_ChromHMM_CRE.svg", height=4, width=3.5)


## TF for diseases

tf <- read.table("processed_data/tf_diseases.tsv", sep="\t")
tf <- tf[tf$padjust<0.05,]
tf$feature_short <- sapply(tf$feature, function(x) strsplit(x, "_")[[1]][1])

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")
a <- acast(tf,  feature_short ~ cell_type , mean, fill=1,
           drop=T, value.var="padjust")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
tf$feature_short <- factor(tf$feature_short,
   levels=levels_motifs)



g1 <- ggplot(tf, aes(x=feature_short, y=cell_type, fill=padjust)) + geom_tile() +
    theme_classic() + 
    scale_fill_gradientn(colours =
                             rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) 

ggsave(g1, file="supp_figures/TF_diseases.svg", height=5, width=7)
