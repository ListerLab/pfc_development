################################################################################
#                                                                              #
#                                                                              #
# Plots concerning peaks and cis regulatory elements                           #
#                                                                              #
#                                                                              #    
################################################################################

library(ggplot2) 
library(reshape2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(readxl)
library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(parallel)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(cba)
library(org.Hs.eg.db)
library(clusterProfiler)
library(viridis)

## 1. ChromHMM comparisons

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

### 1a. ChromHMM peaks across all cell types 

dat <- read.table("snATACseq/processed_data/chrom_hmm_all_peaks.txt", sep="\t")
dat$tissue <- rownames(dat)
dat <- melt(dat, id="tissue")
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
          axis.title = element_blank()) 
ggsave(g1, file="main_figures/chrom_hmm_all_peaks.svg", width=3.5, 
       height=4)


### 1b. ChromHMM peaks for each cell type

dat <- read.table("snATACseq/processed_data/chrom_hmm_peaks_cell_types.txt", 
      sep="\t")
dat <- melt(dat, id=c("tissue", "cell_type"))
dat$tissue <- factor(dat$tissue, levels=rev(c("Brain Angular Gyrus", 
    "Brain Anterior Caudate", "Brain Cingulate Gyrus", "Brain Germinal Matrix",
    "Brain Hippocampus Middle", "Brain Inferior Temporal Lobe", 
    "Brain Dorsolateral Prefrontal Cortex", "Brain Substantia Nigra",
    "Fetal Brain Male", "Fetal Brain Female", "Fetal Lung", "Fetal Thymus",
    "Esophagus", "Left Ventricle", "Liver", "Lung", "Ovary", "Pancreas", 
    "Placenta", "Thymus",  "Spleen")))
dat$cell_type <- factor(dat$cell_type, levels=c("L2_3", "L4",
    "L5_6", "CGE_der", "MGE_der", "Astro", "Oligo", "Micro", "Vas"))


g1 <- ggplot(dat, aes(y=tissue, x=variable, fill=value)) + 
    geom_tile() + theme_classic() + 
    scale_fill_gradientn(colours =
        grDevices::colorRampPalette(colors=blue_colours)(20)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title = element_blank()) + facet_grid(cols=vars(cell_type))

ggsave(g1, file="supp_figures/chrom_hmm_all_cell_types.svg", 
       height=4, width=10)

## 2. Peak information

### 2a. Genomic locations of peaks

peaks <- readRDS(file=
  "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.Rds")

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
    
    res <- reshape2::melt(res)
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

res <- lapply(peaks, function(x) 
    enrich(x, all_annos, chromTotal))

res_all <- res[[1]]
res_all[,3] <- apply(sapply(res, function(x) x[,3]), 1, mean)

g1 <- enrich_plot(res_all,"All peaks")

ggsave(g1, file="supp_figures/genomic_locations_anno.svg", height=5, 
       width=5.5)

### 2b. Peak width

peaks <- readRDS(file=
   "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.Rds")

all_peaks <- as(peaks, "GRangesList")
all_peaks <- unlist(all_peaks, recursive = TRUE)

dat <- data.frame(id=names(all_peaks), width=width(all_peaks),
                   cell_type=sapply(names(all_peaks), function(x) 
                       strsplit(x, ".", fixed=T)[[1]][1]))

g1 <- ggplot(dat, aes(x=width)) + geom_density(fill="black") + 
    theme_classic() + ggtitle("All peaks") + scale_x_log10()
ggsave(g1, file="supp_figures/peak_width.svg", width=3.5, height=2)


### 2c. PC components

peak_mat_list <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.Rds")
mdat <- readxl::read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
mdat <- mdat[mdat$Used=="Yes",]

col <- c(Astro= '#ffc857', 
         `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
         Micro='#484848', Oligo='#255f85', 
         `MGE_der`='#b44622',`CGE_der`='#1c5701')

plot_pca <- function(cell_types, fpkm=TRUE){
    
    cell_dat <- lapply(peak_mat_list[cell_types], function(x)
        x@elementMetadata)
    
    # Extract data
    cell_dat <- lapply(cell_dat, function(x)
        x[,grepl(paste0("_RL"), colnames(x))])
    cell_dat <- lapply(cell_dat, function(x)  x%>% as.data.frame() %>%
        limma::normalizeBetweenArrays(method = "quantile"))
    for(i in 1:length(cell_dat)){
        
        colnames(cell_dat[[i]]) <- gsub(paste0(names(cell_dat)[i], "_"), "", 
            colnames(cell_dat[[i]]))
    }
    
    pr <- lapply(cell_dat, function(x) prcomp(t(x), scale.=TRUE))
    pc1 <- lapply(pr, function(x) 
        (summary(x)$importance[2, 1] * 100) %>% round(digits = 2))
    tmp <- unlist(lapply(1:length(cell_types), function(xx) rep(cell_types[xx],
        length=ncol(cell_dat[[xx]]))))
    
    pca_df <- data.frame(PC1=unlist(lapply(pr, function(xx) xx$x[ ,1])),
        Cell_type=tmp,
        Sample=unlist(lapply(cell_dat, function(x) colnames(x))))

    ind <- match(pca_df$Sample, mdat$Sample)
    stage <- mdat$Stage[ind] 
    arc_age <- mdat[ind,2]
    
    pca_df$Stage <- stage
    pca_df$Age <- arc_age
    pca_df$Stage <- factor(pca_df$Stage, levels=c("Fetal", "Neonatal",
        "Infancy", "Childhood", "Adolescence", "Adult"))
    
    pca_df_mean <- pca_df %>% group_by(Stage, Cell_type) %>% 
        summarise(PC1=mean(PC1))
    
    return(list(mean=pca_df_mean, dat=pca_df))
    
}


df_pc <- plot_pca(c("L2_3", "L4", "L5_6", "CGE_der", "MGE_der",
                    "Astro", "Oligo", "Micro"))

gg_pca <- ggplot(data = df_pc$dat, mapping = aes(x = Stage, y = PC1,
        col=Cell_type)) + geom_point() +
        geom_line(dat=df_pc$mean, aes(x = Stage, y = PC1,
           group=Cell_type), alpha = 1, size = 2) +
        ggtitle("PC1 component per stage") + xlab("") +
        ylab(str_c("")) + theme_classic() +
        scale_color_manual(values=col) + theme(axis.title.y = element_text(
            size=10))

ggsave(gg_pca, file="main_figures/pca_atac.svg", width=6, height=5)


## 3. Information on CREs 

### 3a. Correlation threshold

### 3b. Number of CREs

pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")

peaks_per_gene <- split(pearson_corr_sig$peak_name, pearson_corr_sig$Gene)

find_number_unique_peaks <- function(x, pearson_corr_sig) {
    
    ind <- match(x, pearson_corr_sig$peak_name)
    gr <- GRanges(seqnames = pearson_corr_sig$seqnames[ind],
                  IRanges(start=pearson_corr_sig$start[ind], 
                          end=pearson_corr_sig$end[ind]))
    gr <- GenomicRanges::reduce(gr)
    return(length(gr))
    
}

no_peaks_genes <- mclapply(peaks_per_gene, function(x) 
  find_number_unique_peaks(x, pearson_corr_sig), mc.cores=3)

median(unlist(no_peaks_genes))

hist(unlist(no_peaks_genes), breaks=30)
dat <- data.frame(Freq=table(unlist(no_peaks_genes)))
colnames(dat) <- c("No CREs", "Freq")
dat$`No CREs` <- as.numeric(dat$`No CREs`)
gg <- ggplot(dat, aes(x=`No CREs`, y=Freq)) + geom_col() + theme_classic()

ggsave(gg, file="supp_figures/hist_no_cres.svg",
       height=7, width=7)

### 3c. Distance to TSS

pearson_corr_sig <- readRDS("snATACseq/processed_data/pearson_corr_sig.RDS")


g1 <- ggplot(pearson_corr_sig , aes(x=dist+1)) + geom_density(fill="black") + 
    theme_classic() + scale_x_log10(breaks=c(1,50,1000,10000, 250000),
    labels=c(0,50,1000,10000, 250000)) + xlab("Distance to TSS (bp)") +
    theme(axis.title.y = element_blank())

df <- pearson_corr_sig %>% select("dist", "peak_name")
df$peak_name <- sapply(df$peak_name, function(x) 
    strsplit(x, ".", fixed=T)[[1]][1])
df <- df %>% mutate(type=case_when(dist==0 ~ "Geneic",
                                   (dist>0 & dist<=2000) ~ "Proximal",
                                   dist>1000 ~ "Distal"))

df1 <- df %>% group_by(type) %>% summarise(n=n())
df1$fraction = df1$n/ sum(df1$n)

# Make donut plot
g2 <- ggdonutchart(data=df1, "fraction", label="type",
                   color = "white", fill="type",
                   palette = brewer.pal(3, "Set3")) + guides(fill=FALSE)

gg <- ggdraw(g1) + 
    draw_plot(g2, 0.2, .25, .6, .6) 

ggsave(gg, file="main_figures/cre_distance_tss.svg", width=2.5*3.5, 
       height=2.5*2)


### 3d. Average normalized IPKM across stages

atac_meta <- read_xlsx("annotation/scatac_neuronal_maturation.xlsx")
atac_meta <- atac_meta[atac_meta$Used=="Yes",]

range01 <- function(x){(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm=T))}

paths <- "scrna/processed_data/"

gene_trends <- read.csv(paste0(paths, "gene_cluster_ids.csv"))
sce <- readRDS(paste0(paths, 
    "2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS"))

peaks <- readRDS(
    "snATACseq/processed_data/cell_type_atac_peaks_filtered_anno_gr.rds")
cis <- lapply(peaks, function(x) x[!is.na(x$CRE)])
cis <- as(cis, "GRangesList")

cis <- cis[!names(cis)=="Vas"]

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

gene_trends$gene_name <- rowData(sce)$gene_ids[match(gene_trends$gene_name, 
    rownames(sce))]
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
                       "OPC", "OPC",
                       "Micro", "Micro"
), ncol=2, byrow=T)

gene_trends <- lapply(names(cis),  function(x)
    gene_trends[gene_trends$major_clust %in%
                    cell_types[cell_types[,2]==x,1],])
names(gene_trends) <- names(cis)

for(i in 1:length(cis)){
    
    ind <- sapply(cis[[i]]$CRE, function(x) strsplit(x, "|", fixed=T)[[1]])
    ind <- sapply(ind, function(x) any(x %in% gene_trends[[i]]$gene_name))
    cis[[i]] <- cis[[i]][ind]
}


collect_fpkm <- lapply(cis, function(x) 
    collect_fpkm_fun(x))

for(i in 1:length(collect_fpkm)){
    
    collect_fpkm[[i]]$cell_type <- names(collect_fpkm)[i]
    collect_fpkm[[i]]$variable <- gsub(paste0(names(collect_fpkm)[i], "_"), "",
          collect_fpkm[[i]]$variable)
    collect_fpkm[[i]]$stages <-atac_meta$Stage[
        match(collect_fpkm[[i]]$variable, atac_meta$Sample)]
}

#get associated genes and gene trends

attach_gene_trends <- function(tmp_fpkm, cis_tmp, gene_trends_tmp){
    
    cre_genes <- cis_tmp$CRE[match(tmp_fpkm$peaks, 
            names(cis_tmp))]
    all_trends <- lapply(cre_genes, function(a)
        gene_trends_tmp[grep(a, gene_trends_tmp[,3]),4])
    ind <- lapply(1:length(all_trends), function(x)
        rep(x, length(all_trends[[x]])))
    ind <- unlist(ind)
    tmp_fpkm <- tmp_fpkm[ind,]
    tmp_fpkm$all_trends <- unlist(all_trends)
    return(tmp_fpkm)
    
}

collect_fpkm_new <- mclapply(1:length(cis), function(i) 
    attach_gene_trends(collect_fpkm[[i]], 
        cis[[i]], gene_trends[[i]]), mc.cores=3)


collect_fpkm_new <- do.call(rbind, collect_fpkm_new)
collect_fpkm_new <- collect_fpkm_new[!is.null(unlist
    (collect_fpkm_new$all_trends)),]

stages <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")

collect_fpkm_new$stages <- factor(collect_fpkm_new$stages, levels=stages)

cols <- c(down="#02577F", `interdown`="#D9D9D7",
          `interup`="#FFD433", up="#F3704A")

collect_fpkm_new$stages_num <- as.numeric(collect_fpkm_new$stages)

dlabs <- collect_fpkm_new %>%
    group_by(all_trends) %>% 
    arrange(stages_num) %>% 
    filter(stages_num %in% c(first(stages_num), last(stages_num)))

collect_fpkm_new$major_trends <-  gene_trends$Astro$major_trend[
  match(collect_fpkm_new$all_trends, gene_trends$Astro$gene_trend)]

g1 <- ggplot(collect_fpkm_new, aes(x=stages_num, y=value, 
    group=all_trends, col=major_trends)) +
    stat_summary(fun.y=mean, geom="line", size = 1) +
    theme_classic() + 
    stat_summary(data = dlabs, aes(label = all_trends, 
        hjust = ifelse(stages_num == first(stages_num), 1.2, -.2)), 
        fun = mean, geom = "text", color = "black")  +
    theme(axis.title.x = element_blank()) + ylab("Normalized fpkm") +
    scale_color_manual(values=cols) + scale_x_continuous(labels=stages,
        breaks=1:6,expand = expansion(mult = .1))

ggsave(g1, file="main_figures/scatac_gene_trends.svg", height=4, width=4)
    
### 3e. Plot NMF of CREs

source("snATACseq/R/functions_nmf.R")

plot_list_old <- readRDS("snATACseq/processed_data/plot_atac_rna_corr.RDS")
info_type <- plot_list_old$info$type
info_stages_rna = plot_list_old$info$stages_rna

H_class <- readRDS("snATACseq/processed_data/nmf_sam_H_class.RDS")
tab_cluster_H<- table(info_type, H_class$class0)
tab_stage_H <- table(info_stages_rna, H_class$class0)

# Pie charts

cluster_of_interest <- unique(H_class$class0)

plot_pie <- function(x){
    
    dat_H <- data.frame(cluster=rownames(tab_cluster_H),
                        value=tab_cluster_H[,x])
    dat_stage <- data.frame(stage=rownames(tab_stage_H),
                        value=tab_stage_H[,x])
    
    col_stages_rna <- viridis(6, direction=1)
    names(col_stages_rna) <- c("Fetal", "Neonatal", "Infancy", "Childhood", 
            "Adolescence", "Adult")
    col_types <- c(Astro= '#ffc857', `IN_dev`='#c6a785', 
       `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
        Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN_dev`='#e2cdd0',
       `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')
    
    g1 <- ggplot(dat_stage, aes(x="", y=value, fill=stage))  + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=col_stages_rna) + theme_void() +
        guides(fill=FALSE)
    ggsave(g1, file=paste0("main_figures/pie_charts/", x, "_pie_chart_stage.svg"),
           height=2, width=2)
    g2 <- ggplot(dat_H, aes(x="", y=value, fill=cluster))  + 
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0) +
        scale_fill_manual(values=col_types) + theme_void() +
        guides(fill=FALSE)
    ggsave(g2, file=paste0("main_figures/pie_charts/", x, "_pie_chart_cluster.svg"),
           height=2, width=2)
    
}

lapply(cluster_of_interest, function(x) plot_pie(x))

# Heatmap

plot_list <- readRDS("snATACseq/processed_data/nmf_sam_atac_rna.RDS")
plot_list_atac <- plot_list$sam_atac_new
plot_list_rna <- plot_list$rna_seq

ind_classes <- split(1:nrow(plot_list_rna), plot_list$classes_new)
plot_list_rna_classes <- t(sapply(ind_classes, function(x)
    colSums(plot_list_rna[x,, drop=FALSE])))
rownames(plot_list_rna_classes) <- names(ind_classes)
hc <- hclust(dist(plot_list_rna_classes))
order_row <- order.optimal(dist(plot_list_rna_classes), hc$merge)
classes_order <- rownames(plot_list_rna_classes)[order_row$order]
plot_list$classes_new <- factor(plot_list$classes_new, levels=classes_order)

my_hclust_row <- order(plot_list$classes_new)
dat_rna <- t(apply(plot_list_rna, 1, function(x) scale(x)))
dat_atac <- t(apply(plot_list_atac,1, function(x) scale(x)))

hc <- hclust(dist(t(dat_rna[my_hclust_row,])))
order_col <- order.optimal(dist(t(dat_rna[my_hclust_row,])), hc$merge)
my_hclust_cell <- order_col$order


dat_rna <- dat_rna[my_hclust_row, my_hclust_cell]
dat_atac <- dat_atac[my_hclust_row, my_hclust_cell]

col_stages_rna <- viridis(6, direction=1)
names(col_stages_rna) <- c("Fetal", "Neonatal", "Infancy", "Childhood", 
                           "Adolescence", "Adult")
col_types <- c(Astro= '#ffc857', `IN_dev`='#c6a785', 
      `L2_3`= '#6e0614', L4= '#8b3843', `L5_6`='#a86a72',
       Micro='#484848',OPC='#92afc2', Oligo='#255f85', `PN_dev`='#e2cdd0',
      `MGE_der`='#b44622',`CGE_der`='#1c5701', Vas='#a3a3a3')

ha_left = HeatmapAnnotation(cell = info_type[my_hclust_cell],
                            stages = factor(info_stages_rna[my_hclust_cell],
                                            level=names(col_stages_rna)),
                            annotation_name_side = "left",
                            col = list(stages = col_stages_rna,
                                       cell=col_types))
ha_right = HeatmapAnnotation(cell = info_type[my_hclust_cell],
                             stages = factor(info_stages_rna[my_hclust_cell],
                                             level=names(col_stages_rna)),
                             annotation_name_side = "right",
                             col = list(stages = col_stages_rna,
                                        cell=col_types))

cluster_col<- c('#155145', '#d3d3d3', '#FE1886', '#BDC90F', '#4903FB', '#26FFA8',
    '#A184C0', '#910716', '#06E802', '#0798F1', '#68B167', '#D433FC',
    '#D96F4E', '#FDF853', '#602683', '#72D6F5', '#577A05', '#FDB2AB',
    '#3274A3', '#A8EE91', '#78FA2B', '#7254EB', '#FC3314', '#1C2DCA',
    '#970FCC')[1:length(levels(plot_list$classes_new))]
names(cluster_col) <- levels(plot_list$classes_new)

png("main_figures/cluster_nmf_legend.png",  width=10, height=12, unit="cm",
    res=1200)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =levels(plot_list$classes_new)[1:15], 
       pch=16, pt.cex=1, cex=1, bty='n',
    col = cluster_col[1:15])
legend("topright", legend =levels(plot_list$classes_new)[16:25], 
       pch=16, pt.cex=1, cex=1, bty='n',
    col = cluster_col[16:25])
png()

png("main_figures/cre_heatmap.png",  width=10, height=12, unit="cm",
    res=1200)
f1 = circlize::colorRamp2(seq(-4, 4, length = 3), c("blue", "#EEEEEE", "red"))
rowAnnotation(cluster = as.factor(plot_list$classes_new[my_hclust_row]),
              show_legend=FALSE, 
              col=list(cluster=cluster_col)) +
Heatmap(dat_rna,  top_annotation = ha_left, col=f1,
        cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = F,
        show_column_names = F, name = "snRNAseq", show_heatmap_legend = T) +
Heatmap(dat_atac, top_annotation = ha_right, col=f1,
        cluster_columns = FALSE, cluster_rows = FALSE, show_row_names = F,
        show_column_names = F, name="snATACseq", show_heatmap_legend = T)
dev.off()

# Pathway Analysis

sce <- readRDS("scrna/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS")

all_genes <- lapply(ind_classes, function(x) rownames(plot_list_atac)[x])
all_genes <- all_genes[lengths(all_genes)>10]

res <- lapply(all_genes, function(x) go_term_analysis(x, rowData(sce)$gene_ids))

for (i in 1:length(res)){
    
    res[[i]] <- res[[i]]@result
    res[[i]] <- as.data.frame(res[[i]])
    res[[i]]$cluster <- names(res)[i]
    
}

res <- do.call(rbind, res)

write.table(res, file="snATACseq/processed_data/Go_terms_nmf_cluster.tsv", 
            sep="\t", quote=F, row.names=FALSE)
















