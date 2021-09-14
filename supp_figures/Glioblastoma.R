library(ggpubr)
library(scater)
library(Seurat)
library(viridis)
library(ggplot.multistats)
require(cowplot)


sce <- readRDS("snRNAseq/processed_data/GBM_primary_dataset.RDS")
sce <- logNormCounts(sce)

stem_cell_core <- read.table("annotation/Stem_Cell_Core.txt",
                             sep="\t", skip=1, fill=T)
stem_cell_core <- stem_cell_core[1:140,]

g2m.genes <- cc.genes$g2m.genes
ind.g2m.genes <- match(g2m.genes, rowData(sce)$Symbol)
ind.g2m.genes <- ind.g2m.genes[!is.na(ind.g2m.genes)]

ind <- match(stem_cell_core[,3], rowData(sce)$Symbol)
ind <- ind[!is.na(ind)]
eigen <- prcomp(as.matrix(logcounts(sce[ind,])))

df <- reducedDim(sce, "UMAP")
colnames(df) <- c("UMAP1", "UMAP2")
df <- as.data.frame(df)
df$eigen <- eigen$rotation[,1]
df$mean <- colMeans(as.matrix(logcounts(sce[ind,])))
df$g2m_genes <- colMeans(as.matrix(logcounts(sce[ind.g2m.genes,])))
set.seed(2)
df <- df[sample(1:dim(df)[1]),]


logo_file <- "snRNAseq/processed_data/background.png"
img <- png::readPNG(logo_file)

p1 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = g2m_genes, fill = stat(mean)),
     funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title="G2M Phase"))

p2 <- ggplot(df, aes(x=UMAP1, y=UMAP2)) +
    background_image(img) +
    stat_summaries_hex(aes(z = eigen, fill = stat(mean)),
    funs = c('mean'), bins = 100) + scale_fill_viridis_c() + 
    theme_classic() + scale_x_continuous(limits=c(-6.063561,21.449993)) +
    scale_y_continuous(limits=c(-9.861524,15.003452)) +
    guides(fill=guide_legend(title="Stem Cell Core Network"))
        

ggsave(p1, file="supp_figures/G2M.png",
       height=7, width=7)
ggsave(p2, file="supp_figures/Stem_cell_score.png",
       height=7, width=7)




