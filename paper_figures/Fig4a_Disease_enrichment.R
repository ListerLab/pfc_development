library(data.table)
library(disgenet2r)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(viridis)
library(scater)
library(cba)
library(dplyr)

all_res2 <- read.table("snRNAseq/processed_data/Diseases_enrichment.tsv")

all_res2$padj[all_res2$padj>0.05] <- NA
all_res2$padj[all_res2$overlap<5] <- NA

split_row <- split(all_res2$padj, all_res2$Cell_type_stages)
split_row <- names(split_row)[!sapply(split_row, 
                                      function(x) all(is.na(x)))]
all_res2 <- all_res2[all_res2$Cell_type_stages %in% split_row,]

split_row_order <- unique(all_res2[,c("Cell_type_stages", "Stage")])
split_row_order <- split_row_order$Cell_type_stages[order(split_row_order$Stage)]

all_res2$Cell_type_stages <- factor(all_res2$Cell_type_stages, split_row_order)

tmp <- as.data.frame(all_res2[, c("Cell_type_stages", "pathway", "padj")])
tmp$padj <- as.numeric(tmp$padj)
tmp$padj[is.na(tmp$padj)] <- 1
tmp <- tidyr::spread(tmp, "pathway", "padj")
rownames(tmp) <- tmp$Cell_type_stages
tmp <- tmp[,-1]

hc <- hclust(dist(t(tmp)))
order_col <- order.optimal(dist(t(tmp)), hc$merge)

all_res2$pathway <- factor(all_res2$pathway, 
                           levels=names(order_col$order)[order_col$order])

blue_colours <- c("#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(all_res2, aes(x=pathway, y=Cell_type_stages, fill=padj)) + 
    geom_tile() + theme_classic() +
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust=1),
          panel.background=element_rect(fill="white", colour="white")) +
    scale_fill_gradientn(colours =
                             rev(grDevices::colorRampPalette(colors=blue_colours)(20)), na.value = "white") 
ggsave(g1, file=paste0("paper_figures/Fig4a_diseases_enrichment_reduced.svg"), height=7,
       width=7)