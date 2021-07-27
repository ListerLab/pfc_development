# Heatmap TF Ion Transporter

library(ggplot2)
library(cba)
library(reshape2)

tf_expr <- readRDS("snATACseq/processed_data/tf_cell_type_rna.RDS")
tf_expr <- lapply(tf_expr, function(x) names(x)[x==0])
all_gos <- readRDS("snATACseq/processed_data/Gos_tf_all.RDS")

gos <- c(Ensheatment="GO:0007272", `Cell Cycle`="GO:0007049", 
    `Neuron Migration`="GO:0001764", `Synapse Organization`="GO:0050808",
    `Cytoskeletal Organization`="GO:0051493", Learning="GO:0007612",
    `Cell Death`="GO:0008219", `Ion Transport`="GO:0006811")

enriched <- list()

for(i in 1:length(gos)){
    
    not_exprs <- tf_expr[names(all_gos)[i]]
    tmp_gos <- sapply(names(all_gos[[i]]), function(x) 
        strsplit(x, ".", fixed=T)[[1]][1])
    tmp_index <- split(1:length(tmp_gos), tmp_gos)
    tmp <- lapply(tmp_index, function(x) do.call(rbind, all_gos[[i]][x]))
    for(j in 1:length(tmp)){
        tmp[[j]]$GO <- sapply(rownames(tmp[[j]]), function(x) 
            strsplit(x, ".", fixed=T)[[1]][1])
        tmp[[j]]$trend <- sapply(rownames(tmp[[j]]), function(x) 
            strsplit(x, ".", fixed=T)[[1]][4])
        tmp[[j]]$cell_type <- sapply(rownames(tmp[[j]]), function(x) 
            strsplit(x, ".", fixed=T)[[1]][2])
        tmp[[j]]$pvalue <- 10^(-tmp[[j]]$mlog10p)
        tmp[[j]]$padjust <- p.adjust(tmp[[j]]$pvalue, method="fdr")
        tmp[[j]] <- tmp[[j]][tmp[[j]]$padjust<0.05,]
        tmp[[j]]$feature_short <-  sapply( tmp[[j]]$feature, function(x) 
            strsplit(x, "_", fixed=T)[[1]][1])
        if(dim(tmp[[j]])==0) next
        tmp[[j]] <- tmp[[j]][!(tmp[[j]]$feature_short %in% not_exprs),]
    } 
    
    tmp <- tmp[!sapply(tmp, function(x) dim(x)[1]==0)]

    enriched[[i]] <- tmp    
    
}

names(enriched) <- names(all_gos)
enriched <- unlist(enriched, recursive = FALSE)

enriched_transporter <- do.call(rbind, enriched[grepl("GO:0006811", 
    names(enriched), fixed=T)])
enriched_transporter$name <- paste0(enriched_transporter$cell_type,
    "_", enriched_transporter$trend)

a <- acast(enriched_transporter,  feature_short ~ name , mean, fill=1,
           drop=T, value.var="padjust")
hc <- hclust(dist(a))
order_row <- order.optimal(dist(a), hc$merge)
levels_motifs <- names(order_row$order)[order_row$order]
enriched_transporter$feature_short <- factor(enriched_transporter$feature_short,
    levels=levels_motifs)
enriched_transporter$name <- factor(enriched_transporter$name,
    levels=rev(c("L2_3_up", "L2_3_interup", "L2_3_down", "L4_up", "L4_interup", 
             "L4_down", "L5_6_up", "L5_6_interup", "L5_6_down",
             "CGE_der_up", "MGE_der_up", "MGE_der_interup", "Astro_up",
             "Oligo_up")))

write.table(enriched_transporter, file="Enriched_GOs.tsv", sep="\t", quote=F)

blue_colours <- c("#F7FBFF", "#C6DBEF", "#6BAED6",  "#2171B5",
                  "#08306B")

g1 <- ggplot(enriched_transporter, aes(x=feature_short, y=name, fill=padjust)) +
    geom_tile() + theme_classic() + scale_fill_gradientn(colours =
      rev(grDevices::colorRampPalette(colors=blue_colours)(20))) +
    theme(axis.text.x=element_text(size=8, angle=90, hjust=1), 
          panel.background=element_rect(fill="white", colour="white"),
          axis.title=element_blank()) 

ggsave(g1, file="paper_figures/Fig3k_TF_ion_transporter.svg", 
       height=6, width=10)
