################################################################################
#                                                                              #
#                                                                              #
# Network plots                                                                #
#                                                                              #
#                                                                              #    
################################################################################

library(igraph)
library(ggraph)
library(ggrepel)
library(ggiraph)

cell_types <- c("L4_RORB", "L2_3_CUX2", "L5_6_TLE4", "L5_6_THEMIS")

for(i in 1:length(cell_types)){

    infancy_network <- readRDS(paste0(
    "snATACseq/processed_data/infancy_network_ensheatment_",
    cell_types[i], ".RDS"))
    adol_network <- readRDS(paste0(
    "snATACseq/processed_data/adolescence_network_ensheatment_",
     cell_types[i], ".RDS"))


    g_i <- graph_from_data_frame(infancy_network, 
        directed = FALSE, vertices = NULL)
    comps <- components(g_i)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_i <- delete_vertices(g_i, remove_nodes)
    df_i <- data.frame(betweeness=betweenness(g_i), 
        hub=hub_score(g_i)$vector)
    df_i$gene <- rownames(df_i)

    inf_gg <- ggplot(df_i, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_i[df_i$betweeness>200&df_i$hub>0.5,], 
      aes(label=gene)) + ggtitle(paste0("Infancy: ", cell_types[i]))
    ggsave(inf_gg, file=paste0("supp_figures/network_", cell_types[i],
        "_infancy.svg"), height=7, width=7)

    g_a <- graph_from_data_frame(adol_network, directed = FALSE, 
        vertices = NULL)
    comps <- components(g_a)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_a <- delete_vertices(g_a, remove_nodes)
    df_a <- data.frame(betweeness=betweenness(g_a), hub=hub_score(g_a)$vector)
    df_a$gene <- rownames(df_a)

    a_gg <- ggplot(df_a, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_a[df_a$betweeness>200&df_a$hub>0.5,], 
          aes(label=gene)) + ggtitle(paste0("Adolescence: ", cell_types[i]))
    ggsave(a_gg, file=paste0("supp_figures/network_",
        cell_types[i], "_adolescence.svg"),
       height=7, width=7)

}


### Network plots makers

library(igraph)
library(ggraph)
library(ggrepel)
library(ggiraph)

cell_types <- c("L4_RORB", "L2_3_CUX2", "L5_6_TLE4", "L5_6_THEMIS")

for(i in 1:length(cell_types)){

    brain_network <- readRDS(paste0(
    "snATACseq/processed_data/neonatal_network_brain_",
    cell_types[i], "_degs.RDS"))
    organoids_network <- readRDS(paste0(
    "snATACseq/processed_data/neonatal_network_organoids_",
     cell_types[i], "_degs.RDS"))


    g_i <- graph_from_data_frame(brain_network, 
        directed = FALSE, vertices = NULL)
    comps <- components(g_i)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_i <- delete_vertices(g_i, remove_nodes)
    df_i <- data.frame(betweeness=betweenness(g_i), 
        hub=hub_score(g_i)$vector)
    df_i$gene <- rownames(df_i)

    brain_gg <- ggplot(df_i, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_i[df_i$betweeness>200&df_i$hub>0.5,], 
      aes(label=gene)) + ggtitle(paste0("Brain: ", cell_types[i]))
    ggsave(brain_gg, file=paste0("supp_figures/network_", cell_types[i],
        "_brain_degs.svg"), height=7, width=7)

    g_a <- graph_from_data_frame(organoids_network, directed = FALSE, 
        vertices = NULL)
    comps <- components(g_a)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_a <- delete_vertices(g_a, remove_nodes)
    df_a <- data.frame(betweeness=betweenness(g_a), hub=hub_score(g_a)$vector)
    df_a$gene <- rownames(df_a)

    organoids_gg <- ggplot(df_a, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_a[df_a$betweeness>200&df_a$hub>0.5,], 
          aes(label=gene)) + ggtitle(paste0("Organoids: ", cell_types[i]))
    ggsave(organoids_gg, file=paste0("supp_figures/network_",
        cell_types[i], "_organoids_degs.svg"),
       height=7, width=7)

    }

