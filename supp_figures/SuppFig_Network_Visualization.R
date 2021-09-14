library(igraph)
library(ggraph)
library(ggrepel)
library(ggiraph)

cell_types <- c("L4", "L2_3", "L5_6_TLE4", "L5_6_THEMIS")

for(i in 1:length(cell_types)){

    infancy_network <- readRDS(paste0(
    "snATACseq/processed_data/Infancy_Network_Ensheatment_",
    cell_types[i], ".RDS"))
    adol_network <- readRDS(paste0(
    "snATACseq/processed_data/Adolescence_Network_Ensheatment_",
     cell_types[i], ".RDS"))


    g_i <- graph_from_data_frame(infancy_network, 
        directed = FALSE, vertices = NULL)
    comps <- components(g_i)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_i <- delete_vertices(g_i, remove_nodes)
    g1 <- ggraph(g_i, layout = "centrality", cent=degree(g_i)) + 
        geom_edge_link(aes(color = strength<0)) + 
        geom_node_point(colour="gray") +
        geom_node_text(aes(label = name, filter=degree(g_i)>11), 
                       size = 5, repel = TRUE) +
        theme_graph() 
    ggsave(g1, file=paste0("supp_figures/SuppFig_Network_", cell_types[i],
                               "_Infancy_Net.svg"), height=7, width=7)
    g_i <- delete_vertices(g_i, remove_nodes)
    df_i <- data.frame(betweeness=betweenness(g_i), 
        hub=hub_score(g_i)$vector)
    df_i$gene <- rownames(df_i)

    inf_gg <- ggplot(df_i, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_i[df_i$betweeness>200&df_i$hub>0.5,], 
                    aes(label=gene)) 
    ggsave(inf_gg, file=paste0("supp_figures/SuppFig_Network_", cell_types[i],
        "_Infancy.svg"), height=7, width=7)

    g_a <- graph_from_data_frame(adol_network, directed = FALSE, 
        vertices = NULL)
    comps <- components(g_a)
    keep_clust <- names(which.max(table(comps$membership)))
    remove_nodes <- names(comps$membership)[comps$membership!=keep_clust]
    g_a <- delete_vertices(g_a, remove_nodes)
    g1 <- ggraph(g_a, layout = "centrality", cent=degree(g_a)) + 
        geom_edge_link(aes(color = strength<0)) + 
        geom_node_point(colour="gray") +
        geom_node_text(aes(label = name, filter=degree(g_a)>10), 
                       size = 5, repel = TRUE) +
        theme_graph() 
    ggsave(g1, file=paste0("supp_figures/SuppFig_Network_", cell_types[i],
                           "_Adolescence_Net.svg"), height=7, width=7)
    df_a <- data.frame(betweeness=betweenness(g_a), hub=hub_score(g_a)$vector)
    df_a$gene <- rownames(df_a)

    a_gg <- ggplot(df_a, aes(x=betweeness, y=hub)) + geom_point() + 
        theme_minimal() +
    geom_text_repel(data=df_a[df_a$betweeness>200&df_a$hub>0.5,], 
                    aes(label=gene)) 
    ggsave(a_gg, file=paste0("supp_figures/SuppFig_Network_",
        cell_types[i], "_Adolescence.svg"),
       height=7, width=7)
}
