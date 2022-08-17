
################################################################################
#                                                                              #
#                                                                              #
# Quality control and UMAP plots for scRNA-seq data                            #
#                                                                              #
#                                                                              #    
################################################################################

library(ggplot2) 
library(SingleCellExperiment)
library(scuttle)
library(Polychrome)
library(readxl)
library(rcartocolor)


sce <- readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS")
col <- c(Astro= '#ffc857', CGE_dev='#c6d5c0', ID2= '#558140', 
         `L2/3_CUX2`= '#6e0614', L4_RORB= '#8b3843', `L5/6_THEMIS`='#a86a72',
         `L5/6_TLE4`='#c59ba1', LAMP5_CA1='#8eab80', MGE_dev='#ecd1c8',
         Micro='#484848',OPC='#92afc2', Oligo='#255f85', PN_dev='#e2cdd0',
         PV='#c77459', PV_SCUBE3='#daa290',`Poor-Quality`='#ffffff',
         SST='#b44622',VIP='#1c5701', Vas='#a3a3a3')


colData_sce <- colData(sce)
colData_sce$stage_ids <- factor(colData_sce$stage_ids, levels=
  rev(c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")))
colData_sce$major_clust <- factor(colData_sce$major_clust, levels=
  rev(c("L2/3_CUX2", "L4_RORB", "L5/6_THEMIS", "L5/6_TLE4", "PN_dev",
    "VIP", "ID2", "LAMP5_CA1", "CGE_dev", "SST", "PV", "PV_SCUBE3", "MGE_dev",
    "Astro", "OPC", "Oligo", "Micro", "Vas", "Poor-Quality")))
colData_sce <- as.data.frame(colData_sce)

g1 <- ggplot(data = colData_sce, mapping = aes(y = stage_ids, fill = major_clust)) + 
  geom_bar(position = "fill", col="black", size=0.2) + 
  scale_fill_manual(values = col) + 
  theme_classic() +
  guides(fill=FALSE)


df.summary <- colData_sce %>% 
  count(stage_ids, major_clust, batch)  %>%
  group_by(batch) %>% 
  mutate(per= prop.table(n))
  
df.summary1 <- df.summary %>%
  group_by(stage_ids, major_clust) %>% 
  summarise(
    sd = sd(per, na.rm = TRUE),
    mean = mean(per)
  )

ggplot(df.summary1 %>% filter(major_clust=="Oligo"), aes(mean, stage_ids)) +
  geom_col(fill = "#255f85", color = "black") +
  geom_errorbar(aes(xmin = mean, xmax = mean+sd), width = 0.2) +
  theme_minimal() + theme(axis.title.y = element_blank()) + xlab("Proportion") +
  ggtitle("Oligodendrocytes")

ggplot(df.summary1 %>% filter(major_clust=="Astro"), aes(mean, stage_ids)) +
  geom_col(fill = "#ffc857", color = "black") +
  geom_errorbar(aes(xmin = mean, xmax = mean+sd), width = 0.2) +
  theme_minimal() +  theme(axis.title.y = element_blank()) + xlab("Proportion") +
  ggtitle("Astrocytes")
