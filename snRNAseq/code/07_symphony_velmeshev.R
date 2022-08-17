#Symphony

library(BiocManager)
library(scater)
library(scran)
library(scuttle)
library(reticulate)
library(Seurat)
library(uwot)
library(symphony)

reference = readRDS("snRNAseq/processed_data/2020-12-18_whole-tissue_post-restaged-GABA-clustering.RDS") 
rawref =readRDS('snRNAseq/processed_data/2020-12-18_RAW_whole-tissue_post-restaged-GABA-clustering.RDS')
rownames(rawref) = rowData(rawref)$gene_ids
bycol = 1000/colSums(counts(rawref))            #ensures downsample target is 1000
set.seed(100)
ref_exp = downsampleMatrix(x = counts(rawref), prop = bycol, bycol = TRUE)
ref_metadata = reference@colData[c(71,69,68,66,62,59)]
initembedding = reducedDim(reference, "UMAP")

hvg = subset(rowData(reference)$gene_ids, 
    rowData(reference)$highly_variable == TRUE)
set.seed(0)

new_sym_ref = symphony::buildReference(
  ref_exp,
  ref_metadata,
  K = 100,
  verbose = TRUE,
  do_umap = TRUE,
  do_normalize = TRUE,
  vargenes_method = 'vst',
  topn = 0,
  d = 365,
  additional_genes = hvg,
  save_uwot_path = '.symphony_trial2'
)

krieg = readRDS('snRNAseq/processed_data/transfer_learning/datasets/velmeshev_dataset.RDS') 
#this data is available https://www.science.org/doi/full/10.1126/science.aav8130

ribogeneset = read.csv('annotation/ribo_geneset.txt', header = FALSE)$V1 

ribout = subset(rowData(krieg)$ID, rowData(krieg)$Symbol %in% ribogeneset)

krieg = krieg[!( rownames(krieg) %in% ribout),]
krieg = krieg[rowData(krieg)$Symbol != 'MALAT1',]
krieg =krieg[-which(str_detect(rowData(krieg)$Symbol, '^MT-')),]
krieg = krieg[rowSums(counts(krieg)>0)>=5,]
bycol = 1000/colSums(counts(krieg))            #ensures downsample target is 1000
set.seed(100)
krieg_exp = counts(krieg)
krieg_metadata = krieg@colData[c(1,2,3,8,13,14,15)]

query = mapQuery(krieg_exp, krieg_metadata, new_sym_ref, 
    do_normalize = TRUE, do_umap = TRUE)

reference_umap = umap(t(new_sym_ref$Z_corr), 
        n_neighbors = 25, ret_model = TRUE, init = initembedding, 
        a = 0.58303, b = 1.334167, spread = 1, min_dist = 0.5, 
        verbose = TRUE, n_epochs = 0)
qumap = umap_transform(t(query$Z), reference_umap, n_epochs = 50)

qumap = as.data.frame(qumap)
qumap$dataset = character(nrow(qumap))
qumap$dataset[1:nrow(qumap)] = 'velmeshev'

qumap$pred = knn.reg(train = refmap[,1:2], 
                     test = qumap[,1:2], y = reference$arcsin_ages, k = 10)$pred
qumap$cluster = knn(train = refmap[,1:2], 
                    test = qumap[,1:2], cl = reference$cluster, k = 10)
qumap$age = query$age
qumap$diff = abs(qumap$pred - qumap$age)

qumap$sample = query$sample
qumap$sample = fct_reorder(qumap$sample, qumap$age)
qumap = subset(qumap, !(sample %in% 
                          c('5893_PFC', '5538_PFC_Nova', '5879_PFC_Nova','5936_PFC_Nova', 
                            '5408_PFC_Nova', '5144_PFC', '5278_PFC', '5403_PFC', 
                            '5419_PFC', '5945_PFC')))
smolquery = query[, !(query$sample %in% c('5893_PFC', '5538_PFC_Nova', 
                                          '5879_PFC_Nova','5936_PFC_Nova', '5408_PFC_Nova', '5144_PFC', 
                                          '5278_PFC', '5403_PFC', '5419_PFC', '5945_PFC'))]


myknn = function(ref_embed, ref_class, test_embed, test_class, k){
  
  output = list()
  kest = knn(train = ref_embed, test = test_embed, 
             cl = as.factor(ref_class), prob=TRUE, k =k)
  result = character(nrow(test_embed))
  
  for(i in 1:nrow(test_embed)){
    if(kest[i] == as.factor(test_class)[i]) {     
      result[i] = "correct"
    }
    else {
      result[i] = "incorrect"
    }
  }
  
  output$result = result
  correct = subset(test_class, result == "correct")
  incorrect = subset(test_class, result == "incorrect")
  output$correct = correct
  output$incorrect = incorrect
  output$accuracy = length(correct)/length(result)
  output$corprop = summary(correct)/length(correct)
  output$est = kest
  output$percents = summary(correct)/summary(test_class)
  
  return(output)
  
}


out = myknn(ref_embed = refmap[,1:2], test_embed = qumap[,1:2], 
            test_class = smolquery$cluster, ref_class = reference$cluster, k =10)
qumap$result = out$result
qumap$truecluster = smolquery$cluster

astro = subset(qumap, cluster == 'Astro' & result == 'correct')
l23 = subset(qumap, cluster == 'L2/3_CUX2' & result == 'correct')

plotastro = ggplot(astro, aes(x = pred, y = sample, fill = ..x..))+ geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted astrocyte ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2,3, 4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
                     label = c('ga22', 
                               '34d', 
                               '301d', 
                               '3yr', '10yr',
                               '20yr', '40yr'))+
  scale_fill_viridis_c()



plotl23 = ggplot(l23, aes(x = pred, y = sample, fill = ..x..)) + geom_density_ridges_gradient(bandwidth = 0.25) + 
  xlim(-1,5.5) + theme_classic() +
  ggtitle('Predicted L2/3 ages by sample') + xlab('predicted age') +
  scale_x_continuous(breaks = c(as.numeric(levels(as.factor(reference$arcsin_ages)))[-c(2,3, 4,6,7,8,10,11,12,14,15,16,18, 19, 20,21,23)]),
                     label = c('ga22', 
                               '34d', 
                               '301d', 
                               '3yr', '10yr',
                               '20yr', '40yr'))+
  scale_fill_viridis_c()


ggsave(plotastro, 'supp_figures/velmeshev_symphony_astro_age_prediction.png')
ggsave(plotl23, 'supp_figures/velmeshev_symphony_l23_age_prediction.png')

saveRDS(qumap, "snRNAseq/processed_data/transfer_learning/symphony_velmeshev_umap.RDS")

#process can be repeated to include additional datasets


