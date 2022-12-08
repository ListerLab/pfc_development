# Human prefrontal cortex development map repository

This repository contains the code and supporting information necessary to 
reproduce the analysis accompanying the manuscript "Human prefrontal cortex gene regulatory dynamics from gestation to adulthood at single-cell resolution". 

Note that names of modules and clusters in main and supplementary figures 
of the manuscript have occassionally be changed for readability. A table to 
convert the orginal names to the names used in the manuscript can be found
at `\annotation\Renaming_tables.xlsx`.

## File descriptions and paths

| File Paths | Description |
| ----------- | ----------- |
| `\snATACseq`  | Analysis pertaining to snATAC-seq |
| `\snATACseq\code`  | Code for analysis of snATAC-seq |
| `\snATACseq\browser`  | Code for interactive browser |
| `\snATACseq\R`  | Code for frequently used helper functions |
| `\snRNAseq`  | Analysis pertaining to snRNA-seq |
| `\snRNAseq\code`  | Code for analysis of snRNA-seq |
| `\snRNAseq\notebooks`  | Jupyter notebooks used for analysis of snRNA-seq |
| `\annotation`  | Contains meta data and supporting information |

<span style="font-size:smaller;"> Processed count matrices and metadata, including cell type annotations, can be found [here](https://console.cloud.google.com/storage/browser/neuro-dev/Processed_data). This Google bucket contains count matrices in the form of [anndata object](https://anndata.readthedocs.io/en/latest/) files (.h5ad) for our snRNA and scATAC data. The metadata files (.csv) for barcodes (BCs) and genes are the same as the obs and var data frames, respectively, within the anndata objects. The RNA-seq data has been processed as described in the methods section of the manuscript and detailed in jupyter notebooks prefixed 1-5 on the [publication GitHub](https://github.com/ListerLab/pfc_development/tree/master/snRNAseq/notebooks). Also included in the bucket are count and metadata for all cell types and the GABAergic inhibitory neurons separately. The default count matrix in the RNA anndata objects are the full-raw counts, but the downsampled CPM counts can be found in the anndata layer "ds_norm_cts." The ATAC anndata object contains processed counts as described in the manuscript. The embeddings presented in the publication are located in the obsm of the anndata objects. </span>

