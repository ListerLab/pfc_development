# neuro-dev-atac


## File descriptions and paths

Sample meta-data
```
annotation/scATACseq_neuronal_maturation.xlsx
```

Cell-type Tn5 insertion bed files for each library. One record per insertion event, flanked by 25 base window per side (to facilitate peak calling). Blacklist mapped reads removed. 
```
bedfiles/*_blacklistrm_strand.bed
```

Stage-grouped cell-type Tn5 insertion bed files. One record per insertion event, flanked by 25 base window per side (to facilitate peak calling). 
```
bedfiles/stages/
```

ATAC accessability peaks per cell type
```{bash}
ls peaks/*merge.narrowPeak | grep -v RL
```

```
peaks/Astro_merge.narrowPeak
peaks/CGE_der_merge.narrowPeak
peaks/L2_3_merge.narrowPeak
peaks/L4_merge.narrowPeak
peaks/L5_6_merge.narrowPeak
peaks/MGE_der_merge.narrowPeak
peaks/Micro_merge.narrowPeak
peaks/Oligo_merge.narrowPeak
peaks/OPC_merge.narrowPeak
peaks/Vas_merge.narrowPeak
```

List containing ATAC peak Tn5 insertion counts. 
```
processed_data/atac_peak_counts.Rds
```

The top level items in the list are lists for each cell type. Within these cell type lists are the raw peak instertion counts, and then normalised by insertions per million per kilobase.

For example:

To access the raw counts and normalised counts for Astrocyte peaks
```{r}
peak_mat_list <- readRDS("processed_data/atac_peak_counts.Rds")
astro_counts <- peak_mat_list$Astro$peak_counts
astro_fpkm <- peak_mat_list$Astro$peak_fpkm
```


### ATAC peaks

Per cell type, filtered ATAC peaks
```
processed_data/cell_type_atac_peaks_filtered_gr.Rds
```

