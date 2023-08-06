# BioinformaticsResourcesProject
Project for the Bioinformatics Resources course held by Professor Alessandro Romanel in 2021 for the course of Quantitative and Computational Biology at University of Trento

The project starts with an `RData` file representing RNA-seq count data extracted from different [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) cancer datasets. This file is only used as a sample, so it contains 50 cases (tumor samples) and 50 controls (normal samples).

The data is provided in three dataframes:
 - `raw_counts_df` = contains the raw RNA-seq counts
 - `c_anno_df` = contains the sample name and condition (Case and Control)
 - `r_anno_df` = contains the ENSEMBL gene ids, length of the genes and gene symbols


The analysis I performed involves:
 1. Selecting only **protein coding** genes
 2. Performing a **Differential Gene Expression analysis** using the `edgeR` package
 3. Perform a **Gene set Enrichment Analysis** using the `clusterProfiler` package over Gene Ontology (GO) and KEGG, and visualize one pathway with `pathview`
 4. Identify which transcription Factors (TFs) have enriched scores in the promoters of all up-regulated genes
 5. Identify which up-regulated genes have a region in their promoter with significant binding scores (compared to a previously computed threshold)
 6. Access and use the **STRING** database to find Protein-Protein-Interactions (PPI) among differentially expressed genes and build a network with the `igraph` package


All plots are made with the `plotly` package, making them interactive and publication-ready.
