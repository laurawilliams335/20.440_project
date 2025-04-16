Overview
This repository contains the data and code needed to produce Figures 4 and 5 from the first draft of our report. The matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz were 
downloaded from GEO Accession GSE289268, which contains the single cell sequencing data from Rodrigues et al. (Cell, 2025). Using this datset, we ran principal component analysis 
on the single cell sequencing data and clustered the cells to generate a UMAP plot. We then identified genes upregulated in each of the clusters and, focusing on cluster 6 (Roryt+ APCs), ran
and gene ontology enrichment analysis to produce a dot plot of upregulated pathways in this cluster. 

- Data
At a high level, how was the data generated?
If it’s too large to upload to your GitHub, where can it be
accessed? Include citations, if any.
- Folder structure
At a high level, what is in each folder and subfolder?
- Installation
How do I run your code?
What software and package versions do I need to install?

Reference for dataset source:
Rodrigues, M. E., Moreira, T. G., Canesso, M. C. C., et al. (2025). Rorγt-positive dendritic cells are required for the induction of peripheral regulatory T cells in response to oral 
antigens. Cell. Advance online publication. https://doi.org/10.1016/j.cell.2025.03.020
