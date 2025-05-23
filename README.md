**Overview**

This repository contains the data and code needed to produce all figures from our report 'Single cell transcriptomic analysis of antigen-presenting cells responsible for peripheral tolerance in the gut'.
The file fig4.rds was downloaded from GEO Accession GSE281286, and contains the single cell sequencing data from Figure 4 of Canesso et al. (Science, 2024).
The matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz were downloaded from GEO Accession GSE289268, which contains the single cell sequencing data from Rodrigues et al. (Cell, 2025). 


**Data**

Canesso et al.: Dendritic cells were isolated from the duodenum gLN of 3-4 CD40G5/G5 male mice (8-12 weeks old). The DCs were index sorted as Aqua-CD45.2+CD45.1-Lin-(TCRb-B220-CD64-)CD11chiMHC-IIint/hi Biotin- or Biotin+. Following their modified Smart-seq2 protocols, RNA was extracted from the single cells to generate tagmented and barcoded cDNA using Nextera XT kits, which was sequenced using Illumina NextSeq 500. The DCs from the helminth-infected and uninfected mice were isolated from the duodenal lamina propria, co-stained with hashtag oligonucleotide labeled CD45 and MHC-I antibodies, and sorted as described above. 

Rodrigues et al.: Mesenteric lymph nodes were harvested from 15-day-old RorcE+7kbWT/WT mice and sorted for CD45+CD3−B220−CD11c+MHC-II+ DCs and CD45+CD3−B220−CD11c−CCR6+ ILC3s. Samples were pooled and sequenced using a 10x Genomics Chromium Single Cell Controller.

**Organization**

The code needed to produce all six figures in the report is in the file 'laurawil_lhoorens_20_440.R' in the folder 'code'
The necessary data is all in the folder titled 'data'. Figures 1, 3, 4, and 5 use Fig4.rds, while figures 2 and 6 use the data in matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz. Figure 7 requires both datasets.

**Running the code**

Download R and Rstudio to run the code. 
All of the necessary packages will be installed by running the entire script.
For figures 1-3, 6, and 7, the figures can be generated by running the script. However, for figures 4 and 5, an external application (https://www.gsea-msigdb.org/gsea/index.jsp) was used to run GSEA analysis. To generate the figures, run the script to generate the .cls file, and use this file to run GSEA. GSEA was run using the gene set database m5.all.v2024.1.Mm.symbols.gmt, and the Mouse_Seq_Accession_MSigDB.v2024.1.Mm.chip Chip platform. Gene sets smaller than 15 or larger than 500 were excluded from the analysis. Permutation was done by gene_set because of the small sample number after data aggregation. Then, input the generated .csv file of top 20 upregulated pathways into the remaining code. Alternatively, the .csv files generated by our GSEA analysis are provided in the data folder. top20_enriched_pathways.csv is used for Figure 4, and top20_enriched_pathways_c1.csv was used for Figure 5. 

**Reference for dataset sources:**

Canesso, M. C. C., Rodríguez, M. E., Moreira, T. G., et al. (2025). Identification of antigen-presenting cell–T cell interactions driving immune responses to food. Science, 387, eado5088. https://doi.org/10.1126/science.ado5088

Rodrigues, M. E., Moreira, T. G., Canesso, M. C. C., et al. (2025). Rorγt-positive dendritic cells are required for the induction of peripheral regulatory T cells in response to oral 
antigens. Cell. Advance online publication. https://doi.org/10.1016/j.cell.2025.03.020
