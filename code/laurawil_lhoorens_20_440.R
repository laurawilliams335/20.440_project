## FIGURE 1 - Volcano Plot
library(Seurat)
library(ggplot2)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

# Load the Seurat object, 'fig4' references the figure in Canesso et al. that used this dataset
fig4 <- readRDS("Fig4.rds")

# Fetch specific data from Seurat object for UMAP plots
fig4.data <- FetchData(fig4, vars = c("UMAP_1", "UMAP_2", "RNA_snn_res.0.8", "Conditions", "Biotin", "Treatment", "Cd8a.log10", "Cd103.log10"))

# Set factor levels for the Biotin variable which designates cells as biotin negative or positive
fig4.data$Biotin <- factor(fig4.data$Biotin, levels = c("DCBiop", "DCBion"))

# Find differentially expressed genes between biotin positive and negative
deg_results <- FindMarkers(
  object = fig4,
  ident.1 = "DCBiop",
  ident.2 = "DCBion",
  group.by = "Biotin",
  logfc.threshold = 0,  # Set to 0 if you want all DEGs regardless of fold change
  min.pct = 0.1,        # Minimum percent of cells expressing the gene
  test.use = "wilcox"   # You can also use "MAST", "bimod", etc.
)

# Display volcano plot
EnhancedVolcano(deg_results,
                lab = rownames(deg_results),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Biotin+ DCs vs Biotin- DCs')

########################################################################################################################################################

## FIGURE 2A
install.packages("Matrix")
library(Matrix)

# Set file path
data_dir <- "~/Downloads/20.440"  

# Read the matrix
counts <- readMM(file.path(data_dir, "matrix.mtx.gz"))

# Read features and barcodes
features <- read.delim(gzfile(file.path(data_dir, "features.tsv.gz")), header = FALSE)
barcodes <- read.delim(gzfile(file.path(data_dir, "barcodes.tsv.gz")), header = FALSE)

# Assign row and column names
rownames(counts) <- make.unique(as.character(features$V2))
colnames(counts) <- as.character(barcodes$V1)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = counts, project = "GSE289268")

# Calculate % mitochondrial genes (genes starting with "mt-")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Filter cells 
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 5000 & 
                       percent.mt < 10)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj)

# 1. Run neighbors based on PCA
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# 2. Find clusters,  resolution controls number of clusters, used 0.18 to get 9 clusters matching the paper
seurat_obj <- FindClusters(seurat_obj, resolution = 0.18)  
# 3. Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# 4. Plot UMAP
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

## FIGURE 2B

print(VlnPlot(seurat_obj, features = c("Rorc", "Prdm16", "Cd7", "Esam", "Dtx1", "Il18r1", "Ccr6", "Aire"), stack=TRUE, flip=TRUE, group.by = "seurat_clusters"))
print(VlnPlot(seurat_obj, features = c("Rorc", "Aire"), stack=TRUE, flip=TRUE, group.by = "seurat_clusters"))

########################################################################################################################################################

## FIGURE 3
genes_of_interest <- c("Rorc", "Ccr6", "Prdm16", "Cd7", "Esam", "Dtx1", "Il18r1", "Xcr1", "Cd24a", "Naaa", "Clec9a", "Aire")
dot_plot <- DotPlot(fig4, 
                    features = genes_of_interest,  # Your list of genes
                    group.by = "Conditions")  # Column in metadata specifying groups
dot_plot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for clarity
  labs(title = "Gene Expression in Biop vs Bion", 
       x = "Conditions", 
       y = "Genes") +
  scale_color_gradient(low = "lightblue", high = "darkblue")  # Adjust colors as desired


########################################################################################################################################################

## FIGURE 4 - GSEA on Biotin+ vs Biotin-
# Set identity class to the 'conditions' column, which has infection status and biotin +/-
Idents(fig4) <- "Conditions"

# Compute average expression since we're using a single cell dataset
avg_exp <- AggregateExpression(fig4, return.seurat = FALSE)$RNA

# Subset to just two conditions
conditions_to_keep <- c("OV-DCBiop", "OV-DCBion")
avg_exp <- avg_exp[, conditions_to_keep]

# Add gene names as a column
gct_data <- as.data.frame(avg_exp)
gct_data <- tibble::rownames_to_column(gct_data, var = "Name")

# Add dummy description column (just gene name again)
gct_data$Description <- gct_data$Name

# Reorder columns: Name, Description, expression values
gct_data <- gct_data[, c("Name", "Description", setdiff(names(gct_data), c("Name", "Description")))]

# Write GCT header and data to file
gct_file <- "expression_data.gct"
writeLines(c("#1.2", 
             paste(nrow(gct_data), ncol(gct_data) - 2, sep = "\t")),
           con = gct_file)

# Append data
write.table(gct_data, file = gct_file, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Define classes and class assignments
class_names <- c("OV-DCBiop", "OV-DCBion")  # The two conditions
class_assignments <- c("OV-DCBiop", "OV-DCBion")  # Each pseudo-bulk is labeled as OV-DCBiop or OV-DCBion

# The total number of samples is 2 (one per condition)
num_samples <- length(class_assignments)
num_classes <- length(class_names)

# Define output .cls file
cls_file <- "OV-DCBiop_vs_OV-DCBion.cls"

# Write the .cls file with correct format
writeLines(c(
  paste("#", num_samples, num_classes, "1"),   # Line 1: "# total_samples total_classes 1"
  paste(class_names, collapse = " "),         # Line 2: Class names (separated by space)
  paste(class_assignments, collapse = "\t")   # Line 3: Class assignments (separated by tab)
), con = cls_file)

############################################################################
## Use Broad Institute GSEA program to generate top enriched pathways from .cls file
## OR download top20_enriched_pathways.csv from our GitHub!
############################################################################

## Dot plot
# Load your CSV 
top_pathways <- read.csv("top20_enriched_pathways.csv")

# Sort so most enriched are at the top
top_pathways <- top_pathways[order(top_pathways$NES), ]

# Plot
library(ggplot2)

ggplot(top_pathways, aes(x = NES, y = reorder(pathway_name, NES))) +
  geom_point(aes(color = pvalue), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "p-value") +
  labs(title = "Top Enriched Pathways in Biotin+ DCs",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway") +
  theme_minimal(base_size = 14)

########################################################################################################################################################

## FIGURE 5 - GSEA on cDC2 and cDC1 clusters

Idents(fig4) <- "RNA_snn_res.0.8"
conditions_to_keep <- c("g0", "g1") ## these groups are the cDC1 and cDC2 clusters
avg_exp <- AggregateExpression(fig4, return.seurat = FALSE)$RNA
avg_exp_sv <- avg_exp[, conditions_to_keep]

# Add gene names as a column
gct_data_sv <- as.data.frame(avg_exp_sv)
gct_data_sv <- tibble::rownames_to_column(gct_data_sv, var = "Name")

# Add dummy "Description" column (can just be gene name again)
gct_data_sv$Description <- gct_data_sv$Name

# Reorder columns: Name, Description, expression values
gct_data_sv <- gct_data_sv[, c("Name", "Description", setdiff(names(gct_data_sv), c("Name", "Description")))]

# Write GCT header and data to file
gct_file_sv <- "expression_data_sv_bycluster.gct"
writeLines(c("#1.2", 
             paste(nrow(gct_data_sv), ncol(gct_data_sv) - 2, sep = "\t")),
           con = gct_file_sv)

# Append data
write.table(gct_data_sv, file = gct_file_sv, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)

# Define classes and class assignments
class_names <- c("c0", "c1")  # The two conditions
class_assignments <- c("c0", "c1")  # Each pseudo-bulk is labeled as OV-DCBiop or OV-DCBion

# The total number of samples is 2 (one per condition)
num_samples <- length(class_assignments)
num_classes <- length(class_names)

# Define output .cls file
cls_file <- "c0_vs_c1.cls"

# Write the .cls file with correct format
writeLines(c(
  paste(num_samples, num_classes, "1"),   # Line 1: "# total_samples total_classes 1"
  paste("#", class_names, collapse = " "),         # Line 2: Class names (separated by space)
  paste(class_assignments, collapse = "\t")   # Line 3: Class assignments (separated by tab)
), con = cls_file)

############################################################################
## Use Broad Institute GSEA program to generate top enriched pathways from .cls file
## OR download top20_enriched_pathways_c1.csv from our GitHub!
############################################################################


## Dot plot
# Load your CSV
top_pathways_c1 <- read.csv("top20_enriched_pathways_c1.csv")

# Sort so most enriched are at the top
top_pathways_c1 <- top_pathways_c1[order(top_pathways_c1$NES), ]

# Plot
library(ggplot2)

ggplot(top_pathways_c1, aes(x = NES, y = reorder(pathway_name, NES))) +
  geom_point(aes(color = pvalue), size = 5) +
  scale_color_gradient(low = "red", high = "blue", name = "p-value") +
  labs(title = "Top Enriched Pathways in c1",
       x = "Normalized Enrichment Score (NES)",
       y = "Pathway") +
  theme_minimal(base_size = 14)


########################################################################################################################################################

## FIGURE 6 - GO enrichment analysis on cluster 6 Roryt+ APCs
cluster_markers <- FindAllMarkers(seurat_obj,
                                  only.pos = TRUE,       
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

library(clusterProfiler)
library(org.Mm.eg.db)
cluster_markers_6 <- FindMarkers(seurat_obj, ident.1 = "6", only.pos = TRUE)
gene_list <- rownames(cluster_markers[cluster_markers_6$p_val_adj < 0.05, ])
entrez_ids <- bitr(gene_list, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Run GO enrichment
go_results <- enrichGO(gene         = entrez_ids$ENTREZID,
                       OrgDb        = org.Mm.eg.db,
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       readable      = TRUE)

dotplot(go_results)


########################################################################################################################################################

## FIGURE 7 - AIRE Expression in both datasets
dot_plot <- DotPlot(fig4, 
                    features = c("Aire", "Siglec15"),  # Your list of genes
                    group.by = "Conditions")  # Column in metadata specifying groups
dot_plot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for clarity
  labs(title = "AIRE Expression", 
       x = "Conditions", 
       y = "Genes") +
  scale_color_gradient(low = "lightblue", high = "darkblue")  # Adjust colors as desired

dot_plot <- DotPlot(seurat_obj, 
                    features = c("Aire", "Siglec15"),  # Your list of genes
                    group.by = "seurat_clusters")  # Column in metadata specifying groups
dot_plot + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for clarity
  labs(title = "AIRE Expression", 
       x = "Conditions", 
       y = "Genes") +
  scale_color_gradient(low = "lightblue", high = "darkblue")  # Adjust colors as desired
