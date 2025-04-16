install.packages("Matrix")
install.packages("Seurat")  

library(Matrix)
library(Seurat)

# Set your file path
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

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells (tweak thresholds based on plots)
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
