
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(Matrix)

cd10_mhi.data <- Read10X(data.dir = "data/MhiMlo/Mhi/filtered_gene_bc_matrices")
cd10_mhi.data <- Read10X(data.dir = "data/MhiMlo/Mlo/filtered_gene_bc_matrices")

#Exclude Ig genes
IGHV.mhi.genes <- grep(pattern = "^IGHV", x = rownames(x = cd10_mhi.data), value = FALSE)
IGLV.mhi.genes <- grep(pattern = "^IGLV", x = rownames(x = cd10_mhi.data), value = FALSE)
IGKV.mhi.genes <- grep(pattern = "^IGKV", x = rownames(x = cd10_mhi.data), value = FALSE)
IGV_mhi_genes <- c(IGHV.mhi.genes, IGLV.mhi.genes, IGKV.mhi.genes)

IGHV.mlo.genes <- grep(pattern = "^IGHV", x = rownames(x = cd10_mlo.data), value = FALSE)
IGLV.mlo.genes <- grep(pattern = "^IGLV", x = rownames(x = cd10_mlo.data), value = FALSE)
IGKV.mlo.genes <- grep(pattern = "^IGKV", x = rownames(x = cd10_mlo.data), value = FALSE)
IGV_mlo_genes <- c(IGHV.mlo.genes, IGLV.mlo.genes, IGKV.mlo.genes)

#Exclude HB genes

HBA1.cd19.genes <- grep(pattern = "^HBA1", x = rownames(x = cd19.data), value = FALSE)
HBA2.cd19.genes <- grep(pattern = "^HBA2", x = rownames(x = cd19.data), value = FALSE)
HBB.cd19.genes <- grep(pattern = "^HBB", x = rownames(x = cd19.data), value = FALSE)
cd19Hb_genes <- c(HBA1.cd19.genes, HBA2.cd19.genes, HBB.cd19.genes)

HBA1.cd10_mhi.genes <- grep(pattern = "^HBA1", x = rownames(x = cd10_mhi.data), value = FALSE)
HBA2.cd10_mhi.genes <- grep(pattern = "^HBA2", x = rownames(x = cd10_mhi.data), value = FALSE)
HBB.cd10_mhi.genes <- grep(pattern = "^HBB", x = rownames(x = cd10_mhi.data), value = FALSE)
cd10_mhiHbgenes <- c(HBA1.cd10_mhi.genes, HBA2.cd10_mhi.genes, HBB.cd10_mhi.genes)

HBA1.cd10_mlo.genes <- grep(pattern = "^HBA1", x = rownames(x = cd10_mlo.data), value = FALSE)
HBA2.cd10_mlo.genes <- grep(pattern = "^HBA2", x = rownames(x = cd10_mlo.data), value = FALSE)
HBB.cd10_mlo.genes <- grep(pattern = "^HBB", x = rownames(x = cd10_mlo.data), value = FALSE)
cd10_mloHbgenes <- c(HBA1.cd10_mlo.genes, HBA2.cd10_mlo.genes, HBB.cd10_mlo.genes)

HBA1.rb.genes <- grep(pattern = "^HBA1", x = rownames(x = cd45rb.data), value = FALSE)
HBA2.rb.genes <- grep(pattern = "^HBA2", x = rownames(x = cd45rb.data), value = FALSE)
HBB.rb.genes <- grep(pattern = "^HBB", x = rownames(x = cd45rb.data), value = FALSE)
rbHbgenes <- c(HBA1.rb.genes, HBA2.rb.genes, HBB.rb.genes)

## Exclude IGV and HB genes from each dataset

cd10_mhi.data <- cd10_mhi.data[-IGV_mhi_genes, ]
cd10_mhi.data <- cd10_mhi.data[-cd10_mhiHbgenes, ]
cd10_mlo.data <- cd10_mlo.data[-IGV_mlo_genes, ]
cd10_mlo.data <- cd10_mlo.data[-cd10_mloHbgenes, ]

cd10_mhi <- CreateSeuratObject(counts = cd10_mhi.data, project = "All samples", min.cells = 3, min.features = 200)
cd10_mlo <- CreateSeuratObject(counts = cd10_mlo.data, project = "All samples", min.cells = 3, min.features = 200)

cd10_mhi[["percent.mt"]] <- PercentageFeatureSet(object = cd10_mhi, pattern = "^MT-")
cd10_mlo[["percent.mt"]] <- PercentageFeatureSet(object = cd10_mlo, pattern = "^MT-")

cd10_mhi <- subset(x = cd10_mhi, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 10)
cd10_mlo <- subset(x = cd10_mlo, subset = nFeature_RNA > 500 & nFeature_RNA < 2000 & percent.mt < 10)

# Identify B cell genes
Bcell.genes <- c("CD79A", "CD79B", "CD19", "MS4A1")

cd10_mlo[["percent.Bcell"]] <- PercentageFeatureSet(object = cd10_mlo, features = Bcell.genes)
cd10_mhi[["percent.Bcell"]] <- PercentageFeatureSet(object = cd10_mhi, features = Bcell.genes)

cd10_mlo <- subset(x = cd10_mlo, subset = percent.Bcell > 0.1)
cd10_mhi <- subset(x = cd10_mhi, subset = percent.Bcell > 0.1)

#Add column naming subset
cd10_mlo$Type <- "CD10 IgM low"
cd10_mhi$Type <- "CD10 IgM high"

#Merge datasets
CD10merge <- merge  (x = cd10_mhi, y = cd10_mlo)

#Normalize
NormCD10merge <- CD10merge
NormCD10merge <- NormalizeData(object = NormCD10merge)

NormBcellCD10merge <- NormalizeData(object = NormBcellCD10merge)

#Find Variable genes
NormCD10merge <- FindVariableFeatures(object = NormCD10merge)

#Scale data
NormCD10merge <- ScaleData(object = NormCD10merge)

# Run PCA
NormCD10merge <- RunPCA(object = NormCD10merge)

# PC Elbow plot
ElbowPlot(NormCD10merge)

# UMAP of TS subsets only
NormCD10merge <- FindNeighbors(object = NormCD10merge, dims = 1:30)
NormCD10merge <- FindClusters(object = NormCD10merge, resolution = 1.0)
NormCD10merge <- RunUMAP(object = NormCD10merge, reduction = "pca", dims = 1:30)
DimPlot(object = NormCD10merge, reduction = "umap", label = TRUE)
DimPlot(object = NormCD10merge, reduction = "umap", group.by = "Type", label = TRUE)
