library(Seurat)
library(ggplot2)
library(patchwork)
library(monocle3)
library(RColorBrewer)
library(slingshot)
library(gam)
library(plotly)
library(viridis)
library(zoo)

# ------ REDUCE -----
integrated_data_CD19 <- readRDS("data/processed/integrated.rds")

DefaultAssay(integrated_data_CD19) <- 'integrated'
reduced_data_CD19 <- RunPCA(integrated_data_CD19, assay = 'integrated', 
                            reduction.name = 'pca_rna', 
                            reduction.key = 'RNApca_',
                            verbose=T)
dims = 23
reduced_data_CD19 <- FindNeighbors(reduced_data_CD19, dims=1:dims, reduction = 'pca_rna')
reduced_data_CD19 <- FindClusters(reduced_data_CD19, resolution = 0.8)
  
reduced_data_CD19 <- RunUMAP(reduced_data_CD19, 
                               assay = "integrated", 
                               reduction = 'pca_rna', 
                               dims = 1:dims,
                               n.components = 3L,
                               reduction.name = 'umap3D_RNA',
                               reduction.key = 'UMAP_',
                               verbose=F)

saveRDS(reduced_data_CD19, "data/processed/dim_reduced.rds")
