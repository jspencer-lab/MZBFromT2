library(Seurat)
library(ggplot2)
library(dplyr)
library(plotly)
library(monocle3)
library(RColorBrewer)
library(reticulate)

data <- Read10X(data.dir = "data/MZBTrajectory/Sample1/filtered_feature_bc_matrix")
cd19 <-  CreateSeuratObject(counts = data$`Gene Expression`)
cd19.adt <- data$'Antibody Capture'
cd19[["ADT"]] <- CreateAssayObject(counts = cd19.adt)

## Subset cd19+

feature_data_AHH1 <- FetchData(object = cd19, vars = c('HTO-AHH1-TotalSeqC'), slot='counts')
feature_data_AHH2 <- FetchData(object = cd19, vars = c('HTO-AHH2-TotalSeqC'), slot='counts')

CD19_cell_names =  colnames(cd19)[which(x = feature_data_AHH1 > 30 & feature_data_AHH2 < 4)]
cd19 <- SubsetData(cd19, cells = CD19_cell_names)

#Remove IGHV genes

ig_hl_chain_gene_patterns <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^AC233755")
hashes <-  "^HTO-"
features_to_keep <- c(rownames(cd19@assays$ADT), rownames(cd19@assays$RNA))
for (f in c(ig_hl_chain_gene_patterns, hashes)){
  features_to_keep <- features_to_keep[!grepl(paste0(f, collapse = "|"), features_to_keep)]
}
cd19 <- subset(cd19, features = features_to_keep)

# Filter using B cell genes

Bcell.genes <- c("CD79A", "CD79B","MS4A1")

cd19[["percent.Bcell"]] <- PercentageFeatureSet(object = cd19, features = Bcell.genes)

cd19 <- subset(x = cd19, subset = percent.Bcell > 0.1)


## QC

cd19[["percent.mt"]] <- PercentageFeatureSet(cd19, pattern = "^MT-")


cd19 <- subset(cd19, subset = nFeature_RNA > 200 & nFeature_RNA < 1800 & percent.mt < 13, )

cd19 <- NormalizeData(cd19)

cd19 <- FindVariableFeatures(cd19, selection.method = "vst", nfeatures = 2000)

cd19 <- ScaleData(cd19)

cd19 <- RunPCA(cd19, verbose = FALSE)

cd19 <- FindNeighbors(cd19, dims = 1:15)
cd19 <- FindClusters(cd19, resolution = 5.2)
cd19 <- RunUMAP(cd19, dims = 1:15)

## Identify clusters
memory <- subset(cd19, ident = c(5,11,8))

naive <- subset(cd19, ident = c(1,2,6,7,9,12,13,14,15,16,18,19,20,21,24,25,27,28,29,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45))
memory <- subset(cd19, ident = c(0,3,4,5,10,11,17,22,23,26,30,46,47))

clusters = levels(memory@active.ident)
for (c in clusters){
  cells <- WhichCells(memory, idents = c)
  print (median(as.numeric(memory[['RNA']]@counts['CD1C',cells])))
}

## Rename clusters, IGM high naive being clusters with top 30% Igm expression
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(8)), value = "TS")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(10,11,3)), value = "MZB")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(4,5,22,23,30,46,26)), value = "CSM")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(1,6,9,12,14,15,16,18,20,21,24,25,28,29,31,32,33,34,35,36,38,40,42,44)), value = "Naive")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(2,7,13,19,27,37,39,41,43,45)), value = "Naive Mhi")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(0,17,47)), value = "M ONLY")

## Run 3D UMAP
cd19 <- RunUMAP(cd19, dims = 1:15, reduction.key = 'UMAPRNA_', reduction.name = 'UMAP_rna', n.components = 3L)
