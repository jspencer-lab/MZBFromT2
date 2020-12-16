library(Seurat)
library(ggplot2)
library(dplyr)
library(plotly)
library(monocle3)

data <- Read10X(data.dir = "data/MZBTrajectory/Sample2/filtered_feature_bc_matrix")
cd19 <-  CreateSeuratObject(counts = data$`Gene Expression`)
cd19.adt <- data$'Antibody Capture'
cd19[["ADT"]] <- CreateAssayObject(counts = cd19.adt)

## Subset cd19+

feature_data_AHH1 <- FetchData(object = cd19, vars = c('HTO-AHH1-TotalSeqC'), slot='counts')
feature_data_AHH2 <- FetchData(object = cd19, vars = c('HTO-AHH2-TotalSeqC'), slot='counts')

CD19_cell_names =  colnames(cd19)[which(x = feature_data_AHH1 > 120 & feature_data_AHH2 < 25)]
cd19 <- SubsetData(cd19, cells = CD19_cell_names)

#Remove IGHV genes

ig_hl_chain_gene_patterns <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^IGKC", "^JCHAIN", "^AC233755")
hashes <-  "^HTO-"
features_to_keep <- c(rownames(cd19@assays$ADT), rownames(cd19@assays$RNA))
for (f in c(ig_hl_chain_gene_patterns, hashes)){
  features_to_keep <- features_to_keep[!grepl(paste0(f, collapse = "|"), features_to_keep)]
}
cd19 <- subset(cd19, features = features_to_keep)

## Filter out Non B cells

Bcell.genes <- c("CD79A", "CD79B", "CD19", "MS4A1")

cd19[["percent.Bcell"]] <- PercentageFeatureSet(object = cd19, features = Bcell.genes)

cd19 <- subset(x = cd19, subset = percent.Bcell > 0.1)


## QC 

cd19[["percent.mt"]] <- PercentageFeatureSet(cd19, pattern = "^MT-")
VlnPlot(cd19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


cd19 <- subset(cd19, subset = nFeature_RNA > 300 & nFeature_RNA < 1500 & percent.mt < 9  & nCount_RNA > 0 & nCount_RNA <7500)

cd19 <- NormalizeData(cd19)

cd19 <- FindVariableFeatures(cd19, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(cd19)
plot1 

cd19 <- ScaleData(cd19)

cd19 <- RunPCA(cd19, verbose = FALSE)

cd19 <- FindNeighbors(cd19, dims = 1:20)
cd19 <- FindClusters(cd19, resolution = 3.0)
cd19 <- RunUMAP(cd19, dims = 1:20)


## IdentifY clusters

memory <- subset(cd19, idents = c(0,4,14,20,21,23,24,27,28))

TSandnaive <- subset(cd19, idents = c(1,2,3,5,6,8,9,10,11,12,13,15,16,17,18,19,21,22,25,26,29))
naive <- subset(cd19, idents = c(1,2,3,5,6,8,9,10,11,12,15,16,18,19,20,21,22,25,26,29))


clusters = levels(naive@active.ident)
naive.IgM <- for (c in clusters){
  cells <- WhichCells(naive, idents = c)
  print (median(as.numeric(naive[['ADT']]@counts['ADT-IGM-TotalSeqC',cells])))
}

## Rename clusters, IGM high naive being clusters with top 30% Igm expression

cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(27,28)), value = "DN")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(30)), value = "PB")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(13,17)), value = "TS")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(14)), value = "MZB")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(7,23,24)), value = "CSM")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(2,3,5,6,8,9,11,15,16,18,19,20,22,29)), value = "Naive")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(1,10,12,21,25,26)), value = "Naive Mhi")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(0,4)), value = "M ONLY")

## Run 3D UMAP
cd19 <- RunUMAP(cd19, dims = 1:20, reduction.key = 'UMAPRNA_', reduction.name = 'UMAP_rna', n.components = 3L)
