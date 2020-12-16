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

CD19_cell_names =  colnames(cd19)[which(x = feature_data_AHH1 > 100 & feature_data_AHH2 < 20)]
cd19 <- SubsetData(cd19, cells = CD19_cell_names)

#Remove IGHV genes

ig_hl_chain_gene_patterns <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^AC233755")
hashes <-  "^HTO-"
features_to_keep <- c(rownames(cd19@assays$ADT), rownames(cd19@assays$RNA))
for (f in c(ig_hl_chain_gene_patterns, hashes)){
  features_to_keep <- features_to_keep[!grepl(paste0(f, collapse = "|"), features_to_keep)]
}
cd19 <- subset(cd19, features = features_to_keep)

## Remove non B cells

Bcell.genes <- c("CD79A", "CD79B", "MS4A1")

cd19[["percent.Bcell"]] <- PercentageFeatureSet(object = cd19, features = Bcell.genes)

cd19 <- subset(x = cd19, subset = percent.Bcell > 0.1)

## QC 

cd19[["percent.mt"]] <- PercentageFeatureSet(cd19, pattern = "^MT-")
VlnPlot(cd19, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


cd19 <- subset(cd19, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 15, )

cd19 <- NormalizeData(cd19)

cd19 <- FindVariableFeatures(cd19, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(cd19)
plot1 

cd19 <- ScaleData(cd19)

cd19 <- RunPCA(cd19, verbose = FALSE)

cd19 <- FindNeighbors(cd19, dims = 1:20)
cd19 <- FindClusters(cd19, resolution = 3.0)
cd19 <- RunUMAP(cd19, dims = 1:20)

## Identify  clusters

clusters = levels(cd19@active.ident)
 for (c in clusters){
  cells <- WhichCells(cd19, idents = c)
  print (median(as.numeric(cd19[['RNA']]@counts['CD27',cells])))
}

memory <- subset(cd19, ident = c(1,2,3,5,10,18,20,22,25))

clusters = levels(memory@active.ident)
naive.IgM <- for (c in clusters){
  cells <- WhichCells(memory, idents = c)
  print (median(as.numeric(memory[['ADT']]@counts['ADT-IGD-TotalSeqC',cells])))
}

TS <- subset(cd19, ident = c(14))
TSandnaive <- subset(cd19, ident = c(0,4,6,7,8,9,11,12,13,14,15,16,17,19,21,23,24))
naive <- subset(cd19, ident = c(0,4,6,7,8,9,11,12,13,15,16,17,19,21,23,24))

clusters = levels(naive@active.ident)
for (c in clusters){
  cells <- WhichCells(naive, idents = c)
  print (median(as.numeric(naive[['ADT']]@counts['ADT-IGM-TotalSeqC',cells])))
}

## Rename clusters, IGM high naive being clusters with top 30% Igm expression

cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(26)), value = "PB")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(14,12)), value = "TS")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(2,5)), value = "MZB")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(3,10,18,22,25)), value = "CSM")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(7,8,9,11,13,15,16,17,21,24)), value = "Naive")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(0,4,6,19,23)), value = "Naive Mhi")
cd19 <- SetIdent(cd19, cells = WhichCells(cd19, idents = c(1,20)), value = "M ONLY")

## Run 3D UMAP
cd19 <- RunUMAP(cd19, dims = 1:15, reduction.key = 'UMAPRNA_', reduction.name = 'UMAP_rna', n.components = 3L)
