library(Seurat)
library(ggplot2)
library(destiny)
library(patchwork)
# Necessary cofig change to allow integration
options(future.globals.maxSize = 4000 * 1024^2)

# GLOBAL VARS
# List of sample names
samples = c('Sample1', 'Sample2', 'Sample3')

# -------------------------------------- 1. LOAD DATA -------------------------------------------------
raw_data <- list()
# Load the raw count matrices
for (sample in samples){
  raw_data[[sample]] <- Read10X(data.dir = paste("data/MZBTrajectory/", sample, "/filtered_feature_bc_matrix/", sep=''))
}

QC_nCell_min = 3

# Initialize the Seurat objects with the raw (non-normalized data).
seurat_data <- list()
for (sample in samples){
  # TODO - currently allowing all features, review
  seurat_data[[sample]] <- CreateSeuratObject(counts = raw_data[[sample]][['Gene Expression']], project = sample, 
                                              min.cells = QC_nCell_min)
  # Add antibody capture
  seurat_data[[sample]][["ADT"]] <- CreateAssayObject(raw_data[[sample]][["Antibody Capture"]][, colnames(x = seurat_data[[sample]])])
  cat(sample, '\n', sep='')
  cat(seurat_data[[sample]]@assays[["RNA"]]@counts@Dim[2], ' cells \n', sep='')
  cat(seurat_data[[sample]]@assays[["RNA"]]@counts@Dim[1], ' genes (after gene QC) \n', sep='')
  cat(seurat_data[[sample]]@assays[["ADT"]]@counts@Dim[1], ' antibodies \n', sep='')
  
  # Add the sample name as a meta data column (for identification downstream when data is integrated)
  seurat_data[[sample]][['Sample']] <- sample
}

# --------------------------------------- 1.1 CLASSIFY CD19 vs CD10 -------------------------------------------------
# Split the data  - CD19 vs CD10
# Create scatter plots of HASH1 vs HASH2
for (sample in samples){
  DefaultAssay(seurat_data[[sample]]) <- 'ADT'
  print(FeatureScatter(seurat_data[[sample]], feature1 = 'HTO-AHH2-TotalSeqC', feature2 = 'HTO-AHH1-TotalSeqC'))
}
# Set threshold values. Will end up with 4 populations: CD19 (high hash1, low hash2), CD10 (low hash1, high hash2),
# Doublets (high hash1, high hash2) and Uncertain (low hash1, low hash2).
HASH1_threshold <- list(
  'health_1' = 30,
  'health_2' = 100,
  'health_3' = 120
)
HASH2_threshold <- list(
  'health_1' = 3,
  'health_2' = 30,
  'health_3' = 25
)
# Now plot with the gates
for (sample in samples){
  print(FeatureScatter(seurat_data[[sample]], feature1 = 'HTO-AHH2-TotalSeqC', feature2 = 'HTO-AHH1-TotalSeqC') + 
          geom_vline(xintercept = HASH2_threshold[[sample]]) + geom_hline(yintercept = HASH1_threshold[[sample]]))
}

CD19_data <- list()
CD10_data <- list()
for(sample in samples){
  feature_data_AHH1 <- FetchData(object = seurat_data[[sample]], vars = c('HTO-AHH1-TotalSeqC'), slot='counts')
  feature_data_AHH2 <- FetchData(object = seurat_data[[sample]], vars = c('HTO-AHH2-TotalSeqC'), slot='counts')
  # Set CD19
  CD19_cell_names =  colnames(seurat_data[[sample]])[which(x = feature_data_AHH1 > HASH1_threshold[[sample]] &
                                                               feature_data_AHH2 <= HASH2_threshold[[sample]])]
  CD19_data[[sample]] <- subset(seurat_data[[sample]], cells = CD19_cell_names)
  
  # Set CD10
  CD10_cell_names =  colnames(seurat_data[[sample]])[which(x = feature_data_AHH1 <= HASH1_threshold[[sample]] &
                                                               feature_data_AHH2 > HASH2_threshold[[sample]])]
  if (length(CD10_cell_names) > 0){
    CD10_data[[sample]] <- subset(seurat_data[[sample]], cells = CD10_cell_names)
  }
}

# CHECKPOINT - SAVE
for (sample in samples){
  saveRDS(CD19_data[[sample]], paste('data/processed/checkpoints/1_load_data_CD19_', sample, '.rds', sep=''))
  saveRDS(CD10_data[[sample]], paste('data/processed/checkpoints/1_load_data_CD10_', sample, '.rds', sep=''))
}
# TODO: CD10 from this point onwards

# # ------------------------------------ 2. QUALITY CONTROL -------------------------------------------------
# CHECKPOINT - LOAD
bcell_genes <- c("CD79A", "CD79B", "CD19", "MS4A1")

for (sample in samples){
  DefaultAssay(CD19_data[[sample]]) <- 'RNA'
  # Add mitochondrial percentage
  CD19_data[[sample]][["percent_mt"]] <- PercentageFeatureSet(object = CD19_data[[sample]], pattern = "^MT-")
  # Add B cell gene percentage
  CD19_data[[sample]][["percent_Bcell"]] <- PercentageFeatureSet(object = CD19_data[[sample]], features = bcell_genes)
}

qc_CD19_data <- list()
# # -------------------------------------- 2.1 GENE QC --------------------------------------------------------
# Remove heavy chain genes & Hashes
ig_hl_chain_gene_patterns <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^IGKC", "^JCHAIN", "^AC233755")
# ig_hl_chain_gene_patterns <- c("^IGHV", "^IGLV", "^IGLC", "^IGKV", "^IGKC", "^AC233755")
hashes <-  c("^HTO-")

for (sample in samples){
  # Save JCHAIN and HTO hashes into the meta.data for future reference
  CD19_data[[sample]][['JCHAIN_counts']] <- as.numeric(CD19_data[[sample]]$RNA['JCHAIN',])
  CD19_data[[sample]][['HTO_1_counts']] <- as.numeric(CD19_data[[sample]]$ADT['HTO-AHH1-TotalSeqC',])
  CD19_data[[sample]][['HTO_2_counts']] <- as.numeric(CD19_data[[sample]]$ADT['HTO-AHH2-TotalSeqC',])
  
  features_to_remove <- list()
  # Get all features (genes and antiboides)
  features_to_keep <- c(rownames(CD19_data[[sample]]@assays$ADT), rownames(CD19_data[[sample]]@assays$RNA))
  # Loop through each pattern and remove any features which match the pattern
  for (f in c(ig_hl_chain_gene_patterns, hashes)){
    features_to_remove <- c(features_to_remove, features_to_keep[grepl(paste0(f, collapse = "|"), features_to_keep)])
    features_to_keep <- features_to_keep[!grepl(paste0(f, collapse = "|"), features_to_keep)]
  }
  # Subset to just the kept features
  qc_CD19_data[[sample]] <- subset(CD19_data[[sample]], features = features_to_keep)
}

# # -------------------------------------- 2.2 CELL QC --------------------------------------------------------
# QC thresholds
# QC - Minimum number of unique feature counts
QC_nFeature_min = list(
  'health_1' = 450,
  'health_2' = 250,
  'health_3' = 300
)
# QC - Maximum number of unique feature counts
QC_nFeature_max = list(
  'health_1' = 2000,
  'health_2' = 1500,
  'health_3' = 1500
)

# QC - Minimum number of counts
QC_nCount_min = list(
  'health_1'=0,
  'health_2'=0,
  'health_3'=0
)
# QC - Maximum number of counts
QC_nCount_max = list(
  'health_1'=10000,
  'health_2'=7500,
  'health_3'=7500
)

# QC - Maximum percentage of mitochondrial reads
QC_mt_max = list(
  'health_1' = 30,
  'health_2' = 30,
  'health_3' = 30
)

# QC - Minimum percentage of b cell genes
QC_bcell_min = list(
  'health_1' = 0.01,
  'health_2' = 0.01,
  'health_3' = 0.01
)

# Apply the QC thresholds
for (sample in samples){
  cat(sample, "\n", dim(qc_CD19_data[[sample]])[2], " cells before QC\n", sep='')
  # Get threshold values
  nfmin <- QC_nFeature_min[[sample]]
  nfmax <- QC_nFeature_max[[sample]]
  ncmin <- QC_nCount_min[[sample]]
  ncmax <- QC_nCount_max[[sample]]
  mt <- QC_mt_max[[sample]]
  bc <- QC_bcell_min[[sample]]
  # Apply to data
  qc_CD19_data[[sample]] <- subset(qc_CD19_data[[sample]], 
                                   subset = nFeature_RNA > nfmin & nFeature_RNA < nfmax &
                                     nCount_RNA > ncmin & nCount_RNA < ncmax & 
                                     percent_mt < mt &
                                     percent_Bcell > bc)
  cat(dim(qc_CD19_data[[sample]])[2], " cells after QC\n", sep='')
}

# CHECKPOINT - SAVE
for(sample in samples){
  saveRDS(qc_CD19_data[[sample]], paste('data/processed/quality_control_CD19_', sample, '.rds', sep=''))
}
