library(Seurat)
library(ggplot2)
library(patchwork)

samples = c('Sample1', 'Sample2', 'Sample3')

# Load the data
qc_CD19_data <- list()
for(sample in samples){
  qc_CD19_data[[sample]] <- readRDS(paste('data/processed/checkpoints/quality_control_CD19_', sample, '.rds', sep=''))
}

normalised_CD19_data <- list()
# Normalise data with SCTransform
for (sample in samples){
  # Normalise RNA
  normalised_CD19_data[[sample]] <- SCTransform(qc_CD19_data[[sample]], verbose = TRUE, 
                                                  vars.to.regress = c('percent_mt'))
  
  # Normalise ADT
  normalised_CD19_data[[sample]] <- NormalizeData(normalised_CD19_data[[sample]], assay = "ADT", normalization.method = "CLR",
                                                  verbose = F)
}

saveRDS(normalised_CD19_data, "data/processed/normalised_data.rds")
