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

source("code/convert_seurat_to_monocle.R")
options(future.globals.maxSize = 4000 * 1024^2)

# ----- INTEGRATE ------
normalised_CD19_data <- readRDS("data/processed/normalised.rds")

integrate_features <- SelectIntegrationFeatures(object.list = normalised_CD19_data, nfeatures = 3000)
integrated_data_CD19 <- PrepSCTIntegration(object.list = normalised_CD19_data, anchor.features = integrate_features,
                                           verbose = FALSE, )

integrate_anchors <- FindIntegrationAnchors(object.list = integrated_data_CD19, normalization.method = "SCT",
                                            anchor.features = integrate_features, verbose = FALSE)
integrated_data_CD19 <- IntegrateData(anchorset = integrate_anchors, normalization.method = "SCT",
                                      verbose = FALSE)

saveRDS(integrated_data_CD19, "data/processed/integrated.rds")
