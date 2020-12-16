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

# Better* colour palette
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


# ---- SLINGSHOT ----
reduced_data_CD19 <- readRDS("data/processed/reduced.rds")

# Slingshot - manually-defined dimensions
# See Lavaert, M. et al (2020). Integrated scRNA-Seq Identifies Human Postnatal Thymus Seeding Progenitors and Regulatory Dynamics of 
# Differentiating Immature Thymocytes. Immunity, 1–17. https://doi.org/10.1016/j.immuni.2020.03.019

reduced_data_CD19 <- FindClusters(reduced_data_CD19, resolution = 0.8)

smaller_CD19 <- reduced_data_CD19

# As before
sce <- as.SingleCellExperiment(reduced_data_CD19)

# Now we specify which dimensions we want, namely 1st and 3rd columns of the UMAP3D_RNA
plot_data <- reducedDims(sce)$UMAP3D_RNA[,colnames(reducedDims(sce)$UMAP3D_RNA)[c(1,3)]]
plot(plot_data, col = getPalette(length(levels(reduced_data_CD19@active.ident)))[smaller_CD19@active.ident], pch=16, asp = 1)

# Create slingshot
slingshot_health_CD19 <- slingshot(sce, reduced_data_CD19@active.ident, reducedDim = 'UMAP3D_RNA', start.clus = '11')

# # Get lineages and curves
plot(plot_data, col = getPalette(length(levels(reduced_data_CD19@active.ident)))[reduced_data_CD19@active.ident], pch=16, asp = 1)
a <- SlingshotDataSet(slingshot_health_CD19)
a@curves <- list(SlingshotDataSet(slingshot_health_CD19)@curves$curve2)
lines(a, dims=c(1,3), lwd = 5)

curve <- as.data.frame(SlingshotDataSet(slingshot_health_CD19)@curves$curve2$s[SlingshotDataSet(slingshot_health_CD19)@curves$curve2$ord,])

# Get the 3D curves
reduced_data_CD19@misc$curves <- list()

# Save the curves to the misc section of the seurat object
for (curvename in ls(SlingshotDataSet(slingshot_health_CD19)@curves)){
  curve <- SlingshotDataSet(slingshot_health_CD19)@curves[[curvename]]
  reduced_data_CD19@misc$curves[[curvename]] <- as.data.frame(curve$s[curve$ord,])
}

# Save seurat object for visuals
saveRDS(reduced_data_CD19, "data/processed/reduced_with_curves.rds")
saveRDS(slingshot_health_CD19, "data/processed/slingshot_object.rds")
saveRDS(smaller_CD19, "data/processed/smaller_CD19.rds")

reduced_data_CD19 <- readRDS("data/processed/reduced_with_curves.rds")
slingshot_health_CD19 <- readRDS("data/processed/slingshot_object.rds")
smaller_CD19 <- readRDS("data/processed/smaller_CD19.rds")

sds <- SlingshotDataSet(slingshot_health_CD19)
t <- slingshot_health_CD19$slingPseudotime_2

g1 = c('IGLL5', 'VPREB1')
g2 = c('PLD4', 'SOX4', 'MZB1', 'CD1C')
g3 = c('CD24', 'CD27', 'CD82', 'TNFRSF13B', 'COTL1', 'CLECL1', 'ITGB2')
g4 = c('PLAC8', 'FCER2', 'DAPP1', 'FCRL5', 'MS4A1')
g5 = c('PCDH9', 'EZR', 'S100A11', 'S100A6', 'S100A4', 'CIB1', 'NFKBIA', 'CD69', 'EMP3', 'RALGPS2', 'SAT1', 'LY6E' )
g6 = c('MX1', 'CD74')
g7 = c('ZEB2', 'JUN')

genes = c(g1, g2, g3, g4, g5, g6, g7)

heatdata <- smaller_CD19@assays$integrated@data[genes, order(t, na.last=NA)]

DoHeatmap(smaller_CD19, group.by='test', group.bar = F, disp.max = 4.5, features=genes, cells= colnames(heatdata), assay = 'integrated')

q <- t(rollmean(t(as.matrix(heatdata)), k = 100, align = 'r'))
pheatmap(q, color = viridis(100), cluster_rows = F, cluster_cols = F)

draw_trajectory_heatmap(x = t(as.matrix(heatdata)), time = t[order(t, na.last = NA)], col=viridis(100))

for (i in 1:4){
  f <- (i-1)*9 + 1
  j <- f+8
  if ( j > length(genes)){
    j <- length(genes)
  }
  print(FeaturePlot(reduced_data_CD19, genes[f:j], dims = c(1, 3), max.cutoff = 'q95', order = T))
}

adts <- FeaturePlot(reduced_data_CD19, c('ADT-IGM-TotalSeqC', 'ADT-CD27-TotalSeqC', 'ADT-CD38-TotalSeqC', 'ADT-IGD-TotalSeqC'), dims = c(1, 3), 
                  max.cutoff = 'q95', ncol = 4, combine=F)
for(i in 1:length(adts)) {
  adts[[i]] <- adts[[i]] + NoLegend() + xlab("")  + ylab("") + scale_color_gradientn(colours = c(rgb(224, 224, 224, max = 255), rgb(0, 0, 255, max = 255)))
}
adts[[1]] <- adts[[1]] + ggtitle('ADT IGM')
adts[[2]] <- adts[[2]] + ggtitle('ADT CD27')
adts[[3]] <- adts[[3]] + ggtitle('ADT CD38')
adts[[4]] <- adts[[4]] + ggtitle('ADT IGD')

(adts[[1]] | adts[[2]]) / (adts[[3]] | adts[[4]])

ps <- FeaturePlot(reduced_data_CD19, c('CD1C', 'PLD4', 'SOX4', 'ITGB2', 'PCDH9', 'DAPP1', 'MX1', 'ZEB2'), dims = c(1, 3), 
                  max.cutoff = 'q95', order=T, ncol = 4, combine=F)
for(i in 1:length(ps)) {
  ps[[i]] <- ps[[i]] + NoLegend() + xlab("")  + ylab("") + scale_color_gradientn(colours = c(rgb(224, 224, 224, max = 255), rgb(255, 0, 0, max = 255)))
}
(ps[[1]] | ps[[2]] | ps[[3]] | ps[[4]]) / (ps[[5]] | ps[[6]] | ps[[7]] | ps[[8]])



other_genes = c('IGLL5', 'VPREB1', 'MZB1', 'CD24', 'CD27', 'CD82', 'TNFRSF13B', 'COTL1', 'CLECL1', 'PLAC8', 'FCER2', 'FCRL5', 'MS4A1', 'EZR', 'S100A11', 'S100A6', 
                'S100A4', 'CIB1', 'NFKBIA', 'CD69', 'EMP3', 'RALGPS2', 'SAT1', 'LY6E', 'CD74', 'JUN')

ops <- FeaturePlot(reduced_data_CD19, other_genes, dims = c(1, 3), 
                  max.cutoff = 'q95', order=T, ncol = 4, combine=F)

for(i in 1:length(ops)) {
  ops[[i]] <- ops[[i]] + NoLegend() + xlab("")  + ylab("") + scale_color_gradientn(colours = c(rgb(224, 224, 224, max = 255), rgb(255, 0, 0, max = 255)))
}

cowplot::plot_grid(plotlist = ops)

heatdata <- smaller_CD19@assays$SCT@counts[topgenes, order(t, na.last=NA)]
heatclus <- smaller_CD19@active.ident[order(t, na.last = NA)]

heatmap(as.matrix(heatdata),
        Colv = NA,
        Rowv = NA,
        ColSideColors = getPalette(length(levels(smaller_CD19@active.ident)))[heatclus])

smaller_CD19[['test']] <- 'A'
DoHeatmap(smaller_CD19, group.by='test', features=topgenes, cells= colnames(heatdata), assay = 'SCT')
