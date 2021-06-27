# setwd('~/Xuhang/Projects/spTrans/script/')

library(mclust)
library(ggplot2)
library(patchwork)
library(Seurat)
library(mclust)
options(bitmapType = 'cairo')

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

sp_data <- readRDS(file.path('../output/DLPFC/', sample, '/Seurat/Seurat_final.rds'))

##### SpatialDimPlot
metadata <- read.table(file.path('../data/DLPFC/', sample, 'metadata.tsv'), sep='\t', header=TRUE)

sed_cluster <- read.table(file.path('../output/DLPFC/', sample, 'SEDR/metadata.tsv'), sep='\t', header=TRUE, row.names =1)
sed_cluster$sed_labels <- sed_cluster$leiden_fixed_clusCount

seurat_cluster <- read.table(file.path('../output/DLPFC/', sample, '/Seurat/metadata.tsv'), sep='\t', header=TRUE)
spaGCN_cluster <- read.table(file.path('../output/DLPFC/', sample, '/SpaGCN/metadata.tsv'), sep='\t', header=TRUE)
BayesSpace_cluster <- read.table(file.path('../output/DLPFC/', sample, '/BayesSpace/metadata.tsv'), sep='\t', header=TRUE)
Giotto_cluster <- read.table(file.path('../output/DLPFC/', sample, '/Giotto/meta.tsv'), sep='\t', header=TRUE)
row.names(Giotto_cluster) <- Giotto_cluster$cell_ID
stLearn_cluster <- read.table(file.path('../output/DLPFC/', sample, '/stLearn/metadata.tsv'), sep='\t', header=TRUE)

truth <- as.factor(metadata$layer_guess)
truth <- factor(truth, levels=c('WM', 'nan', 'Layer6', 'Layer5', 'Layer4', 'Layer3', 'Layer2', 'Layer1'))
sp_data <- AddMetaData(sp_data, truth, col.name = 'layer_guess')
sp_data <- AddMetaData(sp_data, sed_cluster[, 'sed_labels', drop=F], col.name = 'SEDR')
sp_data <- AddMetaData(sp_data, seurat_cluster$seurat_clusters, col.name = 'Seurat')
sp_data <- AddMetaData(sp_data, spaGCN_cluster$refined_pred, col.name = 'SpaGCN')
sp_data <- AddMetaData(sp_data, BayesSpace_cluster$spatial.cluster, col.name = 'BayesSpace')
# sp_data <- AddMetaData(sp_data, Giotto_cluster$leiden_clus, col.name = 'Giotto')
sp_data <- AddMetaData(sp_data, Giotto_cluster[, 'HMRF_cluster', drop=F], col.name = 'Giotto')
sp_data <- AddMetaData(sp_data, stLearn_cluster$X_pca_kmeans, col.name = 'stLearn')

SEDR_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$SEDR)
Seurat_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$Seurat)
SpaGCN_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$SpaGCN)
BayesSpace_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$BayesSpace)
Giotto_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$Giotto)
stLearn_ARI = adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$stLearn)

df_clusters <- data.frame(layer_guess = sp_data@meta.data$layer_guess,
                          SEDR = as.factor(sp_data@meta.data$SEDR),
                          Seurat = as.factor(sp_data@meta.data$Seurat),
                          SpaGCN = as.factor(sp_data@meta.data$SpaGCN),
                          BayesSpace = as.factor(sp_data@meta.data$BayesSpace),
                          Giotto = as.factor(sp_data@meta.data$Giotto),
                          stLearn = as.factor(sp_data@meta.data$stLearn))

df <- sp_data@images$slice1@coordinates
df <- cbind(df, df_clusters)

p0 <- ggplot(df, aes(imagecol, imagerow, color=layer_guess)) + geom_point(stroke=0, size=1.1) + ggtitle('layer_guess') +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p1 <- ggplot(df, aes(imagecol, imagerow, color=Seurat)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('Seurat: ARI=', round(Seurat_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p2 <- ggplot(df, aes(imagecol, imagerow, color=Giotto)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('Giotto: ARI=', round(Giotto_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p3 <- ggplot(df, aes(imagecol, imagerow, color=stLearn)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('stLearn: ARI=', round(stLearn_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p4 <- ggplot(df, aes(imagecol, imagerow, color=SpaGCN)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('SpaGCN: ARI=', round(SpaGCN_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p5 <- ggplot(df, aes(imagecol, imagerow, color=BayesSpace)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('BayesSpace: ARI=', round(BayesSpace_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p6 <- ggplot(df, aes(imagecol, imagerow, color=SEDR)) + geom_point(stroke=0, size=1.1) + ggtitle(paste('Sed: ARI=', round(SEDR_ARI, 3))) +
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma", discrete = T)

p0 + p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 4, widths = c(1,1,1,1), heights = c(1,1)) & NoLegend()


ggsave(filename = file.path('../output/DLPFC/', sample, '/comparison.png'), width=11, height=5.5)
ggsave(filename = file.path('../output/DLPFC/', sample, '/comparison.pdf'), width=11, height=5.5)

write.table(df, file.path('../output/DLPFC/', sample, '/comparison.tsv'), sep='\t', quote=FALSE)

