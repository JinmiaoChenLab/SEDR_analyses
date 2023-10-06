args = commandArgs(trailingOnly=TRUE)
sample.name <- args[1]
n_clusters <- args[2]


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
options(bitmapType = 'cairo')


dir.input <- file.path('./data/DLPFC/', sample.name)
dir.output <- file.path('./output/DLPFC/', sample.name, '/Seurat/')

if(!dir.exists(file.path(dir.output))){
  dir.create(file.path(dir.output), recursive = TRUE)
}

### load data
sp_data <- Load10X_Spatial(dir.input, filename = "filtered_feature_bc_matrix.h5")

df_meta <- read.table(file.path(dir.input, 'metadata.tsv'))

sp_data <- AddMetaData(sp_data, 
                       metadata = df_meta$layer_guess,
                       col.name = 'layer_guess')

### Data processing
plot1 <- VlnPlot(sp_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(sp_data, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
ggsave(file.path(dir.output, './Seurat.QC.png'), width = 10, height=5)

# sctransform
sp_data <- SCTransform(sp_data, assay = "Spatial", verbose = FALSE)


### Dimensionality reduction, clustering, and visualization
sp_data <- RunPCA(sp_data, assay = "SCT", verbose = FALSE, npcs = 50)
sp_data <- FindNeighbors(sp_data, reduction = "pca", dims = 1:30)

for(resolution in 50:30){
  sp_data <- FindClusters(sp_data, verbose = F, resolution = resolution/100)
  if(length(levels(sp_data@meta.data$seurat_clusters)) == n_clusters){
    break
  }
}
sp_data <- FindClusters(sp_data, verbose = FALSE, resolution = 0.46)
sp_data <- RunUMAP(sp_data, reduction = "pca", dims = 1:30)

p1 <- DimPlot(sp_data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(sp_data, label = TRUE, label.size = 3)
p1 + p2
ggsave(file.path(dir.output, './Seurat.cell_cluster.png'), width=10, height=5)


##### save data
saveRDS(sp_data, file.path(dir.output, 'Seurat_final.rds'))

write.table(sp_data@reductions$pca@cell.embeddings, file = file.path(dir.output, 'seurat.PCs.tsv'), sep='\t', quote=F)

write.table(sp_data@meta.data, file = file.path(dir.output, './metadata.tsv'), sep='\t', quote=FALSE)


##### 
library(mclust)

print(adjustedRandIndex(sp_data@meta.data$layer_guess, sp_data@meta.data$seurat_clusters))