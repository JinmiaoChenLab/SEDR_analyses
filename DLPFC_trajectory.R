# setwd('~/Xuhang/Projects/spTrans/script')

##### Load library
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(monocle3)
library(SeuratWrappers)
options(bitmapType='cairo')


args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

sp_data <- readRDS(file.path('../output/', sample, '/Seurat/Seurat_final.rds'))

metadata <- read.table(file.path('../output/', sample, '/clustering_results.tsv'), sep='\t', header=T)

sp_data <- AddMetaData(sp_data,
                       metadata = metadata$layer_guess,
                       col.name = 'layer_guess')
sp_data <- AddMetaData(sp_data,
                       metadata = metadata$Sed,
                       col.name = 'Sed')


##### detect the first point for ordering
subset_sp <- subset(sp_data, subset = layer_guess == 'WM')


library(fields)
sub_umap <- subset_sp@reductions$umap@cell.embeddings
df <- data.frame(umap1=sub_umap[,1], umap2=sub_umap[,2])
rownames(df) <- row.names(sub_umap)

sub_mean <- c(mean(sub_umap[,1]), mean(sub_umap[,2]))
sub_dist <- rdist(as.matrix(df), t(as.matrix(sub_mean)))
row.names(sub_dist) <- row.names(sub_umap)
sorted_sub_dist <- sort.DataFrame(sub_dist)
root_rna <- rownames(sorted_sub_dist)[1]



##### RNA-monocle
DefaultAssay(sp_data) <- 'SCT'
rna_cds <- as.cell_data_set(sp_data)

rna_cds <- cluster_cells(rna_cds, reduction_method = 'UMAP')
rna_cds <- learn_graph(rna_cds, use_partition = FALSE)
rna_cds <- order_cells(rna_cds, root_cells = root_rna)


##### Sed-monocle
##### DimPlot 
library(reticulate) 
library(uwot) 
np <- import("numpy")

### load sed 
sed_feat <- np$load(file.path('../output/', sample, '/Sed/SED_result.npz'))$f[['sed_feat']]
rownames(sed_feat) <- row.names(sp_data@meta.data) 
sp_data@reductions$sed <- CreateDimReducObject(embeddings = sed_feat, key = 'sedDimred_')

sed_umap <- tumap(sed_feat, n_neighbors = 15, verbose = TRUE)
colnames(sed_umap) <- c('UMAP_1', 'UMAP_2') 
rownames(sed_umap) <-row.names(sp_data@meta.data) 
sp_data@reductions$sed_umap <-CreateDimReducObject(embeddings = sed_umap, key = 'sedUMAP_')

sed.sp_data <- sp_data
sed.sp_data@reductions$umap <-CreateDimReducObject(embeddings = sed_umap, key = 'UMAP_', assay = 'SCT')


sed_cds <- as.cell_data_set(sed.sp_data)
sed_cds <- cluster_cells(sed_cds, reduction_method = 'UMAP')
sed_cds <- learn_graph(sed_cds, use_partition = FALSE)
sed_cds <- order_cells(sed_cds, root_cells = root_rna)



##### extract pseudotime
rna_pseudotime <- pseudotime(rna_cds)
sp_data <- AddMetaData(sp_data,
                       metadata = rna_pseudotime,
                       col.name = 'rna_pseudotime')
sed_pseudotime <- pseudotime(sed_cds)
sp_data <- AddMetaData(sp_data,
                       metadata = sed_pseudotime,
                       col.name = 'sed_pseudotime')

df <- sp_data@images$slice1@coordinates
df$rna_pseudotime <- rna_pseudotime
df$sed_pseudotime <- sed_pseudotime

df$layer_guess <- sp_data@meta.data$layer_guess




p1 <- DimPlot(sp_data, group.by = 'layer_guess', reduction = 'umap', label=T, repel=T, pt.size = 0.1, label.size = 2.5) + ggtitle('') + theme_bw() + 
  NoLegend() + xlab('rnaUMAP_1') + ylab('rnaUMAP_2') + 
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=8, face='bold'),
        axis.text = element_text(size=6),
        panel.grid = element_blank())

p2 <- DimPlot(sp_data, group.by = 'layer_guess', reduction = 'sed_umap', label = T, repel=T, pt.size=0.1, label.size = 2.5) + ggtitle('') + theme_bw() + 
  NoLegend()+ xlab('sedUMAP_1') + ylab('sedUMAP_2') +
  theme(aspect.ratio = 1, 
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=8, face='bold'),
        axis.text = element_text(size=6),
        panel.grid = element_blank())

p3 <- plot_cells(rna_cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 label_roots = FALSE) + 
  ggtitle('') + theme_bw() + NoLegend() + xlab('rnaUMAP_1') + ylab('rnaUMAP_2')  + 
  theme(aspect.ratio=1, 
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=8, face='bold'),
        axis.text = element_text(size=6),
        panel.grid = element_blank())

p4 <- plot_cells(sed_cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5,
                 label_roots = FALSE) + 
  ggtitle('') + theme_bw() + NoLegend() + xlab('sedUMAP_1') + ylab('sedUMAP_2') +
  theme(aspect.ratio=1, 
        plot.title = element_text(hjust=0.5),
        axis.title = element_text(size=8, face='bold'),
        axis.text = element_text(size=6),
        panel.grid = element_blank())




p5 <- ggplot(df, aes(imagecol, imagerow, color=rna_pseudotime)) + geom_point(stroke=0, size=1.1) + 
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma") + 
  theme(legend.position="bottom") + 
  guides(color = guide_colourbar(ticks = FALSE, barheight = 0.5))

p6 <- ggplot(df, aes(imagecol, imagerow, color=sed_pseudotime)) + geom_point(stroke=0, size=1.1) + 
  coord_fixed() + scale_y_reverse() + theme_void() + viridis::scale_color_viridis(option="plasma") +
  theme(legend.position='bottom') + 
  guides(color = guide_colourbar(ticks = FALSE, barheight = 0.5))

P1 <- ((p1 / p3) | p5) 
P2 <- ((p2 / p4) | p6)
P1
P2
ggsave(file.path('../output/', sample, '/Trajectory_analyses.png'), width=16, height=13, units='cm')
ggsave(file.path('../output/', sample, '/Trajectory_analyses.pdf'), width=16, height=13, units='cm') 

ggsave(plot=P1, file.path('../output/', sample, '/Trajectory_analyses.1.pdf'), width=8, height=6, units='cm', dpi=300) 
ggsave(plot=P2, file.path('../output/', sample, '/Trajectory_analyses.2.pdf'), width=8, height=6, units='cm') 




ggsave(plot=p1, file.path('../output/', sample, '/test.RNA_UMAP.1.pdf'), width=5.5, height=5.5, units='cm')
ggsave(plot=p2, file.path('../output/', sample, '/test.Sed_UMAP.1.pdf'), width=5.5, height=5.5, units='cm')

ggsave(plot=p3, file.path('../output/', sample, '/test.RNA_UMAP.pseudotime.1.pdf'), width=5, height=5, units='cm')
ggsave(plot=p4, file.path('../output/', sample, '/test.Sed_UMAP.pseudotime.1.pdf'), width=5, height=5, units='cm')

ggsave(plot=p5, file.path('../output/', sample, '/test.RNA_spatial.pseudotime.1.pdf'), width=6.1, height=6.1, units='cm')
ggsave(plot=p6, file.path('../output/', sample, '/test.Sed_spatial.pseudotime.1.pdf'), width=6.1, height=6.1, units='cm')

