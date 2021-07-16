args = commandArgs(trailingOnly=TRUE)
sample.name <- args[1]
n_clusters <- as.numeric(args[2])

library(BayesSpace)
library(ggplot2)

dir.input <- file.path('../data/DLPFC/', sample.name)
dir.output <- file.path('../output/DLPFC/', sample.name, '/BayesSpace/')

if(!dir.exists(file.path(dir.output))){
  dir.create(file.path(dir.output), recursive = TRUE)
}


dlpfc <- getRDS("2020_maynard_prefrontal-cortex", sample.name)

set.seed(101)
dec <- scran::modelGeneVar(dlpfc)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(102)
dlpfc <- scater::runPCA(dlpfc, subset_row=top)

## Add BayesSpace metadata
dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)


##### Clustering with BayesSpace
q <- n_clusters  # Number of clusters
d <- 15  # Number of PCs

## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=3, save.chain=TRUE)

labels <- dlpfc$spatial.cluster

## View results
clusterPlot(dlpfc, label=labels, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A", labels = 1:7) +
  labs(title="BayesSpace")

ggsave(file.path(dir.output, 'clusterPlot.png'), width=5, height=5)


write.table(colData(dlpfc), file=file.path(dir.output, 'metadata.tsv'), sep='\t', quote=FALSE)

