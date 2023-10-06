library(dplyr)
library(Giotto)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(ggpubr)
library(mclust)
options(bitmapType = 'cairo')

list.samples <- c("151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676")
list.methods <- c( "Seurat", "Giotto", "stLearn", "SpaGCN", "BayesSpace", "DeepST", "STAGATE", "SEDR")

##### Generate data
c1 <- c()
c2 <- c()
c3 <- c()

for (sample in list.samples) {
  file.results <- file.path('./output/DLPFC/', sample, '/Comparison/comparison.tsv')
  df.results <- read.table(file.results, sep='\t', header=T)
  for (method in list.methods){
    cluster <- df.results  %>% dplyr::select(c(method))
    ARI <- adjustedRandIndex(x = df.results$layer_guess, y = cluster[, 1])

    c1 <- c(c1, method)
    c2 <- c(c2, sample)
    c3 <- c(c3, ARI)
  }
}

df.comp <- data.frame(method = c1,
                      sample = c2,
                      ARI = c3)


##### Plot results
df.comp$method <- as.factor(df.comp$method)
df.comp$method <- factor(df.comp$method,
                         levels = c('Seurat', 'Giotto', 'stLearn', 'SpaGCN', 'BayesSpace','DeepST','STAGATE','SEDR'))


ggplot(df.comp, aes(method, ARI)) +
  geom_boxplot(width=0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, size=0.5) +
  coord_flip() +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=14)
        )

# ggsave('./output/DLPFC/All/ARI.violin_and_jitter.png', width=4, height=4)
# ggsave('./output/DLPFC/All/ARI.violin_and_jitter.pdf', width=4, height=4)
#
