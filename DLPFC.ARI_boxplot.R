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
list.methods <- c( "Seurat", "Giotto", "stLearn", "SpaGCN", "BayesSpace", "SEDR")

# list.methods <- c( "Seurat", "stLearn", "SpaGCN", "BayesSpace", "SEDR")

##### Generate data
c1 <- c()
c2 <- c()
c3 <- c()

for (sample in list.samples) {
  file.results <- file.path('../output/DLPFC/', sample, '/comparison.tsv')
  df.results <- read.table(file.results, sep='\t', header=T)
  for (method in list.methods){
    cluster <- df.results  %>% select(c(method))
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
                         levels = c('Seurat', 'Giotto', 'stLearn', 'SpaGCN', 'BayesSpace', 'SEDR'))


print(mean(subset(df.comp, method == 'SEDR')$ARI))
print(mean(subset(df.comp, method == 'BayesSpace')$ARI))



my_comparisons <- list(c('SEDR', 'BayesSpace'), c('SEDR', 'SpaGCN'), c('SEDR', 'Giotto'), 
                       c('SEDR', 'stLearn'), c('SEDR', 'Seurat'))

compare_means(ARI ~ method,  data = df.comp)

ggboxplot(df.comp, x='method', y='ARI', color='method', palette='jco') + 
  stat_compare_means(comparisons = my_comparisons) + 
  stat_compare_means(labely = 50)

ggsave('../output/DLPFC/Total/pvalue.ARI_boxplot.png', width=6, height=5.5)
ggsave('../output/DLPFC/Total/pvalue.ARI_boxplot.pdf', width=6, height=5.5)


ggboxplot(df.comp, x='method', y='ARI', color='method', palette='jco') + 
  stat_compare_means(label = "p.signif", method = "wilcox.test", 
                     comparisons = my_comparisons, ref.group = ".all.")      

ggsave('../output/DLPFC/Total/sign.ARI_boxplot.png', width=6, height=5.5)
ggsave('../output/DLPFC/Total/sign.ARI_boxplot.pdf', width=6, height=5.5)



ggplot(df.comp, aes(method, ARI)) + 
  geom_boxplot(width=0.5) + 
  geom_jitter(width = 0.1, size=1) +
  coord_flip() +
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text = element_text(size=14)
        )

ggsave('../output/DLPFC/Total/ARI_violin.png', width=4, height=4)
ggsave('../output/DLPFC/Total/ARI_violin.pdf', width=4, height=4)


