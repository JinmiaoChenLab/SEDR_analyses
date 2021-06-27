setwd('~/Xuhang/Projects/spTrans/script/')
library(ggplot2)
library(ggpubr)
library(patchwork)
options(bitmapType = 'cairo')


list_sample <- c('151507','151508','151509','151510','151673','151674','151675','151676')
c_method <- c()
c_score <- c()
for (sample in list_sample){
    df1 <- read.table(file.path('../output/DLPFC/', sample, '/PAGA/SEDR.PAGA_wieght.layer_guess.csv'), sep=',', header=T, row.names = 1)
    df2 <- read.table(file.path('../output/DLPFC/', sample, '/PAGA/seurat.PAGA_wieght.layer_guess.csv'), sep=',', header=T, row.names = 1)

    
    right_connect <- df1['WM', 'Layer6'] + df1['Layer6', 'Layer5'] +  df1['Layer5', 'Layer4'] +  
      df1['Layer4', 'Layer3'] + df1['Layer3', 'Layer2'] +  df1['Layer2', 'Layer1']

    a <- as.matrix(df1)
    total_connect <- sum(a[which(upper.tri(a))])

    score1 <- right_connect/(total_connect)
    c_method <- c(c_method, 'SEDR')
    c_score <- c(c_score, score1)

    right_connect <- df2['WM', 'Layer6'] + df2['Layer6', 'Layer5'] +  df2['Layer5', 'Layer4'] +  
      df2['Layer4', 'Layer3'] + df2['Layer3', 'Layer2'] +  df2['Layer2', 'Layer1']

    a <- as.matrix(df2)
    total_connect <- sum(a[which(upper.tri(a))])
    score2 <- right_connect/(total_connect)
    c_method <- c(c_method, 'Seurat')
    c_score <- c(c_score, score2)
    
    print(score1)
    print(score2)
}

list_sample <- c('151669','151670','151671','151672')
for (sample in list_sample){
  df1 <- read.table(file.path('../output/DLPFC/', sample, '/PAGA/SEDR.PAGA_wieght.layer_guess.csv'), sep=',', header=T, row.names = 1)
  df2 <- read.table(file.path('../output/DLPFC/', sample, '/PAGA/seurat.PAGA_wieght.layer_guess.csv'), sep=',', header=T, row.names = 1)
  
  
  right_connect <- df1['WM', 'Layer6'] + df1['Layer6', 'Layer5'] +  df1['Layer5', 'Layer4'] +  
    df1['Layer4', 'Layer3']
  
  a <- as.matrix(df1)
  total_connect <- sum(a[which(upper.tri(a))])
  
  score1 <- right_connect/(total_connect)
  c_method <- c(c_method, 'SEDR')
  c_score <- c(c_score, score1)
  
  right_connect <- df2['WM', 'Layer6'] + df2['Layer6', 'Layer5'] +  df2['Layer5', 'Layer4'] +  
    df2['Layer4', 'Layer3']
  
  a <- as.matrix(df2)
  total_connect <- sum(a[which(upper.tri(a))])
  score2 <- right_connect/(total_connect)
  c_method <- c(c_method, 'Seurat')
  c_score <- c(c_score, score2)
  
  print(score1)
  print(score2)
}


df <- data.frame(method=c_method, weight_percentage=c_score)

ggplot(df, aes(method, weight_percentage)) + geom_boxplot() + theme_bw() + theme(panel.grid = element_blank())

my_comparisons <- list(c('SEDR', 'Seurat'))
ggboxplot(df, x='method', y='weight_percentage', width = 0.5) + 
  stat_compare_means(comparisons = my_comparisons) + 
  theme(axis.title.y = element_text(size=14, face = 'bold'),
        axis.title.x = element_blank(),
        axis.text = element_text(size=12) )

ggsave('../output/DLPFC/Total/PAGA_weight.comparison.png', width=3.5, height=3.5)
ggsave('../output/DLPFC/Total/PAGA_weight.comparison.pdf', width=3.5, height=3.5)

