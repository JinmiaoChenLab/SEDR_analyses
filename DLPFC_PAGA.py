# Import libraries necessary for this project
from __future__ import division

import numpy as np
import json, os
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from numpy import arange
from time import time
from IPython.display import display # Allows the use of display() for DataFrames
pd.options.display.max_columns = None

import matplotlib.pyplot as pl
from matplotlib import rcParams

# setting visualization/logging parameters
pd.set_option('display.max_columns', None)
sc.set_figure_params(dpi=100, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_versions()

# Pretty display for notebooks
%matplotlib inline



ls_non_match_length = []
# list_samples = ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674', '151675', '151676']

# Reading clusters and embedding files
for sample in tqdm(list_samples):
    df_cluster = pd.read_csv(f'../data/DLPFC/{sample}/metadata.tsv', sep='\t')
    
    # for sed features
    feat = np.load(f'../output/DLPFC/{sample}/SEDR/SED_result.npz')['sed_feat']
    df_feat = pd.DataFrame(feat, 
                           index=df_cluster.index, 
                           columns=['E{}'.format(i) for i in range(0, feat.shape[1])])

#     # for seurat PCs
#     df_feat = pd.read_csv('{}/Seurat/seurat.PCs.tsv'.format(pth), 
#                            index_col=0, sep='\t')

    # Convert to adata and assign Embedding values as X_PCA
    adata = sc.AnnData(df_feat)
    adata.obsm['X_pca'] = df_feat.values
    sc.pp.neighbors(adata, n_pcs=len(df_feat.columns)) # n_pcs=30 for seurat pcs
    

    for col in df_cluster.columns:
        adata.obs['{}'.format(col)] = df_cluster['{}'.format(col)].values
        
    target_label = 'layer_guess'
    adata = adata[adata.obs[target_label].notna()] # dropping rows with nan in layer_guess column

    adata.obs['leiden'] = adata.obs[target_label].astype(str)
    sc.tl.paga(adata)
    sc.pl.paga(adata, fontsize=14, plot=True, edge_width_scale=0.4)
        
    df_paga_edgeWeight = pd.DataFrame(adata.uns['paga']['connectivities'].toarray(), 
                                          columns=sorted(adata.obs['leiden'].unique()),
                                          index=sorted(adata.obs['leiden'].unique()))
        
    print(df_paga_edgeWeight)
    
    df_dist = pd.DataFrame(data = adata.uns['neighbors']['distances'].toarray(),
             columns = adata.obs['layer_guess'].values,
             index = adata.obs['layer_guess'].values)




ls_non_match_length = []
# list_samples = ['151507', '151508', '151509', '151510', '151669', '151670', '151671', '151672', '151673', '151674', '151675', '151676']

# Reading clusters and embedding files
for sample in tqdm(list_samples):
    df_cluster = pd.read_csv(f'../data/DLPFC/{sample}/metadata.tsv', sep='\t')
    
    # for sed features
#     feat = np.load(f'../output/DLPFC/{sample}/Sed/SED_result.npz')['sed_feat']
#     df_feat = pd.DataFrame(feat, 
#                            index=df_cluster.index, 
#                            columns=['E{}'.format(i) for i in range(0, feat.shape[1])])

    # for seurat PCs
    df_feat = pd.read_csv(f'../output/DLPFC/{sample}/Seurat/seurat.PCs.tsv', 
                           index_col=0, sep='\t')

    # Convert to adata and assign Embedding values as X_PCA
    adata = sc.AnnData(df_feat)
    adata.obsm['X_pca'] = df_feat.values
    sc.pp.neighbors(adata, n_pcs=len(df_feat.columns)) # n_pcs=30 for seurat pcs
    

    for col in df_cluster.columns:
        adata.obs['{}'.format(col)] = df_cluster['{}'.format(col)].values
        
    target_label = 'layer_guess'
    adata = adata[adata.obs[target_label].notna()] # dropping rows with nan in layer_guess column

    adata.obs['leiden'] = adata.obs[target_label].astype(str)
    sc.tl.paga(adata)
    sc.pl.paga(adata, fontsize=14, plot=True, edge_width_scale=0.4)
        
    df_paga_edgeWeight = pd.DataFrame(adata.uns['paga']['connectivities'].toarray(), 
                                          columns=sorted(adata.obs['leiden'].unique()),
                                          index=sorted(adata.obs['leiden'].unique()))
        
    print(df_paga_edgeWeight)