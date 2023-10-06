import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys

from sklearn.metrics.cluster import adjusted_rand_score

import STAGATE


# the location of R (used for the mclust clustering)
os.environ['R_HOME'] = '/scbio4/tools/R/R-4.0.3_openblas/R-4.0.3'
os.environ['R_USER'] = '/home/xuhang/python/anaconda3/envs/STAGATE/lib/python3.7/site-packages/rpy2'




section_id = sys.argv[1]
n_clusters = sys.argv[2]


dir_out = f'./output/DLPFC/{section_id}/STAGATE'
os.makedirs(dir_out, exist_ok=True)


input_dir = os.path.join('./data/DLPFC/', section_id)
adata = sc.read_visium(path=input_dir, count_file='filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()


#Normalization
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# read the annotation
Ann_df = pd.read_csv(os.path.join('./data/DLPFC/', section_id, 'metadata.tsv'), sep='\t', index_col=0)
Ann_df['Ground Truth'] = Ann_df['layer_guess']


adata.obs['Ground Truth'] = Ann_df.loc[adata.obs_names, 'Ground Truth']

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, img_key="hires", color=["Ground Truth"])



## Constructing the spatial network
STAGATE.Cal_Spatial_Net(adata, rad_cutoff=150)
STAGATE.Stats_Spatial_Net(adata)


## Runing STAGATE
adata = STAGATE.train_STAGATE(adata, alpha=0)
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)
adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters)


obs_df = adata.obs.dropna()
ARI = adjusted_rand_score(obs_df['mclust'], obs_df['Ground Truth'])
print('Adjusted rand index = %.2f' %ARI)

plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.umap(adata, color=["mclust", "Ground Truth"], title=['STAGATE (ARI=%.2f)'%ARI, "Ground Truth"])


plt.rcParams["figure.figsize"] = (3, 3)
sc.pl.spatial(adata, color=["mclust", "Ground Truth"], title=['STAGATE (ARI=%.2f)'%ARI, "Ground Truth"])


adata.obs['STAGATE'] = adata.obs['mclust']
adata.write(f'{dir_out}/result.h5ad')

adata.obs.to_csv(f'{dir_out}/metadata.tsv', sep='\t')


df = pd.DataFrame(data=adata.obsm['STAGATE'], index=adata.obs.index)
df.to_csv(f'{dir_out}/PCs.tsv', sep='\t')


