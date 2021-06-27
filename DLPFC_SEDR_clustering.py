import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
import sys


sample  = sys.argv[1]
n_cluster = int(sys.argv[2])

def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]
        
        return:
            resolution[int]
    '''
    for res in sorted(list(np.arange(0.2, 2.5, increment)), reverse=True):
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        if count_unique_leiden == fixed_clus_count:
            break
    return res



df_cluster = pd.read_csv(f'../data/DLPFC/{sample}/metadata.tsv', sep='\t')


feat = np.load(f'../output/DLPFC/{sample}/SEDR/SED_result.npz')['sed_feat']

df_feat = pd.DataFrame(feat, 
                        index=df_cluster.index, 
                        columns=['E{}'.format(i) for i in range(0, feat.shape[1])])

adata = sc.AnnData(df_feat)

adata.obsm['X_pca'] = df_feat.values



sc.pp.neighbors(adata, n_pcs=len(df_feat.columns))

adata.obs['layer_guess'] = df_cluster['layer_guess']
fixed_clus_count = n_cluster
res = res_search_fixed_clus(adata, fixed_clus_count=fixed_clus_count, increment=0.05)

sc.tl.leiden(adata, random_state=0, resolution=res)

df_cluster_8clus = pd.DataFrame(adata.obs[['layer_guess', 'leiden']])
df_cluster_8clus.columns = ['layer_guess', 'leiden_fixed_clusCount']

df_cluster_8clus.to_csv(f'../output/DLPFC/{sample}/SEDR/metadata.tsv', sep='\t')


df_cluster_8clus = df_cluster_8clus[~pd.isnull(df_cluster_8clus['layer_guess'])]
ARI= adjusted_rand_score(df_cluster_8clus['layer_guess'], df_cluster_8clus['leiden_fixed_clusCount'])
print(ARI)
