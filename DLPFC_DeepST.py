import os
from DeepST import run
import matplotlib.pyplot as plt
from pathlib import Path
import scanpy as sc
import pandas as pd

data_path = "./data/DLPFC/" #### to your path

sample = sys.argv[1]
data_name = sample

save_root = Path(f'./output/DLPFC/{data_name}/DeepST')
os.makedirs(save_root, exist_ok=True)

# n_domains = 5 if data_name in ['151669','151670','151671','151672'] else 7 ###### the number of spatial domains.
n_domains = sys.argv[2]

deepen = run(
    save_path = save_root,
    task = "Identify_Domain", #### DeepST includes two tasks, one is "Identify_Domain" and the other is "Integration"
    pre_epochs = 800, ####  choose the number of training
    epochs = 1000, #### choose the number of training
    use_gpu = True)


###### Read in 10x Visium data, or user can read in themselves.
adata = deepen._get_adata(platform="Visium", data_path=data_path, data_name=data_name)

###### Segment the Morphological Image
adata = deepen._get_image_crop(adata, data_name=data_name)

###### Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress"
###### is only applicable to 10x visium and the remaining omics selects the other two.
###### "use_morphological" defines whether to use morphological images.
adata = deepen._get_augment(adata, spatial_type="LinearRegress", use_morphological=True)


###### Build graphs. "distType" includes "KDTree", "BallTree", "kneighbors_graph", "Radius", etc., see adj.py
graph_dict = deepen._get_graph(adata.obsm["spatial"], distType = "BallTree")

###### Enhanced data preprocessing
data = deepen._data_process(adata, pca_n_comps = 200)

###### Training models
deepst_embed = deepen._fit(
    data = data,
    graph_dict = graph_dict,
)


###### DeepST outputs
adata.obsm["DeepST_embed"] = deepst_embed

###### Define the number of space domains, and the model can also be customized. If it is a model custom priori = False.
adata = deepen._get_cluster_data(adata, n_domains=n_domains, priori = True)


###### Spatial localization map of the spatial domain
sc.pl.spatial(adata, color='DeepST_refine_domain', frameon = False, spot_size=150)
# plt.savefig(os.path.join(save_root, f'{data_name}_domains.pdf'), bbox_inches='tight', dpi=300)

adata.obs['DeepST'] = adata.obs['DeepST_refine_domain']

import os
from pathlib import Path




adata.write(save_root / 'result.h5ad')
df_PC = pd.DataFrame(data=adata.obsm['DeepST_embed'], index=adata.obs.index)
df_PC.to_csv(save_root / 'PCs.tsv', sep='\t')
adata.obs.to_csv( save_root / 'metadata.tsv', sep='\t')