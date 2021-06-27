# SEDR
## data folder structure
SEDR/  
----data/  
--------DLPFC/  
------------151507/  
----output/  
--------DLPFC/  
------------151507/  
----script/    
------------DLPFC_comp.R   

Please generate folder with structure shown above.  
Download the scripts here and put them in script folder.
DLPFC data can be downloaded from [SpatialLIBD](http://spatial.libd.org/spatialLIBD/). Extract and put all data within data/DLPFC folder.   

Table of DLPFC data:

|Sample_ID|n_cluster|
| ------------- |:-------------:|
|151507|7|
|151508|7|
|151509|7|
|151510|7|
|151669|5|
|151670|5|
|151671|5|
|151672|5|
|151673|7|
|151674|7|
|151675|7|
|151676|7|

Then run following codes to generate results of SEDR and other methods. 

## Run state-of-the-art methods
* Rscript DLPFC_Seurat.R sample
* python DLPFC_stLearn.py sample
* python DLPFC_SpaGCN.py sample n_clusters
* Rscript DLPFC_BayesSpace.R sample n_clusters
* Rscript DLPFC_Giotto.R sample n_clusters

## Do clustering for SEDR
* Follow the instructions in [SEDR](https://github.com/HzFu/SEDR) to run SEDR. 
* python DLPFC_SEDR_clustering.py sample n_clusters

## Do trajectory analyses for SEDR and Seurat
* Rscript DLPFC_trajectory.R sample

## Do PAGA analyses
* python DLPFC_PAGA.py  or run DLPFC_calculate_PAGA.ipynb
* DLPFC_PAGA.weight_comp.R

## Compare SEDR and other methods
* Rscript DLPFC_comp.R sample
* Rscript DLPFC.ARI_boxplot.R

