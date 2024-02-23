
#from examples import *
import pyVIA.core as via
import scanpy as sc
import pandas as pd
import umap
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import numpy as np

foldername="./pyVIA_data/"

adata = sc.read_h5ad(foldername + "Anndata.h5ad")
df_ids = pd.read_csv(foldername + "Anndata_id.csv", delimiter=",")
true_label = df_ids['group_id'].tolist()
ncomps, knn, random_seed, dataset  =30,10, 2,'toy'

v0 = via.VIA(adata.obsm['X_css'][:, 0:ncomps], true_label, jac_std_global=0.15, dist_std_local=1,
             knn=knn, cluster_graph_pruning_std=1, too_big_factor=0.3, preserve_disconnected=True, dataset='group',random_seed=random_seed)
v0.run_VIA()

embedding = adata.obsm['X_umap']
via.via_streamplot(v0, embedding, scatter_size = 5)
plt.axis('equal')
plt.savefig(foldername + 'streamplot_VIA.png', dpi2 = 600)

pt = v0.single_cell_pt_markov
np.savetxt('pt.csv', pt, delimiter=',')
