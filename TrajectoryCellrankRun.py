import sys
import numpy as np
import pandas as pd
import cellrank as cr
import scanpy as sc
import scvelo as scv
import os
scv.settings.verbosity = 3
scv.settings.set_figure_params("scvelo")
cr.settings.verbosity = 2
import warnings
warnings.simplefilter("ignore", category=UserWarning)

def get_files_by_suffix(path, suffix):
	files = []
	for file_name in os.listdir(path):
		if file_name.endswith(suffix):
			files.append(os.path.join(path, file_name))
	return files


def del_files_by_suffix(path, suffix):
	files = []
	for file_name in os.listdir(path):
		if not file_name.endswith(suffix):
			files.append(os.path.join(path, file_name))
	return files

############################### Cellrank ###############################
adata = sc.read_h5ad("path_to_h5ad")
#Preprocess the data
scv.pl.proportions(adata)
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

scv.tl.recover_dynamics(adata, n_jobs=40)
scv.tl.velocity(adata, mode="dynamical")
vk = cr.kernels.VelocityKernel(adata)
vk.compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
print(combined_kernel)
vk.plot_projection(basis='umap',color='celltype',palette=["#FFBE58", "#CEB7A8", "#8C7490"],save="Cellrank.svg",recompute=True,stream=True)
