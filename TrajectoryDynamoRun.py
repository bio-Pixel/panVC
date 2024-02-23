import warnings
warnings.filterwarnings('ignore')
import dynamo as dyn
import scvelo as scv
mport scanpy as sc
import os

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)


############################### Dynamo ###############################
adata=sc.read_h5ad("path_to_h5ad")
preprocessor = dyn.pp.Preprocessor(cell_cycle_score_enable=True)
preprocessor.preprocess_adata(adata)

dyn.tl.dynamics(adata, model='stochastic', cores=3)
dyn.tl.reduceDimension(adata)
dyn.tl.cell_velocities(adata, method='pearson', other_kernels_dict={'transform': 'sqrt'}) 
dyn.tl.cell_wise_confidence(adata)
dyn.tl.confident_cell_velocities(adata, group='celltype')
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltype'], palette=["#FFBE58", "#CEB7A8", "#8C7490"], save='Dynamo.svg', title='')
