import numpy as np
import pandas as pd
import h5py
import os
import argparse
import scipy.sparse as sparse
#import pygwalker as pyg
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm, rcParams
import csv
import scanpy as sc
import scanpy.external as sce
import anndata
import bbknn
#import scanorama
#import scrublet as scr
import scvelo as scv
#import cellrank as cr
import loompy as lp
from scipy.stats import pearsonr
#import palantir
#import doubletdetection
from igraph import *
#from MulticoreTSNE import MulticoreTSNE as TSNE #faster TSNE alternative
from anndata import read_h5ad
from anndata import read_csv
from matplotlib import rcParams
sc.settings.verbosity = 3
sc.logging.print_versions()

%matplotlib inline

#read_adata
adata = sc.read('omicron_02_umap_adata.h5ad')

scv.pl.proportions(adata)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata, n_jobs = 64)#

scv.pl.velocity_embedding_stream(adata, basis='umap',save = 'velo_01_neu_01_cluster.pdf')

scv.pl.velocity(adata, ['ISG15',  'IFITM3', 'GBP1', 'IL1B'], ncols=2,save = 'velo_01_neu_02_marker.pdf')

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

df = adata.obs.groupby('clusters')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)

