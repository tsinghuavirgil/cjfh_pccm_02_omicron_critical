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

sc.settings.set_figure_params(300)#
#read_adata
adata = sc.read('adata_01_raw.h5ad')

sc.pl.highest_expr_genes(adata, n_top=20, )

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

adata = adata[adata.obs.pct_counts_mt < 10, :].copy()

sc.pp.normalize_total(umap_adata, target_sum=1e4)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.pl.highly_variable_genes(adata)

adata.raw = adata

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, use_highly_variable=True)

sc.pl.pca_scatter(adata, color=['tissue'])

sc.pl.pca_variance_ratio(adata, log = True)

sc.pl.pca_loadings(adata)

neighbor = 30
resolution = 1
pcs = 15

sce.pp.harmony_integrate(adata, 'sample_id')

sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors = neighbor, n_pcs = pcs)

sc.tl.leiden(adata, resolution = 1)

sc.tl.umap(adata)

sc.pl.umap(adata, color=['leiden'])

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

deg_adata = pd.DataFrame(umap_adata.uns['rank_genes_groups']['names']).head(50)
deg_adata

adata.write('omicron_02_umap_adata.h5ad')