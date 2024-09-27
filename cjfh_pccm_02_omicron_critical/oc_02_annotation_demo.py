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
umap_adata = sc.read('omicron_02_umap_adata.h5ad')
#############################
#manual_cell_type_annotation#
#############################
sc.pl.umap(umap_adata, color=['EPCAM','KRT8',#Epithelial cells
                              'CD79A','CD79B',#B cells
                              'CST3','CD74',#DC cells
                              'APOE','MRC1',#Macrophages
                              'CPA3','KIT',#Mast cell
                              'PF4','PPBP',#Megakaryocyte
                              'C1QB','LYZ',#Monocyte
                              'S100A9','CCRL2',#Neutrophil
                              'CD3D','CD3E','CD4','CD8A','CD8B'#T cells
                             ], cmap = 'Reds', vmax = 10)
umap_adata.obs['major_cluster'] = umap_adata.obs['leiden'].map(
    {
        '0': 'Neutrophil',
        '1': 'TNK',
        '2': 'Monocyte'
        ...
    }
)
###############
#marker_visual#
###############
ax = sc.pl.heatmap(umap_adata, celltype_marker, groupby='celltype', cmap='cividis',
                   dendrogram=False, vmax = 1, figsize=(2,5),
                   save = '_celltype_marker.pdf'
                  )#
###############
#density_adata#
###############
sc.tl.embedding_density(umap_adata, basis='umap', groupby='diagnosis')
#'oc_01_control', 'oc_02_survival', 'oc_03_death'
sc.pl.embedding_density(umap_adata, basis='umap',
                        key='umap_density_diagnosis',
                        group='oc_01_control',colorbar_loc=None,
                        save = 'density_01_control.pdf')
#############
#dotplot_ecm#
#############
umap_adata.layers['scaled'] = sc.pp.scale(umap_adata, copy=True).X
marker_dict = {'AEC': ['EPCAM','KRT8'],
               'BPC': ['CD79A','CD79B'],
               'DC': ['CST3','CD74'],
               'Ma': ['APOE','MRC1'],
               'MC': ['CPA3','KIT'],
               'Mega': ['PF4','PPBP'],
               'Mono': ['C1QB','LYZ'],
               'Neu': ['S100A9','CCRL2'],
               'T&NK': ['CD3D','CD3E','CD4','CD8A','CD8B'
                             ]
                }
x = sc.pl.stacked_violin(umap_adata, marker_dict, groupby='celltype',
                          cmap = 'Spectral_r', rotation=90, vmax = 1,
                          save = 'marker_adata.pdf'
                         )#
######
#Save#
######
umap_adata.write('omicron_02_umap_adata.h5ad')