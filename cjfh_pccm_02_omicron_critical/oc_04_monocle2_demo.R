rm(list=ls())#
#######
library('BiocManager')
library('devtools')
library(remotes)
library(pheatmap)
library(data.table)
library(dplyr)
library(limma)
library(monocle)
library(clusterProfiler)
library(Seurat)
library(reticulate)
library(sceasy)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)
library(reshape2)
library(ggpubr)
library(colorRamps)
###################
cds_neu <- readRDS('cds_neu.RDS')
cds_neu
cds_neu <- estimateSizeFactors(cds_neu)
cds_neu <- estimateDispersions(cds_neu)
Idents(cds_neu)='cluster'
deg_cluster <- FindAllMarkers(cds_neu)
table(deg_cluster$cluster)
express_genes <- subset(deg_cluster,p_val_adj<0.05)$gene

cds_neu <- setOrderingFilter(cds_neu,express_genes)
cds_neu
plot_ordering_genes(cds_neu)
colnames(pData(cds_neu))
#############
cds_neu <- reduceDimension(cds_neu,max_components = 2,method = 'DDRTree')
colnames(pData(cds_neu))
cds_neu <- orderCells(cds_neu)
plot_cell_trajectory(cds_neu,color_by='Pseudotime',size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_neu,color_by='cluster',size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_neu, color_by = 'State', size=1,show_backbone=TRUE)
plot_cell_trajectory(cds_neu, color_by = 'State') + facet_wrap('~State', nrow = 1)
saveRDS(cds_neu, file = 'cds_neu.RDS')
