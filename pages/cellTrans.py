import pandas as pd
import numpy as np
import scanpy as sc
import plotly.express as px
import io
from skimage.io import imread
from matplotlib import pyplot as plt
#import matplotlib
#matplotlib.use("agg")
from plotly import graph_objects as go


'''
def processData(adata, n_top_genes, n_pcs, n_neighbors):
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata
'''






# 请不要删除，请不要删除，请不要删除！
'''
def filterData(umi, ngene, mito, ncell, sample, data):
    adata = sc.read_10x_h5(data)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=ngene)
    sc.pp.filter_genes(adata, min_cells=ncell)
    sc.pp.filter_cells(adata, min_counts=umi)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #adata.var['ribo'] = adata.var[adata.var_names.str.startswith('RP')]
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
    return adata

'''




'''
def plotFilter(adata):
    #fig = sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #                   jitter=0.4, multi_panel=True, show=False)
    #buf = io.BytesIO()
    #plt.savefig(buf, format="png", dpi=300)
    #fig = px.imshow(imread(buf))
    #fig = fig.update_layout(xaxis={'showticklabels': False}, yaxis={'showticklabels': False})
    v1 = go.Violin(y=adata.obs.n_genes_by_counts, xaxis='x1', yaxis='y1', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v2 = go.Violin(y=adata.obs.total_counts, xaxis='x2', yaxis='y2', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v3 = go.Violin(y=adata.obs.pct_counts_mt, xaxis='x3', yaxis='y3', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    layout = go.Layout(xaxis={'domain': [0, 0.3], 'anchor': 'x1', 'showticklabels': False, 'title':'aaaa'},
                       xaxis2={'domain': [0.35, 0.65], 'anchor': 'x2', 'showticklabels': False, 'title':'aaaa'},
                       xaxis3={'domain': [0.7, 1], 'anchor': 'x3', 'showticklabels': False, 'title':'aaaa'},
                       yaxis2={'position': 0.34},
                       yaxis3={'position': 0.69}
                    )
    violin = go.Figure([v1, v2, v3], layout=layout)
    scatter = (px.scatter(adata.obs, x='total_counts', y='pct_counts_mt')
              .update_traces(marker={'size': 2}))

    return violin, scatter
'''