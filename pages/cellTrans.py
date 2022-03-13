import pandas as pd
import numpy as np
import scanpy as sc
import plotly.express as px
import io
from skimage.io import imread
from matplotlib import pyplot as plt

def filterData(umi, ngene, mito, ncell, sample, data):
    adata = sc.read_10x_h5(data)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=ngene)
    sc.pp.filter_genes(adata, min_cells=ncell)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    #adata.var['ribo'] = adata.var[adata.var_names.str.startswith('RP')]
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True)
    #sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #                   jitter=0.4, multi_panel=True, show=False)
    #buf = io.BytesIO()
    #plt.savefig(buf, format="png")
    #data = base64.b64encode(buf.getbuffer()).decode("utf8")
    #plt.close()
    #return "data:image/png;base64,{}".format(data)
    #buf.seek(0)
    #fig = px.imshow(imread(buf))
    x_values = list(range(1, 1001))  # 含x值的列表
    y_values = [x ** 2 for x in x_values]  # 含y值的列表
    plt.scatter(x_values, y_values, c=y_values, cmap=plt.cm.Blues, s=40)
    fig = px.bar([1, 2, 3])
    return fig