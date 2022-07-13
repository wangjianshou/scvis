import dash
import os
import numpy as np
import pandas as pd
import json
import dash_bootstrap_components as dbc
import dash_uploader as du
from flask_caching import Cache
from dash import Input, Output, State, html, dcc, callback, no_update

from os import path


from plotly import express as px
from plotly import graph_objects as go


from pages.setting import DATA_PATH, PLOTLY_LOGO
from pages import processPage, visualPage, testPage

import scanpy as sc
from copy import deepcopy

app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)
server = app.server

CACHE_CONFIG = {
    # try 'FileSystemCache' if you don't want to setup redis
    'CACHE_TYPE': 'redis',
    'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'redis://localhost:6379')
}
cache = Cache()
cache.init_app(app.server, config=CACHE_CONFIG)


sumData = dict(filterData=dict(), processData=dict(), visData=dict())


@cache.memoize()
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


@cache.memoize()
def processData(adata, n_top_genes, n_pcs, n_neighbors):
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=n_top_genes)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    #adata.layers['norm'] = sc.pp.normalize_total(adata, target_sum=1e4, copy=True).X
    #adata.layers['lognorm'] = sc.pp.log1p(adata, layer='norm', copy=True).X
    #adata.layers['scaled'] = sc.pp.scale(adata, layer='lognorm', copy=True).X
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    return adata


@cache.memoize()
def clusterData(adata, res):
    sc.tl.umap(adata)
    sc.tl.tsne(adata)
    sc.tl.leiden(adata, resolution=res)
    return adata




@cache.memoize()
def plotFilter(adata):
    v1 = go.Violin(y=adata.obs.n_genes_by_counts, xaxis='x1', yaxis='y1', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v2 = go.Violin(y=adata.obs.total_counts, xaxis='x2', yaxis='y2', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v3 = go.Violin(y=adata.obs.pct_counts_mt, xaxis='x3', yaxis='y3', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    layout = go.Layout(xaxis={'domain': [0, 0.3], 'anchor': 'x1', 'showticklabels': False, 'title':'n_genes_by_counts'},
                       xaxis2={'domain': [0.35, 0.65], 'anchor': 'x2', 'showticklabels': False, 'title':'total_counts'},
                       xaxis3={'domain': [0.7, 1], 'anchor': 'x3', 'showticklabels': False, 'title':'pct_counts_mt'},
                       yaxis2={'position': 0.34},
                       yaxis3={'position': 0.69}
                    )
    violin = go.Figure([v1, v2, v3], layout=layout)
    scatter = (px.scatter(adata.obs, x='total_counts', y='pct_counts_mt')
              .update_traces(marker={'size': 2}))
    return violin, scatter


@cache.memoize()
def plotFilterViolin(adata, x1, x2, x3):
    v1 = go.Violin(y=adata.obs[x1], xaxis='x1', yaxis='y1', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v2 = go.Violin(y=adata.obs[x2], xaxis='x2', yaxis='y2', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    v3 = go.Violin(y=adata.obs[x3], xaxis='x3', yaxis='y3', points='all', pointpos=0,
                   jitter=0.8, marker={'size': 1}, opacity=0.5, showlegend=False)
    layout = go.Layout(xaxis={'domain': [0, 0.3], 'anchor': 'x1', 'showticklabels': False, 'title':x1},
                       xaxis2={'domain': [0.35, 0.65], 'anchor': 'x2', 'showticklabels': False, 'title':x2},
                       xaxis3={'domain': [0.7, 1], 'anchor': 'x3', 'showticklabels': False, 'title':x3},
                       yaxis2={'position': 0.34},
                       yaxis3={'position': 0.69}
                    )
    violin = go.Figure([v1, v2, v3], layout=layout)
    return violin

@cache.memoize()
def plotFilterScatter(adata, x, y):
    scatter = (px.scatter(adata.obs, x=x, y=y)
              .update_traces(marker={'size': 2}))
    return scatter





du.configure_upload(app, folder=DATA_PATH,)

dropdownSinglecell = dbc.DropdownMenu(
            [
                dbc.DropdownMenuItem("Process", href="/singlecell/process",),
                dbc.DropdownMenuItem("Visual", href="/singlecell/visual",),
                dbc.DropdownMenuItem("Analysis", href="/singlecell/analysis",),
            ],
            label="Singlecell", in_navbar=True, nav=True, size='md',
)

dropdownSpatial = dbc.DropdownMenu(
            [
                dbc.DropdownMenuItem("Data", href="/spatial/data"),
                dbc.DropdownMenuItem("Process", href="/spatial/process"),
                dbc.DropdownMenuItem("Analysis", href="/spatial/analysis"),
            ],
            label="Spatial", in_navbar=True, nav=True, size='md',
)

nav = dbc.Nav(
    [
        dbc.NavLink("Login",  href="/login", style={'fontSize':18}),
        #dbc.NavLink("MyData", href="/mydata", style={'fontSize':18}),
        dbc.NavItem(dbc.Button("seeMyData", id='sample-info-button', n_clicks=0, size='md',
                               className='me-1', color='info', outline=False)),
        #dbc.NavLink("Another", href="/another", style={'fontSize':18}),
    ]
)

navbar = dbc.Navbar(
        dbc.Container(
            [
                dbc.Row(
                    [
                         dbc.Col(html.A(
                                 html.Img(src=PLOTLY_LOGO, height='30px'),
                                 href='http://www.biomarker.com.cn',
                                 style={"textDecoration": "none"},
                         )),
                        #dbc.Col(dcc.Link(
                        #    dbc.NavbarBrand("Home", className="ms-2"),
                        #    href='/home', style={"textDecoration": "none"},
                        #), ),
                        dbc.Col(dcc.Link(
                            dbc.NavItem("Home", className="ms-2"),
                            href='/', style={"textDecoration": "none"},
                        ), ),
                        #dbc.Col(dcc.Link(
                        #    dbc.NavbarBrand("Singlecell", className="ms-2"),
                        #    href='/singlecell', style={"textDecoration": "none"},
                        #    ),),
                        dbc.Col(
                            dbc.NavItem(dropdownSinglecell, className="ms-2"),
                        ),
                        dbc.Col(
                            dbc.NavItem(dropdownSpatial, className="ms-2"),
                        ),
                        dbc.Col(dcc.Link(
                            dbc.NavItem("Another", className="ms-2"),
                            href='/another', style={"textDecoration": "none"},
                        ),)
                    ],
                align='center',
                className='g-0'),
                nav,
            ],
            fluid=False,
        ),
    color='dark',
    dark=True,
)

sampleCanvas = dbc.Offcanvas(
                         html.Div(id='sampleinfo-table'),
                         id="offcanvas-sampleinfo",
                         title="sample information",
                         is_open=False,
                         placement="top")



app.layout = html.Div([dcc.Location(id="url",  refresh=False),
                       dcc.Store('meta-data', storage_type='session', data=dict()),
                       dcc.Store('meta-data2', storage_type='session', data=dict()),
                       sampleCanvas,
                       navbar,
                       dbc.Container(id='page-temp', className='pt-4'),
                       dbc.Container(id="page-content", className="pt-4"),
                       ])



@app.callback(Output("page-content", "children"),
              Input("url", "pathname"))
def render_page_content(pathname):
    if pathname == "/":
        return dash.no_update
    elif pathname == "/singlecell/process":
        return processPage.layout
    elif pathname == "/singlecell/visual":
        return visualPage.layout
    elif pathname == "/another":
        return testPage.layout

@app.callback(
    Output("offcanvas-sampleinfo", "is_open"),
    Input("sample-info-button", "n_clicks"),
    State("offcanvas-sampleinfo", "is_open"),
)
def toggle_offcanvas(n1, is_open):
    if n1:
        return not is_open
    return is_open



@callback(output=[
              Output('violin', 'figure'),
              Output('scatter', 'figure'),
              Output('filter-sample-fig', 'value'),
              Output('filter-sample-fig', 'options'),
              Output('meta-data2', 'data')
          ],
          inputs = dict(
              n=Input('run-filter-data', 'n_clicks'),
              metadata=Input("meta-data", 'data'),
              v1=Input('filter-violin-1', 'value'),
              v2=Input('filter-violin-2', 'value'),
              v3=Input('filter-violin-3', 'value'),
              scatter_x=Input('filter-scatter-x', 'value'),
              scatter_y=Input('filter-scatter-y', 'value'),
              sample=Input('filter-sample-fig', 'value'),

          ),
          state = dict(
              umi=State('filter-umi', 'value'),
              ngene=State('filter-ngene', 'value'),
              mito=State('filter-mito', 'value'),
              ncell=State('filter-ncell', 'value'),
              samples=State('checklist-samples', 'value'),
              metadata2=State('meta-data2', 'data'),
          ),
          prevent_initial_call=True,
          )
def filterDataPipeline(n, metadata, v1, v2, v3, scatter_x, scatter_y, sample,
                       umi, ngene, mito, ncell, samples, metadata2):
    r_violin = r_scatter = r_sample = r_options = r_metadata2 = no_update
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = 'No clicks yet'
        return no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if 'sampleInfo' not in metadata.keys():
        return no_update
    if samples is None:
        return no_update
    info = deepcopy(metadata2.get('sampleInfo', metadata['sampleInfo']))
    if button_id == 'meta-data':
        metadata2['sampleInfo'] = metadata2.get('sampleInfo', dict())
        metadata2['sampleInfo'].update(metadata['sampleInfo'])
        r_metadata2 = deepcopy(metadata2)
    elif button_id == 'run-filter-data':
        data = {s: path.join(DATA_PATH, info[s]['uuid'], info[s]['fileName']) for s in samples}
        scData = {s: filterData(umi, ngene, mito, ncell, s, data[s]) for s in samples}
        sumData['filterData'].update(scData)
        for s in samples:
            #info[s]['filterParameter'] = f"umi={umi};ngene={ngene};mito={mito};ncell={ncell}"
            info[s]['filterParameter'] = json.dumps(dict(umi=umi,ngene=ngene, mito=mito, ncell=ncell))
        metadata2['sampleInfo'] = info
        r_violin, r_scatter = plotFilter(scData[samples[0]])
        r_options = [{'label': i, 'value': i} for i in info.keys() if info[s].get('filterParameter', 0)]
        r_metadata2 = deepcopy(metadata2)
        r_sample = samples[0]
    else:
        p = json.loads(info[sample]['filterParameter'])
        dpath = path.join(DATA_PATH, info[sample]['uuid'], info[sample]['fileName'])
        adata = filterData(p['umi'], p['ngene'], p['mito'], p['ncell'], sample, dpath)
    if button_id=='filter-sample-fig':
        r_violin = plotFilterViolin(adata, v1, v2, v3)
        r_scatter = plotFilterScatter(adata, scatter_x, scatter_y)
    elif button_id in ['filter-violin-1', 'filter-violin-2', 'filter-violin-3']:
        r_violin = plotFilterViolin(adata, v1, v2, v3)
    elif button_id in ['filter-scatter-x', 'filter-scatter-y']:
        r_scatter = plotFilterScatter(adata, scatter_x, scatter_y)
    return r_violin, r_scatter, r_sample, r_options, r_metadata2




@callback(output=[
        Output('pca-main', 'figure'),
        Output('pca-color-attr', 'options'),
        Output('pca-color-attr', 'value'),
    ],
    inputs=dict(
        n=Input('run-process-data', 'n_clicks'),
        gene=Input('pca-color-attr', 'value'),
        sample=State('filter-sample-fig', 'value'),
    ),
    state=dict(
        samples=State('process-checklist-samples', 'value'),
        nTopGene=State('n_top_genes', 'value'),
        nPCA=State('n_pcs', 'value'),
        nNeighbor=State('n_neighbors', 'value'),
        #metadata2=State('meta-data2', 'data')
    ),
    prevent_initial_call=True,
)
def processDataPipeline(n, gene, sample, samples, nTopGene, nPCA, nNeighbor): #, metadata2):
    r_fig = r_val = r_opt = no_update
    if samples is None:
        return no_update
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = 'No clicks yet'
        return no_update
    else:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    #info = deepcopy(metadata2['sampleInfo'])
    #tmpData = dict()
    #for s in samples:
    #    p = json.loads(info[s]['filterParameter'])
    #    dpath = path.join(DATA_PATH, info[s]['uuid'], info[s]['fileName'])
    #    tmpData[s] = filterData(p['umi'], p['ngene'], p['mito'], p['ncell'], s, dpath)
    #scData = {s:processData(tmpData[s], nTopGene, nPCA, nNeighbor) for s in samples}
    if 'run-process-data' == button_id:
        scData = {s: processData(sumData['filterData'][s], nTopGene, nPCA, nNeighbor) for s in samples}
        sumData['processData'].update(scData)
    adata = sumData['processData'][sample]
    df = pd.DataFrame(adata.obsm['X_pca'][:, :3], columns=['pca_' + str(i) for i in range(1, 4)])
    if button_id != 'pca-color-attr':
        r_opt = [{'label': i, 'value': i} for i in adata.var_names]
        r_val = adata.var_names[0]
        gene = r_val
    r_fig = (px.scatter(df, 'pca_1', 'pca_2',
                        color=sc.get.obs_df(adata, [gene, ]).values.squeeze(),
                        color_continuous_scale=['white', 'blue'])
             .update_traces(marker={'size': 3}))
    return r_fig, r_opt, r_val






#可视化回调（不要删除！不要删除！！不要删除！！！）
'''
@callback(output=[
            Output('visual-tsne', 'figure'),
            Output('visual-tsne-gene', 'figure'),
            Output('visual-umap', 'figure'),
            Output('visual-umap-gene', 'figure'),
            Output('visual-sample', 'options'),
            Output('visual-sample', 'value'),
            Output('tsne-gene-single', 'options'),
            Output('tsne-gene-single', 'value'),
            Output('tsne-gene-multi', 'options'),
            Output('tsne-gene-multi', 'value'),
            Output('umap-gene-single', 'options'),
            Output('umap-gene-single', 'value'),
            Output('umap-gene-multi', 'options'),
            Output('umap-gene-multi', 'value'),
          ],
          inputs=dict(
              n=Input('run-visual-data', 'n_clicks'),
              sample=Input('visual-sample', 'value'),
              tsne_gene_single=Input('tsne-gene-single', 'value'),
              tsne_gene_multi=Input('tsne-gene-multi', 'value'),
              umap_gene_single=Input('umap-gene-single', 'value'),
              umap_gene_multi=Input('umap-gene-multi', 'value'),
          ),
          state=dict(
              samples=State('checklist-visual-samples', 'value'),
              res=State('resolution', 'value'),
              sample_opt = State('visual-sample', 'options'),
          )
)
def visualPipeline(n, sample, tsne_gene_single, tsne_gene_multi, umap_gene_single, umap_gene_multi,
                   samples, res, sample_opt):
    r_tsne = r_tsne_gene = r_umap = r_umap_gene = no_update
    r_val = r_opt = r_gene = r_gene_val = r_gene_leiden = r_gene_single_opt = no_update
    if samples is None:
        return no_update
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        button_id = [i['prop_id'].split('.')[0] for i in ctx.triggered]
    if 'run-visual-data' in button_id:
        scData = {s: clusterData(sumData['processData'][s], res) for s in samples}
        sumData['visData'].update(scData)
        r_val = samples[0]
        if sample_opt is None:
            r_opt = [{'label':i,'value':i} for i in samples]
        else:
            tmp = [s['value'] for s in sample_opt]
            r_opt = [{'label':i,'value':i} for i in samples if i not in tmp] + sample_opt
        sample = r_val
    adata = sumData['visData'][sample]
    if not (set(button_id) & {'tsne-gene-single', 'tsne-gene-multi',
                              'umap-gene-single', 'umap-gene-multi'}):
        r_gene = [{'label':i, 'value':i} for i in adata.var.index]
        r_gene_leiden = 'leiden'
        umap_gene_single = tsne_gene_single = r_gene_leiden
        umap_gene_multi = tsne_gene_multi = adata.var_names[:3].tolist()
        r_gene_val = tsne_gene_multi
        r_gene_single_opt = [{'label':'leiden', 'value':'leiden'}] + r_gene
    df = pd.DataFrame(adata.obsm['X_tsne'], columns=['tsne-1', 'tsne-2'], index=adata.obs.index)
    dff = pd.concat([df, sc.get.obs_df(adata, tsne_gene_multi)], axis=1)
    dff = pd.melt(dff, value_vars=tsne_gene_multi, id_vars=['tsne-1', 'tsne-2'],
                  var_name='Gene', value_name='exp')
    r_tsne_gene = px.scatter(dff, x='tsne-1', y='tsne-2', color='exp',
                     facet_col='Gene', facet_col_wrap=3,
                     range_color=list(dff.exp.quantile([0.01, 0.99])),
                     color_continuous_scale='Blues'
                     )
    r_tsne_gene.update_traces(marker={'size':3})
    #for i in range(len(r_tsne_gene.data)):
    #    r_tsne_gene.data[i].update(marker={'size': 2})
    df = pd.concat([df, sc.get.obs_df(adata, [tsne_gene_single])], axis=1)
    df['leiden'] = adata.obs['leiden'].values
    r_tsne = (px.scatter(df, x='tsne-1', y='tsne-2', color=tsne_gene_single,
                         color_continuous_scale='Blues')
             .update_traces(marker={'size':3}))

    umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['umap-1', 'umap-2'], index=adata.obs.index)
    umap_dff = pd.concat([umap_df, sc.get.obs_df(adata, umap_gene_multi)], axis=1)
    umap_dff = pd.melt(umap_dff, value_vars=umap_gene_multi, id_vars=['umap-1', 'umap-2'],
                       var_name='Gene', value_name='exp')
    r_umap_gene = px.scatter(umap_dff, x='umap-1', y='umap-2', color='exp',
                             facet_col='Gene', facet_col_wrap=3,
                             range_color=list(umap_dff.exp.quantile([0.01, 0.99])),
                             color_continuous_scale='Blues'
                             )
    r_umap_gene.update_traces(marker={'size': 3})

    umap_df = pd.concat([umap_df, sc.get.obs_df(adata, [umap_gene_single])], axis=1)
    umap_df['leiden'] = adata.obs['leiden'].values
    r_umap = (px.scatter(umap_df, x='umap-1', y='umap-2', color=umap_gene_single,
                         color_continuous_scale='Blues')
              .update_traces(marker={'size': 3}))
    return [r_tsne, r_tsne_gene, r_umap, r_umap_gene, r_opt, r_val, r_gene_single_opt,
            r_gene_leiden, r_gene, r_gene_val, r_gene_single_opt, r_gene_leiden,
            r_gene, r_gene_val]
'''



@callback(output=[
            Output('tsne-gene-single', 'options'),
            Output('tsne-gene-single', 'value'),
            Output('tsne-gene-multi', 'options'),
            Output('tsne-gene-multi', 'value'),
            Output('umap-gene-single', 'options'),
            Output('umap-gene-single', 'value'),
            Output('umap-gene-multi', 'options'),
            Output('umap-gene-multi', 'value'),
            Output('visual-sample', 'options'),
            Output('visual-sample', 'value'),
          ],
          inputs=dict(
              n=Input('run-visual-data', 'n_clicks'),
              sample=Input('visual-sample', 'value'),
          ),
          state=dict(
              samples=State('checklist-visual-samples', 'value'),
              res=State('resolution', 'value'),
              sample_opt = State('visual-sample', 'options'),
          )
)
def visualOpt(n, sample, samples, res, sample_opt):
    r_opt = r_val = r_gene = r_gene_single_opt = r_gene_val = r_single_val = no_update
    if samples is None:
        return no_update
    scData = {s: clusterData(sumData['processData'][s], res) for s in samples}
    sumData['visData'].update(scData)
    if sample_opt is None:
        r_opt = [{'label':i,'value':i} for i in samples]
        r_val = samples[0]
        sample = r_val
    else:
        tmp = [s['value'] for s in sample_opt]
        r_opt = [{'label':i,'value':i} for i in samples if i not in tmp] + sample_opt
    adata = sumData['visData'][sample]
    r_gene = [{'label':i, 'value':i} for i in adata.var.index]
    r_gene_single_opt = [{'label':'leiden', 'value':'leiden'}] + r_gene
    r_gene_val = adata.var_names[:3].tolist()
    r_gene_single_val = 'leiden'
    return [r_gene_single_opt, r_gene_single_val, r_gene, r_gene_val, r_gene_single_opt,
            r_gene_single_val, r_gene, r_gene_val, r_opt, r_val]






@callback(output=[
            Output('visual-tsne', 'figure'),
            Output('visual-tsne-gene', 'figure'),
            Output('visual-umap', 'figure'),
            Output('visual-umap-gene', 'figure'),
          ],
          inputs=dict(
              tsne_gene_single=Input('tsne-gene-single', 'value'),
              tsne_gene_multi=Input('tsne-gene-multi', 'value'),
              umap_gene_single=Input('umap-gene-single', 'value'),
              umap_gene_multi=Input('umap-gene-multi', 'value'),
          ),
          state=dict(
              sample=State('visual-sample', 'value'),
          )
)
def visualPipeline(tsne_gene_single, tsne_gene_multi, umap_gene_single, umap_gene_multi,
                   sample):
    r_tsne = r_tsne_gene = r_umap = r_umap_gene = no_update
    adata = sumData['visData'][sample]

    df = pd.DataFrame(adata.obsm['X_tsne'], columns=['tsne-1', 'tsne-2'], index=adata.obs.index)
    dff = pd.concat([df, sc.get.obs_df(adata, tsne_gene_multi)], axis=1)
    dff = pd.melt(dff, value_vars=tsne_gene_multi, id_vars=['tsne-1', 'tsne-2'],
                  var_name='Gene', value_name='exp')
    r_tsne_gene = px.scatter(dff, x='tsne-1', y='tsne-2', color='exp',
                     facet_col='Gene', facet_col_wrap=3,
                     range_color=list(dff.exp.quantile([0.01, 0.99])),
                     color_continuous_scale='Blues'
                     )
    r_tsne_gene.update_traces(marker={'size':3})
    df = pd.concat([df, sc.get.obs_df(adata, [tsne_gene_single])], axis=1)
    df['leiden'] = adata.obs['leiden'].values
    r_tsne = (px.scatter(df, x='tsne-1', y='tsne-2', color=tsne_gene_single,
                         color_continuous_scale='Blues')
             .update_traces(marker={'size':3}))
    umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['umap-1', 'umap-2'], index=adata.obs.index)
    umap_dff = pd.concat([umap_df, sc.get.obs_df(adata, umap_gene_multi)], axis=1)
    umap_dff = pd.melt(umap_dff, value_vars=umap_gene_multi, id_vars=['umap-1', 'umap-2'],
                       var_name='Gene', value_name='exp')
    r_umap_gene = px.scatter(umap_dff, x='umap-1', y='umap-2', color='exp',
                             facet_col='Gene', facet_col_wrap=3,
                             range_color=list(umap_dff.exp.quantile([0.01, 0.99])),
                             color_continuous_scale='Blues'
                             )
    r_umap_gene.update_traces(marker={'size': 3})

    umap_df = pd.concat([umap_df, sc.get.obs_df(adata, [umap_gene_single])], axis=1)
    umap_df['leiden'] = adata.obs['leiden'].values
    r_umap = (px.scatter(umap_df, x='umap-1', y='umap-2', color=umap_gene_single,
                         color_continuous_scale='Blues')
              .update_traces(marker={'size': 3}))
    return [r_tsne, r_tsne_gene, r_umap, r_umap_gene, ]










if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8050, debug=True)