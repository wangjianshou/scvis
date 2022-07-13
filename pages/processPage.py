import dash
import uuid
import dash_uploader as du
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, callback, State, no_update

from os import path
from copy import deepcopy
from .tools import generateTable
from .setting import DATA_PATH


upload_data_tab = dbc.Tab(
    dbc.Card(
        [dbc.CardBody(
            [
            du.Upload(
                id='upload-file', text='upload file', text_completed='Uploaded',
                cancel_button=True, pause_button=True, filetypes=None, max_files=1,
                default_style={
                    'background-color': '#fafafa',
                    'font-weight': 'bold',
                    'min-height':'10px',
                    'line-height':'15px',
                },
                upload_id=uuid.uuid1(),
            ),

            dbc.Input(id="input-samplename", placeholder="SampleNames...", type="text", style={'margin-top':'20px'}),
            dbc.Input(id="input-groupname", placeholder="GroupNames...", type="text", style={'margin-top':'10px'}),
            html.Br(),
            dbc.Row(dbc.Col(dbc.Button("Add", id='upload-data', className="me-1", color='dark', outline=True),
                            width={'offset':8, 'width':4})),
            dbc.Alert(id="alert-upload-data", dismissable=True, is_open=False, color='success',
                                      style={"height":'60%', 'verticalAlign':'middle',},),
            #dbc.Row(dbc.Col(id='sampleInfo')),
            ],)
        ],
        #style={"width": "18rem"},
    ),
    label='Upload',
)

filter_data_tab = dbc.Tab(
    dbc.Card(
        [dbc.CardBody(
            [
                html.H6('Filter Cells Parameters:'),
                dbc.InputGroup([dbc.InputGroupText("nUMI", style={'width':'40%'}), dbc.Input(value=500, id='filter-umi')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                dbc.InputGroup([dbc.InputGroupText("nGene", style={'width':'40%'}), dbc.Input(value=200, id='filter-ngene')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                dbc.InputGroup([dbc.InputGroupText("mitoRatio", style={'width':'40%'}), dbc.Input(value=0.2, id='filter-mito')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                html.Br(),
                html.H6("Filter Genes Parameter:"),
                dbc.InputGroup([dbc.InputGroupText("nCell", style={'width':'40%'}), dbc.Input(value=10, id='filter-ncell')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                html.Br(),
                html.H6("Apply These Parameters to:"),
                dbc.Button(html.H6('Apply These Parameters to:'), id='refresh-samples-list', size='sm', color='white',
                           outline=False, style={'padding':'0px', 'cursor':'default', "visibility": "hidden"}),
                dbc.Checklist(id="checklist-samples", persistence_type='session', inline=True, style={"textAlign":'center'}),
                html.Br(),
                dbc.Row(dbc.Col(dbc.Button("Run", id='run-filter-data', className="me-1", color='dark', outline=True),
                                width={'offset': 8, 'width': 4})),
            ],
        ),],
    #style={"width": "18rem"},
    ),
    label='Filter',
)


process_data_tab = dbc.Tab(
    dbc.Card(
        [dbc.CardBody(
            [
                html.H6('Variable Features Parameters:'),
                dbc.InputGroup([dbc.InputGroupText("nTopGene", style={'width':'40%'}), dbc.Input(value=2000, id='n_top_genes')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                dbc.InputGroup([dbc.InputGroupText("nPCA", style={'width':'40%'}), dbc.Input(value=50, id='n_pcs')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                dbc.InputGroup([dbc.InputGroupText("nNeighbor", style={'width':'40%'}), dbc.Input(value=10, id='n_neighbors')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                html.Br(),
                html.H6("Apply These Parameters to:"),
                dbc.Button(html.H6('Apply These Parameters to:'), id='refresh-process-samples-list', size='sm', color='white',
                           outline=False, style={'padding':'0px', 'cursor':'default', "visibility": "hidden"}),
                dbc.Checklist(id="process-checklist-samples", persistence_type='session', inline=True, style={"textAlign":'center'}),
                html.Br(),
                dbc.Row(dbc.Col(dbc.Button("Run", id='run-process-data', className="me-1", color='dark', outline=True),
                                width={'offset': 8, 'width': 4})),
            ],
        ),],
    #style={"width": "18rem"},
    ),
    label='Process',
)



tmp_x = [{'label':'n_genes_by_counts', 'value':'n_genes_by_counts'},
         {'label':'total_counts', 'value':'total_counts'},
         {'label':'pct_counts_mt', 'value':'pct_counts_mt'},]

filterDataFig = dbc.Tab(
    [
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='filter-violin-1', value='n_genes_by_counts', options=tmp_x, clearable=False), md=3),
                dbc.Col(dcc.Dropdown(id='filter-violin-2', value='total_counts', options=tmp_x, clearable=False), md=3),
                dbc.Col(dcc.Dropdown(id='filter-violin-3', value='pct_counts_mt', options=tmp_x, clearable=False), md=3),
            ],
            justify='evenly', style={'margin-top':'20px'},
        ),
        dcc.Graph(id='violin'),
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='filter-scatter-x', value='total_counts', options=tmp_x, clearable=False), md=3),
                dbc.Col(dcc.Dropdown(id='filter-scatter-y', value='pct_counts_mt', options=tmp_x, clearable=False), md=3),
            ],
            justify='evenly', style={'margin-top': '20px'},
        ),
        dcc.Graph(id='scatter'),
    ],
    label='FilterFigure'
)


processDataFig = dbc.Tab(
    [
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='pca-color-attr', clearable=False, searchable=True), md=3),
                #dbc.Col(dcc.Dropdown(id='pca-', value='', options=[], clearable=False), md=3),
            ],
            justify='evenly', style={'margin-top': '20px'},
        ),
        dcc.Graph(id='pca-main'),
        #dcc.Graph(id='pca-scree'),
    ],
    label='ProcessFigure'
)




figParameter = dbc.Card(
    [dbc.CardBody(
        [
            html.H6('Visualization Parameters:'),
            dcc.Dropdown(id='filter-sample-fig', placeholder="Select sample"),
        ],
    ),],
    style={"margin-top": "20px"},
)



layout = html.Div(
    dbc.Row(
        [
            dbc.Col([dbc.Tabs([upload_data_tab, filter_data_tab, process_data_tab], ), figParameter], md=3,),
            dbc.Col(dbc.Tabs([filterDataFig, processDataFig, ]), md=9),
        ]
    )
)



@callback(Output('alert-upload-data', 'children'),
          #Output('sampleInfo', 'children'),
          Output('alert-upload-data', 'is_open'),
          Output('meta-data', 'data'),
          Input('upload-data', 'n_clicks'),
          State('input-samplename', 'value'),
          State('input-groupname', 'value'),
          State('upload-file', 'upload_id'),
          State('upload-file', 'fileNames'),
          State('alert-upload-data', 'is_open'),
          State('meta-data', 'data')
        )
def uploadData(n, sampleName, groupName, upload_id, fileNames, is_open, metadata):
    metadata = deepcopy(metadata)
    if not n or sampleName is None:
        return no_update, no_update, no_update
    if 'sampleInfo' not in metadata.keys():
        metadata['sampleInfo'] = dict()
    metadata['sampleInfo'][sampleName] = dict(sampleName=sampleName, groupName=groupName, uuid=upload_id, fileName=fileNames[0])
    #table = generateTable(metadata['sampleInfo'], '100px', 10)
    if n and sampleName is not None:
        return f"Data of {sampleName} in Group {groupName} is added!", is_open and no_update or True, metadata
    elif n and sampleName is None:
        return "Please assign sample name to data", is_open and no_update or True, metadata



@callback(Output('sampleinfo-table', 'children'),
          Input('meta-data', 'data'),
          Input('meta-data2', 'data'),
          )
def sampleinfoTable(metadata, metadata2):
    if 'sampleInfo' not in metadata.keys():
        return no_update
    elif metadata2 == dict():
        return generateTable(metadata['sampleInfo'], '200px', 10)
    else:
        return generateTable(metadata2['sampleInfo'], '200px', 10)



@callback(Output('checklist-samples', 'options'),
          Output('process-checklist-samples', 'options'),
          Input('refresh-samples-list', 'n_clicks'),
          Input("refresh-process-samples-list", 'n_clicks'),
          Input('meta-data', 'data'),
          Input('meta-data2', 'data')
          )
def sampleChecklist(n_filter, n_process, metadata, metadata2):
    r_filter = r_process = no_update
    ctx = dash.callback_context
    if not ctx.triggered:
        button_id = None
    else:
        #button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        button_id = [i['prop_id'].split('.')[0] for i in ctx.triggered]
    if button_id is None:
        r_filter = ('sampleInfo' in metadata.keys() and
                    [{'label': i, 'value': i} for i in metadata['sampleInfo'].keys()]
                    or no_update)
        r_process = ('sampleInfo' in metadata2.keys() and
                     [{'label': i, 'value': i}
                        for i in metadata2['sampleInfo'].keys()
                        if "filterParameter" in metadata2['sampleInfo'][i].keys()]
                     or no_update)
    elif set(button_id) & {'meta-data', 'refresh-samples-list'}:
        r_filter =  [{'label':i, 'value':i} for i in metadata['sampleInfo'].keys()]
    elif set(button_id) & {'meta-data2', "refresh-process-samples-list"}:
        r_process = [{'label':i, 'value':i}
                     for i in metadata2['sampleInfo'].keys()
                     if "filterParameter" in metadata2['sampleInfo'][i].keys()]
    return r_filter, r_process


'''
@callback(Output('checklist-samples', 'options'),
          Input('refresh-samples-list', 'n_clicks'),
          Input('meta-data', 'data'),
          )
def sampleChecklist(n_filter, metadata):
    if 'sampleInfo' not in metadata.keys():
        return no_update
    r_filter =  [{'label':i, 'value':i} for i in metadata['sampleInfo'].keys()]
    return r_filter


@callback(Output('process-checklist-samples', 'options'),
          Input("refresh-process-samples-list", 'n_clicks'),
          Input('meta-data2', 'data')
          )
def sampleChecklist(n, metadata2):
    if 'sampleInfo' not in metadata2.keys():
        return no_update
    r_process = [{'label':i, 'value':i}
                   for i in metadata2['sampleInfo'].keys()
                   if "filterParameter" in metadata2['sampleInfo'][i].keys()]
    return r_process
'''

