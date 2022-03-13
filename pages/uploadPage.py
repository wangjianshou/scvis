import dash
import uuid
import dash_uploader as du
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, callback, State, no_update




from dash.exceptions import PreventUpdate
from os import path
from .tools import generateTable
from .cellTrans import filterData
from .setting import DATA_PATH


'''
layout = dbc.Container(
        [
            dbc.Row([dbc.Col(du.Upload(
                id='upload-file',
                text='upload file',
                text_completed='Uploaded',
                cancel_button=True,
                pause_button=True,
                filetypes=None,
                max_files=1,
                default_style={
                    'background-color': '#fafafa',
                    'font-weight': 'bold',
                    'min-height':'10px',
                    'line-height':'15px',
                },
                upload_id=uuid.uuid1()
                ), width=4),
                dbc.Col(dbc.Input(id="input-samplename", placeholder="SampleNames...", type="text",), width=2),
                dbc.Col(dbc.Input(id="input-groupname", placeholder="GroupNames...", type="text"), width=2),
                dbc.Col(dbc.Button("Submit", id='upload-data', className="me-1", color='dark', outline=True), width=1),
                dbc.Col(dbc.Alert(id="alert-upload-data", dismissable=True, is_open=False, color='success',
                                  style={"height":'60%', 'verticalAlign':'middle',},
                                  )
                        ),
                ],
                align='baseline',
                justify='baseline',
            ),
            dbc.Row(dbc.Col(id='sampleInfo', md=4)),
        ],
        style={'margin-top':'10px'},
)
'''

upload_data_tab = dbc.Tab(
    dbc.Card(
        [ dbc.CardBody(
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
    label='FilterData',
)

violinFig = dbc.Tab(dcc.Graph(id='violin'), label='Fig1')



layout = html.Div(
    dbc.Row([
                dbc.Col(dbc.Tabs([upload_data_tab, filter_data_tab], ), md=3,),
                dbc.Col(dbc.Tabs([violinFig, dbc.Tab(label='Fig2')]), md=9),
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
    if not n or sampleName is None:
        return no_update, no_update, no_update
    if 'sampleInfo' not in metadata.keys():
        metadata['sampleInfo'] = dict()
    metadata['sampleInfo'][sampleName] = (sampleName, groupName, upload_id, fileNames[0])
    table = generateTable(metadata['sampleInfo'], '100px', 10)
    if n and sampleName is not None:
        return f"Data of {sampleName} in Group {groupName} is added!", is_open and no_update or True, metadata
    elif n and sampleName is None:
        return "Please assign sample name to data", is_open and no_update or True, metadata


@callback(Output('sampleinfo-table', 'children'),
          Input('meta-data', 'data'),)
def sampleinfoTable(metadata):
    return generateTable(metadata['sampleInfo'], '100px', 10)


@callback(Output('checklist-samples', 'options'),
          Input('meta-data', 'data'),
          Input('refresh-samples-list', 'n_clicks'),
          )
def sampleChecklist(metadata, n):
    if 'sampleInfo' not in metadata.keys():
        return no_update
    return [{'label':i, 'value':i} for i in metadata['sampleInfo'].keys()]

@callback(Output('violin', 'figure'),
          Input('run-filter-data', 'n_clicks'),
          State('filter-umi', 'value'),
          State('filter-ngene', 'value'),
          State('filter-mito', 'value'),
          State('filter-ncell', 'value'),
          State('checklist-samples', 'value'),
          State("meta-data", 'data'),
          )
def filterDataPipeline(n, umi, ngene, mito, ncell, samples, metadata):
    if not n:
        return no_update
    info = metadata['sampleInfo']
    for s in samples:
        data = path.join(DATA_PATH, info[s][2], info[s][3])
        fig = filterData(umi, ngene, mito, ncell, s, data)
    return fig






