import dash
import uuid
import dash_uploader as du
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, callback, State, no_update



cluster_tab = dbc.Tab(
    dbc.Card(
        [dbc.CardBody(
            [
                html.H6('Cluster Parameters:'),
                dbc.InputGroup([dbc.InputGroupText("resolution", style={'width':'40%'}), dbc.Input(value=0.8, id='resolution')],
                               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                #dbc.InputGroup([dbc.InputGroupText("nGene", style={'width':'40%'}), dbc.Input(value=200, id='filter-ngene')],
                #               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                #dbc.InputGroup([dbc.InputGroupText("mitoRatio", style={'width':'40%'}), dbc.Input(value=0.2, id='filter-mito')],
                #               size='sm', style={'width':'80%', 'margin-left':'auto'}),
                html.Br(),
                html.H6("Apply These Parameters to:"),
                dbc.Button(html.H6('Apply These Parameters to:'), id='refresh-visual-samples-list', size='sm', color='white',
                           outline=False, style={'padding':'0px', 'cursor':'default', "visibility": "hidden"}),
                dbc.Checklist(id="checklist-visual-samples", persistence_type='session', inline=True, style={"textAlign":'center'}),
                html.Br(),
                dbc.Row(dbc.Col(dbc.Button("Run", id='run-visual-data', className="me-1", color='dark', outline=True),
                                width={'offset': 8, 'width': 4})),
            ],
        ),],
    #style={"width": "18rem"},
    ),
    label='Visual',
)

figParameter = dbc.Card(
    [dbc.CardBody(
        [
            html.H6('Visualization Parameters:'),
            dcc.Dropdown(id='visual-sample', placeholder="Select sample"),
        ],
    ),],
    style={"margin-top": "20px"},
)


tsneFig = dbc.Tab(
    [
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='tsne-gene-single', clearable=False, searchable=True), md=3),
            ],
            justify='evenly', style={'margin-top':'20px'},
        ),
        dcc.Graph(id='visual-tsne'),
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='tsne-gene-multi', clearable=False, searchable=True, multi=True), md=6),
            ],
            justify='evenly', style={'margin-top':'20px'},
        ),
        dcc.Graph(id='visual-tsne-gene'),
        html.Br(),
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='umap-gene-single', clearable=False, searchable=True), md=3),
            ],
            justify='evenly', style={'margin-top':'20px'},
        ),
        dcc.Graph(id='visual-umap'),
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(id='umap-gene-multi', clearable=False, searchable=True, multi=True), md=6),
            ],
            justify='evenly', style={'margin-top':'20px'},
        ),
        dcc.Graph(id='visual-umap-gene'),
    ],
    label='VisualFigure'
)


layout = html.Div(
    dbc.Row(
        [
            dbc.Col([dbc.Tabs([cluster_tab, ], ), figParameter], md=3,),
            dbc.Col(dbc.Tabs([tsneFig, ]), md=9),
        ]
    )
)


@callback(Output("checklist-visual-samples", 'options'),
          Input('meta-data2', 'data'),)
def refreshSamplesChecklist(metadata2):
    opt = [{'label': i, 'value': i}
     for i in metadata2['sampleInfo'].keys()
     if "filterParameter" in metadata2['sampleInfo'][i].keys()]
    return opt
