import dash
import uuid
import dash_bootstrap_components as dbc
import dash_uploader as du
from collections import defaultdict
from dash import Input, Output, State, html, dcc, callback
from dash.exceptions import PreventUpdate

from pages.setting import DATA_PATH, PLOTLY_LOGO
from pages import uploadPage, processPage, testPage



app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)
server = app.server
du.configure_upload(app, folder=DATA_PATH,)

dropdownSinglecell = dbc.DropdownMenu(
            [
                dbc.DropdownMenuItem("Data", href="/singlecell/data",),
                dbc.DropdownMenuItem("Process", href="/singlecell/process",),
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



@app.callback(Output("page-content", "children"),
              Input("url", "pathname"))
def render_page_content(pathname):
    if pathname == "/":
        return dash.no_update
    elif pathname == "/singlecell/data":
        return uploadPage.layout
    elif pathname == "/singlecell/process":
        return processPage.layout
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


app.layout = html.Div([dcc.Location(id="url",  refresh=False),
                       dcc.Store('meta-data', storage_type='session', data=dict()),
                       sampleCanvas,
                       navbar,
                       dbc.Container(id='page-temp', className='pt-4'),
                       dbc.Container(id="page-content", className="pt-4"),
                       ])


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8050, debug=True)