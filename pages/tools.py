import dash_bootstrap_components as dbc
from dash import html

def generateTable(data, height, fontSize=10):
    table_header = [html.Thead(html.Tr([html.Th(i) for i in data[list(data)[0]].keys()]))]
    table_body = [html.Tbody(list(map(lambda i: html.Tr([html.Td(j) for j in i.values()]), data.values())))]
    table = dbc.Table(table_header + table_body, bordered=True, responsive=True, size='sm', hover=True,
                      style={'display': 'block', 'overflow': 'scroll', 'whiteSpace': 'nowrap',
                             'textAlign': 'center', 'height': height, 'fontSize': fontSize})
    return table



