import dash
import uuid
import json
import dash_uploader as du
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, callback, State

from .tools import generateTable

card = dbc.Card(
    [
        dbc.CardImg(src="/static/images/placeholder286x180.png", top=True),
        dbc.CardBody(
            [
                html.H4("Card title", className="card-title"),
                html.P(
                    "Some quick example text to build on the card title and "
                    "make up the bulk of the card's content.",
                    className="card-text",
                ),
                dbc.Button("Go somewhere", color="primary"),
            ]
        ),
    ],
    style={"width": "18rem"},
)


layout = dbc.Container(
    [
        #generateTable(META_DATA['sampleInfo'], '100px', 10),
        card,

    ]
)

