import numpy as np
import pandas as pd
import plotly.express as px
import dash
import plotly.graph_objects as go
from dash import Dash, html, dcc, Input, Output, callback, clientside_callback
from qdrant_client import QdrantClient
import wrenlab.normalize
from qdrant_client.conversions.common_types import Record
from qdrant_client.conversions.conversion import grpc_to_payload, json_to_value
from qdrant_client.models import Filter, FieldCondition, Range, PointStruct, HasIdCondition, PointIdsList, \
    VectorParams, IsEmptyCondition, PayloadField, \
    SearchRequest, RecommendRequest, TextIndexParams, TokenizerType, MatchText, \
    PayloadSchemaType, MatchValue, Distance, CreateAliasOperation, CreateAlias, OptimizersConfigDiff
from qdrant_client.uploader.grpc_uploader import payload_to_grpc
from joblib import Memory
import collections


# Template for querying biodata from qdrant vector search engine
# with small Dash-based webapp



# cachedir = '/data/tmp/qdrant_tmp'

#dash.register_page(__name__, path='/')


cx = QdrantClient(host="localhost", port=6333)

columns = cx.retrieve(collection_name="metadata",
            ids=[570],
            with_payload=True)[0].payload["columns"]

query_filter = Filter(
    must=[
        FieldCondition(
            key='Age',  # Condition based on values of `rand_number` field.
            range=Range(
                gte=25,
                lte=100
            ),
        ),
    ],
    must_not=[
        IsEmptyCondition(
            is_empty=PayloadField(key="TissueName")
        )
    ]
)

hits = cx.scroll(
        collection_name="GPL570",
        scroll_filter=query_filter,
        with_payload=True,
        with_vectors=False,
       limit=500)[0]    

def fetch_df(query_filter):
    hits = cx.scroll(
        collection_name="GPL570",
        scroll_filter=query_filter,
        with_payload=True,
        with_vectors=True,
        limit=500
    )[0]
    return pd.DataFrame.from_dict({hit.id: np.array(hit.vector) for hit in hits}, 
                                  orient="index",
                                  columns=columns)

df = fetch_df(query_filter)
essa_df = fetch_df(query_filter)

def process_df():
    df_mtx = fetch_df(query_filter)
    
    df_mtx.reset_index(inplace=True)
    df_mtx = df_mtx.rename(columns={'index': 'GSM'})
    df_mtx.set_index('GSM', inplace=True)
    return df_mtx


mtx_final_df = process_df()


def tis_id():
    columns_tis = []
    rec_id = []

    for hit in hits:
        rec_id.append(hit.id)
        tis = hit.payload['TissueName']
        columns_tis.append(tis)
    zipped = zip(rec_id, columns_tis)
    essa_dict = dict(zipped)
    return pd.DataFrame.from_dict(essa_dict, orient='index')

df_tis = tis_id()
df_tis.rename(columns = {0:"TissueName"}, inplace=True)
tis_col = df_tis.pop('TissueName')    
df.insert(0, 'TissueName', tis_col)  


def sex_id():

    l_hit = []
    columns_sex = []
    rec_id = []

    for hit in hits:
        rec_id.append(hit.id)
        l_hit.append(hit)
    for i in l_hit:
        i = list(i)
        sex = (i[1][1]['Sex'])
        columns_sex.append(sex)
    zipped = zip(rec_id, columns_sex)
    sex_dict = dict(zipped)
    return pd.DataFrame.from_dict(sex_dict, orient='index')

df_sex = sex_id()
df_sex.rename(columns = {0:"Sex"}, inplace=True)
sex_col = df_sex.pop('Sex')    
df.insert(1, 'Sex', sex_col)    


def age_id():
    columns_age = []
    rec_id = []

    for hit in hits:
        rec_id.append(hit.id)
        age = hit.payload['Age']
        columns_age.append(age)           
    zipped = zip(rec_id, columns_age)
    age_dict = dict(zipped)
    return pd.DataFrame.from_dict(age_dict, orient='index')

df_age = age_id()
df_age.rename(columns = {0:"Age"}, inplace=True)
age_col = df_age.pop('Age')    
df.insert(1, 'Age', age_col)    
age_bins_col = pd.cut(x=df['Age'], bins=[10, 19, 29, 39, 49, 59, 69, 79, 89, 99])    # Not the best idea lol
df.insert(2, 'Age_bins', age_bins_col)
age_deciles_col = pd.cut(x=df['Age'], bins=[10, 19, 29, 39, 49, 59, 69, 79, 89, 99], 
                         labels=['10s', '20s', '30s', '40s', '50s', '60s', '70s', '80s', '90s'])
df.insert(3, 'Decades', age_deciles_col)
df = df.set_index(['TissueName', 'Age', 'Age_bins', 'Decades', 'Sex'])

df = df.dropna()
df = df.T

final_df = wrenlab.normalize.quantile(df)

def fn():
    path = "data/GPL570.map.tsv"
    df = pd.read_table(path)
    o = {}
    for probe_id, symbols, gene_id in df.dropna().to_records(index=False):
        o[probe_id] = symbols.split(" /// ")[0]
    return o


M = fn()

final_df = final_df.T
final_df = final_df.groupby(M, axis=1).mean()

decades = list(final_df.index.get_level_values('Decades').unique().sort_values())
tissues = list(final_df.index.get_level_values('TissueName').unique())

mean_age = final_df.groupby('Decades').mean()
mean_tissue = final_df.groupby('TissueName').mean()


# Dropdown options for age/tissue

tis_age_opt = {
    'Age': decades,
    'Tissue': tissues
}


layout = html.Div(
    children=[
        dcc.Graph(
            id='graph',
        ),
        dcc.Store(
            id='clientside-figure-store',
        ),
        'Gene',
        dcc.Dropdown(
            options=[{'label': i, 'value': i} for i in final_df.columns],
            value='A1BG',
            placeholder="Select a gene",
            clearable=True,
            id='gene_dropdown',
        ),
        'Summarize expression by:',
        dcc.Dropdown(
            options=[{'label': i, 'value': i} for i in tis_age_opt.keys()],
            value='Tissue',
            placeholder="Select Age/Tissue",
            clearable=True,
            id='age_tissue',
        ),
        'Download expression matrix',
        dcc.Dropdown(
            options=[{'label': 'GPL570', 'value': 'final_df'}],
            value='final_df',
            clearable=False,
            id='dow_matrix',
        ),
        html.Button("Download matrix csv", id="btn_csv"),
        dcc.Download(id="download_df_csv"),
    ]
)

@callback(
    Output('clientside-figure-store', 'data'),
    Input('gene_dropdown', 'value'),
    Input('age_tissue', 'value')
)
def update_store_data(sel_gene, sel_age_tis):
    if sel_age_tis == 'Tissue':
        px_df = mean_tissue
    else:
        px_df = mean_age
    
    sd = px_df[sel_gene].std()
    
    print(sd)

    fig = px.histogram(px_df, x=sel_gene, y=px_df.index)

    fig.update_layout(title='Mean expression histogram', height=600)
    
    fig.update_yaxes(title='Tissue')
    fig.update_xaxes(title='Mean expression')
    
    return fig

@callback(
        Output('download_df_csv', 'data'),
        Input('btn_csv', 'n_clicks'),
        prevent_initial_call=True
)
def func(n_clicks):
    return dcc.send_data_frame(mtx_final_df.to_csv, "exp_matrix.csv")


clientside_callback(
      """
    function(figure, scale) {
        if(figure === undefined) {
            return {'data': [], 'layout': {}};
        }
        const fig = Object.assign({}, figure, {
            'layout': {
                ...figure.layout,
                'yaxis': {
                    ...figure.layout.yaxis, type: scale
                }
             }
        });
        return fig;
    }
    """,
    Output('graph', 'figure'),
    Input('clientside-figure-store', 'data')
)

