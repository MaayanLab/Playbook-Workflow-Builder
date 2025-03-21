import json
import typing
import pandas as pd
import scanpy as sc
import plotly.express as px
from sklearn.feature_extraction.text import TfidfVectorizer

class DescriptionSet(typing.TypedDict):
  set: list[str]
  description: typing.Optional[str]

def idf_umap(x: dict[str, DescriptionSet]):
  terms, sets = zip(*x.items())
  matrix = TfidfVectorizer(analyzer=lambda s: s['set']).fit_transform(sets)
  df = sc.AnnData(matrix)
  sc.pp.neighbors(df, use_rep='X')
  sc.tl.umap(df, n_components=2, random_state=42)
  sc.tl.leiden(df, key_added="leiden")

  umap = pd.DataFrame({
    'Leiden':df.obs['leiden'].values,
    'UMAP-1':df.obsm['X_umap'][:,0],
    'UMAP-2':df.obsm['X_umap'][:,1],
  }, index=terms).sort_values('Leiden')
  fig = px.scatter(
    umap,
    x='UMAP-1',
    y='UMAP-2',
    color='Leiden',
    hover_data=[umap.index],
    height=1000,
  )
  return json.loads(fig.to_json())
