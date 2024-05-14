import json
import typing
import pandas as pd
from umap import UMAP
import plotly.express as px
from sklearn.feature_extraction.text import TfidfVectorizer

class DescriptionSet(typing.TypedDict):
  set: list[str]
  description: typing.Optional[str]

def idf_umap(x: dict[str, DescriptionSet]):
  terms, sets = zip(*x.items())
  matrix = TfidfVectorizer(analyzer=lambda s: s['set']).fit_transform(sets)
  umap = pd.DataFrame(UMAP(random_state=42).fit_transform(matrix), columns=['UMAP-1', 'UMAP-2'], index=terms)
  fig = px.scatter(
    umap,
    x='UMAP-1',
    y='UMAP-2',
    hover_data=[umap.index],
    height=1000,
  )
  return json.loads(fig.to_json())
