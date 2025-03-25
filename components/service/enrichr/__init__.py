import io
import json
import requests
import pandas as pd
import plotly.graph_objects as go
from urllib.parse import quote

def resolveGeneScoredLibraryUMAP(enrichrset):
  ''' Use the preprocessed Enrichr-Viz-Appyter coordinates to display UMAPs for the gene set libraries
  '''
  req = requests.get(f"https://raw.githubusercontent.com/MaayanLab/Enrichr-Viz-Appyter/master/Enrichr-Processed-Library-Storage/Clustered_Scatterplots/{quote(enrichrset['background'])}.csv")
  req.raise_for_status()
  df = pd.read_csv(io.BytesIO(req.content)).merge(pd.DataFrame.from_records(
    enrichrset['scored'], index=enrichrset['terms']
  ), right_index=True, left_on='term', how='outer')
  import scipy.stats as st
  df['p'] = st.norm.cdf(df['zscore'])
  df['p'] = df['p'].fillna(0.1)**20

  fig = go.Figure(data=[
    go.Scattergl(
      name=cluster,
      x=d['x'],
      y=d['y'],
      mode='markers',
      marker=dict(
        opacity=d['p'],
      ),
      hovertext=d['term'],
    )
    for cluster, d in df.groupby('cluster')
  ])
  # for cluster, d in df.groupby('cluster'):
  #   fig.update_traces(marker_opacity=d['p'], selector=dict(type='scattergl', name=cluster))

  return json.loads(fig.to_json())

def resolveGenesetLibraryUMAP(enrichrset):
  ''' Use the preprocessed Enrichr-Viz-Appyter coordinates to display UMAPs for the gene set libraries
  '''
  req = requests.get(f"https://raw.githubusercontent.com/MaayanLab/Enrichr-Viz-Appyter/master/Enrichr-Processed-Library-Storage/Clustered_Scatterplots/{quote(enrichrset['background'])}.csv")
  req.raise_for_status()
  df = pd.read_csv(io.BytesIO(req.content))
  d_in = df[df['term'].isin(enrichrset['terms'])]
  d_out = df[~df['term'].isin(enrichrset['terms'])]
  fig = go.Figure(data=[
    go.Scattergl(
      name='Relevant Sets',
      x=d_in['x'],
      y=d_in['y'],
      mode='markers',
      marker=dict(
        color='black',
      ),
      hovertext=d_in['term'],
    )
  ]+[
    go.Scattergl(
      name=cluster,
      x=d['x'],
      y=d['y'],
      mode='markers',
      opacity=0.1,
      hovertext=d['term'],
    )
    for cluster, d in d_out.groupby('cluster')
  ])
  return json.loads(fig.to_json())
