import json
import pandas as pd
import plotly.express as px

def tissue_barplot(x):
  x = pd.DataFrame.from_records(x)
  fig = px.bar(
    x,
    y='zscore',
    x='term',
    orientation='v',
    title=f"Tissues with significant expression",
    height=1000,
  )
  return json.loads(fig.to_json())
