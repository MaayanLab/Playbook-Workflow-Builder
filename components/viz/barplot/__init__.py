import json
import pandas as pd
import plotly.express as px

def barplot(x, terms='Tissues'):
  x = pd.DataFrame.from_records(x)
  fig = px.bar(
    x,
    y='zscore',
    x='term',
    orientation='v',
    title=f"{terms} with significant expression",
    height=1000,
  )
  return json.loads(fig.to_json())
