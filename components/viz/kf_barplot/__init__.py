import json
import pandas as pd
import plotly.graph_objects as go


def barplot(tumor_expr, tissue_expr, terms='Tumor Expression'):
  tumor_expr = pd.DataFrame.from_records(tumor_expr)
  tumor_expr['tumor_zscore'] = (tumor_expr['TPM_median'] - tumor_expr['TPM_median'].mean()) / tumor_expr['TPM_median'].std()
  tumor_expr = tumor_expr.sort_values(by='tumor_zscore',ascending=False)

  tissue_expr = pd.DataFrame.from_records(tissue_expr)
  tissue_expr['tissue_zscore'] = tissue_expr['zscore']
  tissue_expr = tissue_expr.sort_values(by='tissue_zscore',ascending=True)

  fig = go.Figure(data=[
    go.Bar(name='tumors',x=tumor_expr['Disease'].head(10),y=tumor_expr['tumor_zscore'].head(10)),
    go.Bar(name='tissues',x=tissue_expr['term'].head(10),y=tissue_expr['tissue_zscore'].head(10))
  ])

  fig.update_layout(barmode='group')

  return json.loads(fig.to_json())