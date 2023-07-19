import json
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import norm

def createbarchart(x, terms='Tissues'):
    x = pd.DataFrame.from_records(x)

    # Calculate p-values
    x['pvalue'] = 1 - norm.cdf(x['zscore'])

    # Sort by p-value in descending order and select top 10
    x = x.sort_values(by='pvalue', ascending=False).drop_duplicates(subset='pvalue').head(10)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x['pvalue'],
        y=x['term'],
        orientation='h',
        text=x.apply(lambda row: f"{row['term']} (p={row['pvalue']:.4f})", axis=1),
        textposition='inside',
        insidetextanchor='middle',
        textfont=dict(color='white'),
        marker=dict(line=dict(color='black', width=1))
    ))

    fig.update_layout(
        title=f"Top 10 {terms} with Highest Unique P-Values",
        height=600,
        yaxis=dict(showticklabels=False),
        autosize=False,  # Disable automatic resizing
        margin=dict(t=100),  # Add top margin for the title
        title_font=dict(size=15)  # Customize the title font size
    )

    bar_width = 0.8
    fig.update_traces(width=bar_width)

    return json.loads(fig.to_json())


