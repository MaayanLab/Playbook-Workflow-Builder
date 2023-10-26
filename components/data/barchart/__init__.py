import json
import pandas as pd
import plotly.graph_objects as go
from scipy.stats import norm

def createbarchart(x, terms='Tissues', pvalue_threshold=0.05):
    x = pd.DataFrame.from_records(x)
    x['zscore'] = x['zscore'].astype(float)

    # Calculate p-values from z-scores
    x['pvalue'] = norm.sf(abs(x['zscore'])) * 2

    # Sort by the absolute values of z-score in descending order
    x = x.loc[x['zscore'].abs().sort_values(ascending=False).index[:10][::-1]]

    # Define colors based on p-value significance
    x['bar_color'] = x['pvalue'].apply(lambda pval: 'blue' if pval <= pvalue_threshold else 'grey')

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=x['zscore'],
        y=x['term'],
        orientation='h',
        text=x.apply(lambda row: f"{row['term']} (p={row['pvalue']:.4f})", axis=1),
        textposition='inside',
        insidetextanchor='middle',
        textfont=dict(color='white'),
        marker=dict(color=x['bar_color'], line=dict(color='black', width=1))
    ))

    fig.update_layout(
        title=f"Top 10 {terms} with Highest Z-Scores",
        height=600,
        yaxis=dict(showticklabels=False),
        autosize=False,  # Disable automatic resizing
        margin=dict(t=100),  # Add top margin for the title
        title_font=dict(size=15)  # Customize the title font size
    )

    bar_width = 0.8
    fig.update_traces(width=bar_width)

    return json.loads(fig.to_json())

