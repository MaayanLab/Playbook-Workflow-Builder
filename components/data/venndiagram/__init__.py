import plotly.graph_objects as go
import json
import base64
from io import BytesIO
import matplotlib.pyplot as plt
from venn import venn
import sys

GeneSet = dict[str,list[str]]

def createvenn(*sets:GeneSet):
    genesets = {}
    for gs in sets:
        if set:
            genesets[gs["description"]] = set(gs["set"])


    plt.figure(figsize=(12,12), dpi=300)
    venn(genesets)

    buf = BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=300)
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close()

    fig = go.Figure()
    fig.add_layout_image(dict(
        source=f"data:image/png;base64,{img_base64}",
        xref="paper", yref="paper",
        x=0, y=1, sizex=1, sizey=1,
        xanchor="left", yanchor="top"
    ))

    fig.update_layout(
        width=1250,
        height=1250,
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )

    return json.loads(fig.to_json())