import plotly.graph_objects as go
import numpy as np
import json
from components.data.gene_signature import signature_from_file

def run(signature, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive'):
    # Loop through signature
    color = []
    text = []
    for index, rowData in signature.iterrows():
        # Text
        text.append('<b>'+str(index)+'</b><br>LogFC = '+str(round(rowData['LogFC'], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData['Pval'])+('<br>Adj P = '+'{:.2e}'.format(rowData['AdjPval']) if 'AdjPval' in rowData else ''))
        # Color
        if rowData['Pval'] < pvalue_threshold:
            if rowData['LogFC'] < -logfc_threshold:
                color.append('blue')
            elif rowData['LogFC'] > logfc_threshold:
                color.append('red')
            else:
                color.append('black')
        else:
            color.append('black')
    # Results
    volcano_plot_results = {'x': signature['LogFC'], 'y': -np.log10(signature['Pval']), 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results

def plot(volcano_plot_results):
    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=volcano_plot_results['x'],
        y=volcano_plot_results['y'],
        mode='markers',
        marker=dict(color=volcano_plot_results['color'])
    ))
    fig.update_layout(
        xaxis_title='log2FC',
        yaxis_title='-log10P',
        title_text='Volcano Plot',
        title_font=dict(family='Arial', size=16),
        
    )
    # Add annotation
    fig.add_annotation(
        x=1.0,  
        y=1.0,  
        text='Red: Upregulated<br>Blue: Downregulated',
        showarrow=False,
        font=dict(family='Arial', size=12),
        xref='paper',
        yref='paper'
    )
    return fig

def createvolcano(File):
    dataset = signature_from_file(File)
    data = run(dataset, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive')
    fig = plot(data)
    return json.loads(fig.to_json())
