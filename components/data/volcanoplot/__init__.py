import matplotlib.pyplot as plt
import numpy as np
import json
from components.data.gene_signature import signature_from_file

def run(signature, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive'):
    # Loop through signature
    color = []
    text = []
    for index, rowData in signature.iterrows():
        # Text
        text.append('<b>'+str(index)+'</b><br>Avg Expression = '+str(round(rowData['AveExpr'], ndigits=2))+'<br>logFC = '+str(round(rowData['logFC'], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData['P.Value'])+'<br>FDR = '+'{:.2e}'.format(rowData['adj.P.Val']))
        # Color
        if rowData['P.Value'] < pvalue_threshold:
            if rowData['logFC'] < -logfc_threshold:
                color.append('blue')
            elif rowData['logFC'] > logfc_threshold:
                color.append('red')
            else:
                color.append('black')
        else:
            color.append('black')
    # Results
    volcano_plot_results = {'x': signature['logFC'], 'y': -np.log10(signature['P.Value']), 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results

def plot(volcano_plot_results):
    fig, ax = plt.subplots()
    ax.scatter(volcano_plot_results['x'], volcano_plot_results['y'], c=volcano_plot_results['color'])
    ax.set_xlabel('log2FC')
    ax.set_ylabel('-log10P')
    ax.set_title('Volcano Plot', fontweight='bold')
    

    # Create a legend with color labels
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Blue: High logFC'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Red: Low logFC'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=8, label='Black: Not significant')
    ]
    ax.legend(handles=legend_elements, title='Legend', loc='upper right')
    
    plt.show()
    return fig

def createvolcano(File):
    dataset = signature_from_file(File)
    data = run(dataset, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive')
    fig = plot(data)
    return json.loads(fig.to_json())

