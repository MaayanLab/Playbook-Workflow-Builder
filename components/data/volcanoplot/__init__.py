import numpy as np
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import iplot
from IPython.display import display, Markdown
import sys, os
import json
from components.data.gene_count_matrix import anndata_from_path
import matplotlib.pyplot as plt

# ##### 2. Other libraries #####
# sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), 'core_scripts', 'shared', 'shared.py'))
# from shared import *


# def run(signature, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive'):

# 	# Loop through signature
# 	color = []
# 	text = []
# 	for index, rowData in signature.iterrows():

# 		# Text
# 		text.append('<b>'+index+'</b><br>Avg Expression = '+str(round(rowData['AveExpr'], ndigits=2))+'<br>logFC = '+str(round(rowData['logFC'], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData['P.Value'])+'<br>FDR = '+'{:.2e}'.format(rowData['adj.P.Val']))

# 		# Color
# 		if rowData['P.Value'] < pvalue_threshold:
# 		# if rowData['adj.P.Val'] < 0.05:
# 			if rowData['logFC'] < -logfc_threshold:
# 				color.append('blue')
# 			elif rowData['logFC'] > logfc_threshold:
# 				color.append('red')
# 			else:
# 				color.append('black')

# 		else:
# 			color.append('black')
	
# 	# results
# 	volcano_plot_results = {'x': signature['logFC'], 'y': -np.log10(signature['P.Value']), 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
# 	return volcano_plot_results

# #plot

# def plot(volcano_plot_results):
# 	spacer = ' '*50
# 	plot_2D_scatter(
# 		x=volcano_plot_results['x'],
# 		y=volcano_plot_results['y'],
# 		text=volcano_plot_results['text'],
# 		color=volcano_plot_results['color'],
# 		symmetric_x=True,
# 		xlab='log2FC',
# 		ylab='-log10P',
# 		title='<b>{volcano_plot_results[signature_label]} Signature | Volcano Plot</b>'.format(**locals()),
# 		labels=volcano_plot_results['signature_label'].split(' vs '),
# 		plot_type=volcano_plot_results['plot_type'],
# 		de_type='volcano'
# 	)

# 	return plot

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.offline import iplot
import json

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
    fig = plt.figure()
    plt.scatter(volcano_plot_results['x'], volcano_plot_results['y'], c=volcano_plot_results['color'])
    plt.xlabel('log2FC')
    plt.ylabel('-log10P')
    plt.title('<b>{}</b> Signature | Volcano Plot'.format(volcano_plot_results['signature_label']))
    plt.show()
    return fig

def createvolcano(gene_count_matrix):
	dataset = anndata_from_path(gene_count_matrix['url'])
	dataset = dataset.to_df()
	dataset = pd.DataFrame(dataset.values)
	data = run(dataset)
	fig = plot(data)

	return json.loads(fig.to_json())

