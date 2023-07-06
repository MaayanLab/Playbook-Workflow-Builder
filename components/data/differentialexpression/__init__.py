#loadlibraries
from scipy import datasets
from sklearn.decomposition import PCA
import plotly.graph_objs as go
from plotly.offline import iplot, plot_mpl
import scipy.stats as ss
import warnings, sys, os
import pandas as pd
from IPython.display import display, Markdown
import chart_studio.plotly as py
import json
from components.data.gene_signature import signature_from_path
import numpy as np

from sklearn.decomposition import PCA
from sklearn.preprocessing import FunctionTransformer
import numpy as np

import pandas as pd

import chart_studio.plotly as py

import sys, os
from matplotlib import figure
import plotly.graph_objs as go
from plotly.offline import iplot
from IPython.display import display, Markdown


import pandas as pd
import sys
import os
from IPython.display import display, Markdown, HTML

##### 2. Other libraries #####
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), 'core_scripts', 'shared', 'shared.py'))
from shared import *




#run

def run(gene_signature, signature_label=''):
	sig = signature_from_path(gene_signature['url'])
	
	return sig
	#return {'signature': gene_signature.copy(), 'signature_label': signature_label}



#plot

def plot(signature_table_results):
	signature = signature_table_results['signature']
	signature_txt = signature.to_csv(sep='\t')
	# signature.index = ['<a href="http://www.genecards.org/cgi-bin/carddisp.pl?gene={x}" target="_blank">{x}</a>'.format(**locals()) for x in signature.index] # human
	signature.index = ['<a href="https://maayanlab.cloud/Harmonizome/gene/'+x.upper()+'" target="_blank">'+x+'</a>' for x in signature.index] # mouse
	signature.index = [index.replace('target="_blank">', 'target="_blank"><b>*').replace('</a>', '</b></a>') if rowData['adj.P.Val'] < 0.05 else index for index, rowData in signature.iterrows()] # mouse
	signature.index.name = 'Gene'
	signature['logFC'] = round(signature['logFC'], ndigits=2)
	signature['AveExpr'] = round(signature['AveExpr'], ndigits=2)
	html_table = signature.rename(columns={'gene_symbol': 'Gene', 'P.Value': 'P-value', 'adj.P.Val': 'FDR'}).drop(['t', 'B'], axis=1).sort_values('P-value').head(100).to_html(escape=False, classes='w-100')
	html_results = '<div style="max-height: 200px; overflow-y: scroll;">{}</div>'.format(html_table)
	display(HTML('<style>.w-100{width: 100%;}</style>'))
	# display(qgrid.show_grid(signature.rename(columns={'gene_symbol': 'Gene', 'P.Value': 'P-value', 'adj.P.Val': 'FDR'}).drop(['t', 'B'], axis=1), grid_options={'maxVisibleRows': 4}))
	display(HTML(html_results))




	# Add download button
	#download_button(signature_txt, 'Download Signature', '{signature_label}.txt'.format(**signature_table_results))
	# display(HTML('<textarea id="textbox" style="display: none;">{}</textarea> <button id="create">Download Results</button> <a download="signature.txt" id="downloadlink" style="display: none">Download</a>'.format(signature_txt)))
	# display(HTML('<script type="text/javascript">!function(){var e=null,t=document.getElementById("create"),n=document.getElementById("textbox");t.addEventListener("click",function(){var t,l,c=document.getElementById("downloadlink");c.href=(t=n.value,l=new Blob([t],{type:"text/plain"}),null!==e&&window.URL.revokeObjectURL(e),e=window.URL.createObjectURL(l)),c.click()},!1)}();</script>'))
	# download_button(signature_txt, 'Download Results', 'signature.txt')

	# Figure Legend
	#display(Markdown('** Table '+plot_counter('table')+' | Differential Expression Table.** The figure displays a browsable table containing the gene expression signature generated from a differential gene expression analysis. Every row of the table represents a gene; the columns display the estimated measures of differential expression. Links to external resources containing additional information for each gene are also provided'.format(**locals())))
