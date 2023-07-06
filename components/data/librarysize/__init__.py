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
from components.data.gene_count_matrix import anndata_from_path
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

##### 2. Other libraries #####
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))), 'core_scripts', 'shared', 'shared.py'))
import shared as s


import pandas as pd
import matplotlib.pyplot as plt

# Assuming your gene count matrix is stored in a pandas DataFrame called 'gene_counts'
# The columns represent samples and the rows represent genes


def createlibrarysize(gene_count_matrix):
	dataset = anndata_from_path(gene_count_matrix['url'])
	dataset = dataset.to_df()
	dataset = pd.DataFrame(dataset.values)

	# Calculate the library sizes by summing the counts for each sample
	library_sizes = dataset.sum()

	# Create a bar graph using Plotly
	fig = go.Figure(data=go.Bar(x=library_sizes.index, y=library_sizes))

	# Customize the graph layout
	fig.update_layout(
		xaxis_title='Samples',
		yaxis_title='Library Size',
		title='Library Sizes of Samples'
	)


	return json.loads(fig.to_json())

