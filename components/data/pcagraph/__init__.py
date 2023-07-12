#loadlibraries
from sklearn.decomposition import PCA
import plotly.graph_objs as go
import sys, os
import pandas as pd
import json
from components.data.gene_count_matrix import anndata_from_file
import numpy as np
from sklearn.preprocessing import FunctionTransformer

def run(dataset):

# Assuming you have your data stored in a variable called `data`

	# Define a function to perform logCPM transformation
	logcpm_transform = FunctionTransformer(np.log1p, validate=True)

	# Apply logCPM transformation to your data
	expression_dataframe = logcpm_transform.transform(dataset)

	# Convert the NumPy array to a pandas DataFrame
	expression_dataframe = pd.DataFrame(expression_dataframe, columns=dataset.columns)


	# Run PCA
	pca=PCA(n_components=3)
	pca.fit(expression_dataframe)
	#pca.fit_transform(expression_dataframe)

	# Get Variance
	var_explained = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]


	# Return
	#pca_results = {'pca': pca, 'var_explained': var_explained, 'sample_metadata': dataset['sample_metadata'].loc[expression_dataframe.columns], 'color_by': color_by, 'color_type': color_type, 'nr_genes': nr_genes, 'normalization': normalization, 'signature_metadata': dataset.get('signature_metadata'), 'plot_type': plot_type}
	return pca

#Plot

def plot(pca):


	fig = go.Figure(data=[go.Scatter3d(
	x=pca.components_[0],
	y=pca.components_[1],
	z=pca.components_[2],
	mode='markers'
)])

	# Set the axis labels
	fig.update_layout(
		scene=dict(
			xaxis_title='PC1',
			yaxis_title='PC2',
			zaxis_title='PC3'
		)
	)

	return fig

def createpca(gene_count_matrix):
	dataset = anndata_from_file(gene_count_matrix)
	dataset = dataset.to_df()
	dataset = pd.DataFrame(dataset.values)
	data = run(dataset)
	fig = plot(data)

	return json.loads(fig.to_json())

