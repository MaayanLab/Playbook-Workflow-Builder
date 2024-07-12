#loadlibraries
from sklearn.decomposition import PCA
import plotly.graph_objs as go
import pandas as pd
import json
from components.data.gene_count_matrix import anndata_from_file
import numpy as np
from sklearn.preprocessing import FunctionTransformer

def run(dataset):

    # Define a function to perform logCPM transformation
    logcpm_transform = FunctionTransformer(np.log1p, validate=True)

    # Apply logCPM transformation to your data
    expression_dataframe = logcpm_transform.transform(dataset)

    # Convert the NumPy array to a pandas DataFrame
    expression_dataframe = pd.DataFrame(expression_dataframe, columns=dataset.columns)

    # Run PCA
    pca = PCA(n_components=3)
    pca.fit(expression_dataframe)

    # Get Variance explained for each principal component
    var_explained = ['PC'+str((i+1))+' ('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]

    # Store variance explained with the PCA object
    pca.variance_explained = var_explained

    # Return PCA object
    return pca


# Code for making PCA with no metadata

def plotnometa(pca):
    fig = go.Figure(data=[go.Scatter3d(
        x=pca.components_[0],
        y=pca.components_[1],
        z=pca.components_[2],
        mode='markers'
    )])

    # Set the axis labels with variance explained and font size for readability
    fig.update_layout(
        scene=dict(
            xaxis_title=pca.variance_explained[0],
            yaxis_title=pca.variance_explained[1],
            zaxis_title=pca.variance_explained[2],
            xaxis_title_font=dict(size=10),
            yaxis_title_font=dict(size=10),
            zaxis_title_font=dict(size=10)
        )
    )

    return fig

def createpcanometa(gene_count_matrix):
    dataset = anndata_from_file(gene_count_matrix)
    dataset = dataset.to_df()
    dataset = pd.DataFrame(dataset.values)
    data = run(dataset)
    fig = plotnometa(data)

    return json.loads(fig.to_json())



# Code for making PCA with metadata

def createmetapcagraph(anndata):
    dataset = anndata_from_file(anndata)
    
    # Extract the relevant columns from .obs
    if 'Type: Control or Perturbation' in dataset.obs.columns:
        col = 'Type: Control or Perturbation'
        grp_ids = ['Control', 'Perturbation']
    else:
        col = dataset.obs.columns[0]
        grp_ids = dataset.obs[col].unique()

    ctrl_ids = dataset.obs[dataset.obs[col] == grp_ids[0]].index.tolist()
    ctrl_mask = [x in ctrl_ids for x in dataset.obs_names]
    case_ids = dataset.obs[dataset.obs[col] == grp_ids[1]].index.tolist()
    case_mask = [x in case_ids for x in dataset.obs_names]
    
    # Get PCA data
    pca = run(dataset.to_df().transpose())
    
    # Assign colors based on control and case masks
    colors = ['blue' if mask else 'red' for mask in (ctrl_mask + case_mask)]
    
    # Create the scatter plot with colored points
    fig = go.Figure(data=[go.Scatter3d(
        x=pca.components_[0],
        y=pca.components_[1],
        z=pca.components_[2],
        mode='markers',
        marker=dict(
            color=colors,
            size=4
        )
    )])

    # Set the axis labels with variance explained and font size for readability
    fig.update_layout(
        scene=dict(
            xaxis_title=pca.variance_explained[0],
            yaxis_title=pca.variance_explained[1],
            zaxis_title=pca.variance_explained[2],
            xaxis_title_font=dict(size=10),
            yaxis_title_font=dict(size=10),
            zaxis_title_font=dict(size=10)
        )
    )
    
    # Add legend/key that labels colors as control or perturbation
    fig.update_layout(
        annotations=[
            go.layout.Annotation(
                x=0.25,
                y=-0.1,
                showarrow=False,
                text='Control (Blue)',
                xref='paper',
                yref='paper',
                font=dict(color='blue')
            ),
            go.layout.Annotation(
                x=0.6,
                y=-0.1,
                showarrow=False,
                text='Perturbation (Red)',
                xref='paper',
                yref='paper',
                font=dict(color='red')
            )
        ]
    )

    return json.loads(fig.to_json())
