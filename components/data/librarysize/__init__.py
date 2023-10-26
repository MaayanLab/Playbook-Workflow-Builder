import json
from components.data.gene_count_matrix import anndata_from_file
import plotly.graph_objs as go
import pandas as pd
import plotly.graph_objects as go

# Code for making library size bar plot from gene count matrix

def createlibrarysize(gene_count_matrix):
    
    cell_labels = gene_count_matrix["columns"]
    
    dataset = anndata_from_file(gene_count_matrix)
    dataset = dataset.to_df()
    dataset = pd.DataFrame(dataset.values)

    # Calculate the library sizes by summing the counts for each sample
    library_sizes = dataset.sum()

    # Create a horizontal bar graph using Plotly
    fig = go.Figure(data=go.Bar(y=cell_labels, x=library_sizes, orientation='h'))

    # Customize the graph layout
    fig.update_layout(
        xaxis_title='Library Size',
        yaxis_title='Samples',
        title='Library Sizes of Samples',
    )

    return json.loads(fig.to_json())


# Code for making library size bar plot from AnnData file

def createlibrarysizefromanndata(anndata):
    
    # Read the AnnData file and create the AnnData object
    dataset = anndata_from_file(anndata)
    
    # Access the row names (sample names) from the AnnData object's index
    cell_labels = dataset.obs_names

    # Calculate the library sizes by summing the counts for each sample (row)
    library_sizes = dataset.X.sum(axis=1)

    # Create a horizontal bar graph using Plotly
    fig = go.Figure(data=go.Bar(y=cell_labels, x=library_sizes, orientation='h'))

    # Customize the graph layout
    fig.update_layout(
        xaxis_title='Library Size',
        yaxis_title='Samples',
        title='Library Sizes of Samples',
    )

    return json.loads(fig.to_json())