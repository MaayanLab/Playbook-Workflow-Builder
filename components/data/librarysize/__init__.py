#loadlibraries
import json
from components.data.gene_count_matrix import anndata_from_file
import plotly.graph_objs as go
import pandas as pd
import plotly.graph_objects as go

# Assuming your gene count matrix is stored in a pandas DataFrame called 'gene_counts'
# The columns represent samples and the rows represent genes

def createlibrarysize(gene_count_matrix):
    
    cell_labels = gene_count_matrix["columns"]
    
    dataset = anndata_from_file(gene_count_matrix)
    dataset = dataset.to_df()
    dataset = pd.DataFrame(dataset.values)

    

    # Calculate the library sizes by summing the counts for each sample
    library_sizes = dataset.sum()

    # Create a bar graph using Plotly
    fig = go.Figure(data=go.Bar(x=cell_labels, y=library_sizes))

    # Customize the graph layout
    fig.update_layout(
        xaxis_title='Samples',
        yaxis_title='Library Size',
        title='Library Sizes of Samples'
    )

    # Add column labels as x-axis tick labels
    text = ['<b>{}</b><br>'.format(label) for label in cell_labels]

    fig.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(len(cell_labels))),
            ticktext=list(text)
        )
    )

    return json.loads(fig.to_json())

