import scanpy as sc
import pandas as pd
from bokeh.embed import json_item

from components.data.gene_count_matrix import GeneCountMatrix, anndata_from_file
from components.viz.bokeh import interactive_circle_plot

def umap_transformation(gene_count_matrix: GeneCountMatrix):
  df = anndata_from_file(gene_count_matrix)

  # leiden coloring
  sc.pp.neighbors(df, use_rep='X')
  sc.tl.umap(df, n_components=3)
  sc.tl.leiden(df, key_added="leiden")

  df_for_visualization = pd.DataFrame({'leiden':df.obs['leiden'].values,
                                       'x':df.obsm['X_umap'][:,0],
                                       'y':df.obsm['X_umap'][:,1],
                                       'z':df.obsm['X_umap'][:,2]})

  plot = interactive_circle_plot(df_for_visualization, "UMAP-1", "UMAP-2", df_for_visualization.columns[0])

  return json_item(plot)
