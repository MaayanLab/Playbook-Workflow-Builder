from bokeh.embed import json_item
import pandas as pd
import scanpy as sc

from components.data.gene_count_matrix import GeneCountMatrix, anndata_from_file
from components.viz.bokeh import interactive_circle_plot

def pca_transformation(gene_count_matrix: GeneCountMatrix):
  df = anndata_from_file(gene_count_matrix)

  # PCA
  sc.pp.pca(df)

  # leiden coloring
  sc.pp.neighbors(df)
  sc.tl.leiden(df, key_added="leiden")

  df_for_visualization = pd.DataFrame({'leiden':df.obs['leiden'].values,
                                       'x':df.obsm['X_pca'][:,0],
                                       'y':df.obsm['X_pca'][:,1],
                                       'z':df.obsm['X_pca'][:,2],})

  plot = interactive_circle_plot(df_for_visualization, "PC-1", "PC-2", str(df_for_visualization.columns[0]))

  return json_item(plot)
