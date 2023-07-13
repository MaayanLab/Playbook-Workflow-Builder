import scanpy as sc
import pandas as pd
from bokeh.embed import json_item
from sklearn.manifold import TSNE

from components.data.gene_count_matrix import GeneCountMatrix, anndata_from_file
from components.viz.bokeh import interactive_circle_plot

def tsne_transformation(gene_count_matrix: GeneCountMatrix):
  df = anndata_from_file(gene_count_matrix)
  df = df.transpose()

  # t-SNE
  sc.pp.pca(df)
  tsne = TSNE(n_components=3)
  df.obsm['X_tsne'] = tsne.fit_transform(df.obsm['X_pca'])

  # leiden coloring
  sc.pp.neighbors(df)
  sc.tl.leiden(df, key_added="leiden")

  df_for_visualization = pd.DataFrame({'leiden':df.obs['leiden'].values,
                                       'x':df.obsm['X_tsne'][:,0],
                                       'y':df.obsm['X_tsne'][:,1],
                                       'z':df.obsm['X_tsne'][:,2]})

  plot = interactive_circle_plot(df_for_visualization, "t-SNE-1", "t-SNE-2", df_for_visualization.columns[0])

  return json_item(plot)
