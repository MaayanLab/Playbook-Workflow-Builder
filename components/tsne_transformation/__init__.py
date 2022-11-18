import scanpy as sc
import pandas as pd
import multiprocessing as mp
from bokeh.embed import json_item
from sklearn.manifold import TSNE

from components.gene_count_matrix import anndata_from_path
from components.bokeh import interactive_circle_plot

def tsne_transformation(gene_count_matrix):
  df = anndata_from_path(gene_count_matrix['url'])
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
