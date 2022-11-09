from bokeh.embed import json_item
import pandas as pd
from sklearn.decomposition import PCA
import scanpy as sc
import numpy as np

from components.gene_count_matrix import anndata_from_path
from components.bokeh import interactive_circle_plot

def pca_transformation(gene_count_matrix):
  anndata_input = anndata_from_path(gene_count_matrix['url'])
  df = anndata_input.to_df()
  df = df.transpose()

  # PCA
  pca = PCA()
  transformed_pca = pca.fit_transform(df)

  # leiden coloring
  leiden_df = sc.AnnData(df,dtype=np.float32)
  sc.pp.pca(leiden_df)
  sc.pp.neighbors(leiden_df)
  sc.tl.leiden(leiden_df, key_added="leiden")

  df_for_visualization = pd.DataFrame({'leiden':leiden_df.obs['leiden'].values,
                                       'x':transformed_pca[:,0],
                                       'y':transformed_pca[:,1],
                                       'z':transformed_pca[:,2],})

  plot = interactive_circle_plot(df_for_visualization, "PC-1", "PC-2", df_for_visualization.columns[0])

  return json_item(plot)
