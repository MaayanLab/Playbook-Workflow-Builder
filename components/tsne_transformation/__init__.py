from MulticoreTSNE import MulticoreTSNE as TSNE
import scanpy as sc
import numpy as np
import pandas as pd
import multiprocessing as mp
from bokeh.embed import json_item

from components.gene_count_matrix import anndata_from_path
from components.bokeh import interactive_circle_plot

def tsne_transformation(gene_count_matrix):
  anndata_input = anndata_from_path(gene_count_matrix['url'])
  df = anndata_input.to_df()
  df = df.transpose()

  # t-SNE
  n_jobs=mp.cpu_count()
  tsne = TSNE(n_jobs = n_jobs, n_components=3)
  transformed_tsne = tsne.fit_transform(df)

  # leiden coloring
  leiden_df = sc.AnnData(df,dtype=np.float32)
  sc.pp.pca(leiden_df)
  sc.pp.neighbors(leiden_df)
  sc.tl.leiden(leiden_df, key_added="leiden")

  df_for_visualization = pd.DataFrame({'leiden':leiden_df.obs['leiden'].values,
                                       'x':transformed_tsne[:,0],
                                       'y':transformed_tsne[:,1],
                                       'z':transformed_tsne[:,2]})

  plot = interactive_circle_plot(df_for_visualization, "t-SNE-1", "t-SNE-2", df_for_visualization.columns[0])

  return json_item(plot)
