import numpy as np
from scipy.stats import zscore
from components.data.gene_count_matrix import GeneCountMatrix, gene_count_matrix, anndata_from_file
from components.core.file import upsert_file

def z_score_normalize_gene_count_matrix(m: GeneCountMatrix):
  df = anndata_from_file(m)

  # z-score normalization
  df.X = zscore(df.X, axis=0)

  # filter out genes with any null variances
  df = df[:, ~np.isnan(df.X).any(axis=0)]

  with upsert_file('.h5ad') as f:
    df.write_h5ad(f.file)

  return gene_count_matrix(f)

