import numpy as np
from components.data.gene_count_matrix import gene_count_matrix, anndata_from_file
from components.core.file import upsert_file

def log2_normalize(x, offset=1.):
  return np.log2(x + offset)

def log_normalize_gene_count_matrix(m):
  df = anndata_from_file(m)

  # log2 normalization
  df.X = log2_normalize(df.X, offset=1.)

  with upsert_file('.h5ad') as f:
    df.write_h5ad(f.file)

  return gene_count_matrix(f)

