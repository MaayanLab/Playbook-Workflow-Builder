import numpy as np
from scipy.stats import zscore
from components.data.metabolite_count_matrix import MetaboliteCountMatrix, metabolite_count_matrix, metanndata_from_file
from components.core.file import upsert_file

def z_score_normalize_metabolite_count_matrix(m: MetaboliteCountMatrix):
  df = metanndata_from_file(m)
  
  # z-score normalization
  df.X = zscore(df.X, axis=0, ddof=1, nan_policy = 'omit')

  # filter out metabolites with any null variances
  df = df[:, ~np.isnan(df.X).any(axis=0)]

  with upsert_file('.h5ad') as f:
    df.write_h5ad(f.file)

  return metabolite_count_matrix(f)

