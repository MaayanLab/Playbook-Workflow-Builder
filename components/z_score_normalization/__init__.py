from scipy.stats import zscore
import pandas as pd
import scanpy as sc
from components.gene_count_matrix import gene_count_matrix, anndata_from_path
from components.file import upsert_file

def z_score_normalize_gene_count_matrix(m):
  anndata_input = anndata_from_path(m['url'])

  # z-score normalization
  norm_data = anndata_input.to_df()
  normalized_df = pd.DataFrame(zscore(norm_data, axis=0), index=norm_data.index, columns=norm_data.columns)

  with upsert_file('.csv') as f:
    normalized_df.to_csv(f.file)

  return gene_count_matrix(f.url)

