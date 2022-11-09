import qnorm
import scanpy as sc
from components.gene_count_matrix import gene_count_matrix, anndata_from_path
from components.file import upsert_file

def quantile_normalize_gene_count_matrix(m):
  anndata_input = anndata_from_path(m['url'])

  # quantile normalize
  df = anndata_input.to_df()
  norm_data = qnorm.quantile_normalize(df, axis=1)

  with upsert_file('.csv') as f:
    norm_data.to_csv(f.file)

  return gene_count_matrix(f.url)

