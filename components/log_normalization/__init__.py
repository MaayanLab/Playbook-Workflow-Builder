import anndata as ad
import pandas as pd
import scanpy as sc
from maayanlab_bioinformatics.normalization import log2_normalize
from components.gene_count_matrix import gene_count_matrix, anndata_from_path
from components.file import upsert_file

def log_normalize_gene_count_matrix(m):
  anndata_input = anndata_from_path(m['url'])

  # log2 normalization
  df = anndata_input.to_df()
  norm_data = log2_normalize(df, offset=1)

  with upsert_file('.csv') as f:
    norm_data.to_csv(f.file)

  return gene_count_matrix(f.url)

