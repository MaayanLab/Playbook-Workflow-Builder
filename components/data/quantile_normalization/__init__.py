import qnorm
from components.data.gene_count_matrix import GeneCountMatrix, gene_count_matrix, anndata_from_file
from components.core.file import upsert_file

def quantile_normalize_gene_count_matrix(m: GeneCountMatrix):
  df = anndata_from_file(m)

  # quantile normalize
  df.X = qnorm.quantile_normalize(df.X, axis=1)

  with upsert_file('.h5ad') as f:
    df.write_h5ad(f.file)

  return gene_count_matrix(f)

