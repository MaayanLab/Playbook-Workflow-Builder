import typing
import pandas as pd
from components.core.file import File, upsert_file
from components.data.gene_count_matrix import GeneCountMatrix, gene_count_matrix, anndata_from_file
from components.data.metadata_matrix import MetadataMatrix, metadata_from_file

class AnnDataFile(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def anndata(file: File) -> AnnDataFile:
  
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  assert file['filename'].endswith('h5ad')
  return gene_count_matrix(file)

def merge_metadata(obs: pd.DataFrame, df_meta: pd.DataFrame):
  ''' Set the metadata's values on the samples (like `obs.update(df_meta)`: where the metadata has
  no value, the existing one is kept), each column keeping its own type: `update` writes into the
  existing columns' types, which pandas 3 refuses for e.g. text metadata in a new (float) column.
  '''
  meta = df_meta.reindex(obs.index)
  for c in df_meta.columns:
    if c not in obs.columns:
      obs[c] = meta[c]
    elif pd.api.types.is_numeric_dtype(meta[c]) and pd.api.types.is_numeric_dtype(obs[c]):
      obs[c] = meta[c].where(meta[c].notna(), obs[c])
    else:
      obs[c] = meta[c].astype(object).where(meta[c].notna(), obs[c].astype(object))

def anndata_from_gene_count_matrix_and_metadata_matrix(gene_count_matrix: GeneCountMatrix, metadata_matrix: MetadataMatrix):
  df = anndata_from_file(gene_count_matrix)
  df_meta = metadata_from_file(metadata_matrix)
  merge_metadata(df.obs, df_meta)
  with upsert_file('.h5ad', description=gene_count_matrix.get('description')) as f:
    df.write_h5ad(f.file)
  return anndata(f)
