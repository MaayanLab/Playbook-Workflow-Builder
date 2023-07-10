import typing
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

def anndata_from_gene_count_matrix_and_metadata_matrix(gene_count_matrix: GeneCountMatrix, metadata_matrix: MetadataMatrix):
  df = anndata_from_file(gene_count_matrix).T
  df_meta = metadata_from_file(metadata_matrix)
  df.obs = df_meta
  with upsert_file('.h5ad') as f:
    df.write_h5ad(f.file)
  return anndata(f)
