import typing
from components.core.file import File, upsert_file
from components.data.metabolite_count_matrix import MetaboliteCountMatrix, metabolite_count_matrix, metanndata_from_file
from components.data.metadata_matrix import MetadataMatrix, metadata_from_file

class MetAnnDataFile(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def metanndata(file: File) -> MetAnnDataFile:
  
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  assert file['filename'].endswith('h5ad')
  return metabolite_count_matrix(file)

def metanndata_from_metabolite_count_matrix_and_metadata_matrix(metabolite_count_matrix: MetaboliteCountMatrix, metadata_matrix: MetadataMatrix):
  df = metanndata_from_file(metabolite_count_matrix)
  df_meta = metadata_from_file(metadata_matrix)
  df.obs = df_meta
  with upsert_file('.h5ad', description=metabolite_count_matrix.get('description')) as f:
    df.write_h5ad(f.file)
  return metanndata(f)
