import typing
import numpy as np
import pandas as pd
from components.core.file import File, file_as_path, file_as_stream
from components.data.gene_count_matrix import np_jsonifyable

class MetadataMatrix(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def metadata_from_file(file: File):
  ''' Read from a bunch of different formats, get a metadata table
  '''
  if file['filename'].endswith('.csv'):
    with file_as_stream(file, 'r') as fr:
      return pd.read_csv(fr, index_col=0)
  elif file['filename'].endswith('.tsv'):
    with file_as_stream(file, 'r') as fr:
      return pd.read_csv(fr, sep='\t', index_col=0)
  elif file['filename'].endswith('.txt') or file['filename'].endswith('.tab') or file['filename'].endswith('.data'):
    with file_as_stream(file, 'r') as fr:
      return pd.read_csv(fr, sep=None, index_col=0, engine='python')
  elif file['filename'].endswith('.xlsx'):
    with file_as_path(file, 'r') as fr:
      return pd.read_excel(fr, index_col=0)
  else:
    raise NotImplementedError

def metadata_matrix(file: File) -> MetadataMatrix:
  ''' Read the metadata file
  '''
  d = metadata_from_file(file)
  if d.shape[0] >= 10:
    top = 5
    bottom = 5
  elif d.shape[0] > 5:
    top = 5
    bottom = d.shape[0] - top
  else:
    top = d.shape[0] - 1
    bottom = 1

  if d.shape[1] != 1:
    raise Exception("Metadata file should contain exactly one column \
                    indicating the class to which each sample belongs.")
  col = d.columns[0]
  if len(d[col].unique()) != 2:
    raise Exception("Sample class column should only have two unique values, \
                    identifying the control group and perturbation group.")

  index = np.concatenate([d.index[:top], d.index[-bottom:]]).tolist()
  columns = d.columns.tolist()
  values = np_jsonifyable(np.concatenate([
    d.iloc[:top], d.iloc[-bottom:]
  ]))
  ellipses = [
    top if len(index) != d.shape[0] else None,
    None,
  ]
  return dict(
    file,
    shape=d.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )
