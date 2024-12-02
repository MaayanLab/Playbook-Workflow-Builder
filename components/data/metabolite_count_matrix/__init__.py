import typing
import numpy as np
import anndata as ad
from components.core.file import File, file_as_path, file_as_stream, upsert_file
from components.data.gene_count_matrix import np_jsonifyable

class MetaboliteCountMatrix(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def metanndata_from_gctx(file: File):
  with file_as_path(file, 'r') as fr:
    import h5py
    f = h5py.File(fr, 'r')
    return ad.AnnData(
      X=f['0']['DATA']['0']['MATRIX'],
      obs=f['0']['META']['0']['ROW'],
      var=f['0']['META']['0']['COL'],
    )

def metanndata_from_gct(file: File):
  with file_as_stream(file, 'r') as fr:
    import pandas as pd
    version = fr.readline()
    shape = list(map(int, fr.readline().split('\t')))
    columns = fr.readline().split('\t')
    df = pd.read_csv(fr, sep='\t', header=columns)
    return ad.AnnData(
      X=df.iloc[-shape[1]:],
      obs=df.iloc[:-shape[1]],
    )

def metanndata_from_file(file: File):
  ''' Read from a bunch of different formats, get an anndata file
  '''
  print("I AM HERE")
  if file['filename'].endswith('.h5ad'):
    with file_as_path(file, 'r') as fr:
      return ad.read_h5ad(fr)
  elif file['filename'].endswith('.csv'):
    with file_as_stream(file, 'r') as fr:
      return ad.read_text(fr, delimiter=',').transpose()
  elif file['filename'].endswith('.tsv'):
    with file_as_stream(file, 'r') as fr:
      return ad.read_text(fr, delimiter='\t').transpose()
  elif file['filename'].endswith('.txt') or file['filename'].endswith('.tab') or file['filename'].endswith('.data'):
    with file_as_stream(file, 'r') as fr:
      return ad.read_text(fr, delimiter=None).transpose()
  elif file['filename'].endswith('.xlsx'):
    with file_as_path(file, 'r') as fr:
      return ad.read_excel(fr).transpose()
  elif file['filename'].endswith('.gctx'):
    return metanndata_from_gctx(file)
  elif file['filename'].endswith('.gct'):
    return metanndata_from_gct(file)
  elif file['filename'].endswith('.h5'):
    with file_as_path(file, 'r') as fr:
      return ad.read_hdf(fr).transpose()
  else:
    raise NotImplementedError

def metabolite_count_matrix(file: File) -> MetaboliteCountMatrix:
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  d = metanndata_from_file(file)
  # we want to see the gene count matrix preview in transposed form
  d = d.transpose()
  if d.shape[0] >= 10:
    top = 5
    bottom = 5
  elif d.shape[0] > 5:
    top = 5
    bottom = d.shape[0] - top
  else:
    top = d.shape[0] - 1
    bottom = 1
  if d.shape[1] >= 10:
    left = 5
    right = 5
  elif d.shape[1] > 5:
    left = 5
    right = d.shape[1] - left
  else:
    left = d.shape[1] - 1
    right = 1

  index = np.concatenate([d.obs_names[:top], d.obs_names[-bottom:]]).tolist()
  columns = np.concatenate([d.var_names[:left], d.var_names[-right:]]).tolist()
  values = np_jsonifyable(np.concatenate([
    np.concatenate([d.X[:top, :left], d.X[:top, -right:]], axis=1),
    np.concatenate([d.X[-bottom:, :left], d.X[-bottom:, -right:]], axis=1),
  ]))
  ellipses = [
    top if len(index) != d.shape[0] else None,
    left if len(columns) != d.shape[1] else None,
  ]
  return dict(
    file,
    shape=d.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )

def transpose(m: File):
  d = metanndata_from_file(m)
  d = d.T
  with upsert_file('.h5ad') as f:
    d.write_h5ad(f.file)
  return metabolite_count_matrix(f)
