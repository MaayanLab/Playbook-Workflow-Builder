import numpy as np
import anndata as ad
from components.core.file import file_as_path, file_as_stream, upsert_file





def anndata_from_path(path):
  ''' Read from a bunch of different formats, get an anndata file
  '''
  if path.endswith('.csv'):
    with file_as_stream(path, 'r') as fr:
      return ad.read_text(fr, delimiter=',')
  elif path.endswith('.tsv'):
    with file_as_stream(path, 'r') as fr:
      return ad.read_text(fr, delimiter='\t')
  elif path.endswith('.txt') or path.endswith('.tab') or path.endswith('.data'):
    with file_as_stream(path, 'r') as fr:
      return ad.read_text(fr, delimiter=None)
  elif path.endswith('.xlsx'):
    with file_as_path(path, 'r') as fr:
      return ad.read_excel(fr)
  else:
    raise NotImplementedError

def np_jsonifyable(x):
  x_ = x.astype('object')
  x_[np.isnan(x)] = 'nan'
  x_[np.isposinf(x)] = 'inf'
  x_[np.isneginf(x)] = '-inf'
  return x_.tolist()

def metabolite_count_matrix(url):
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  d = anndata_from_path(url)
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
    url=url,
    shape=d.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )

def transpose(m):
  d = anndata_from_path(m['url'])
  d = d.T
  with upsert_file('.txt') as f:
    d.write_text(f.file)
  return metabolite_count_matrix(f.url)
