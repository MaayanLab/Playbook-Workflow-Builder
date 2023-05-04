import numpy as np
import anndata as ad
from components.core.file import fsspec_open_as_path, fsspec_open_as_iterator, upsert_file

def anndata_from_gctx(path):
  with fsspec_open_as_path(path) as fr:
    import h5py
    f = h5py.File(fr, 'r')
    return ad.AnnData(
      X=f['0']['DATA']['0']['MATRIX'],
      obs=f['0']['META']['0']['ROW'],
      var=f['0']['META']['0']['COL'],
    )

def anndata_from_gct(path):
  with fsspec_open_as_iterator(path, 'r') as fr:
    import pandas as pd
    version = fr.readline()
    shape = list(map(int, fr.readline().split('\t')))
    columns = fr.readline().split('\t')
    df = pd.read_csv(fr, sep='\t', header=columns)
    return ad.AnnData(
      X=df.iloc[-shape[1]:],
      obs=df.iloc[:-shape[1]],
    )

def anndata_from_path(path):
  ''' Read from a bunch of different formats, get an anndata file
  '''
  if path.endswith('.h5ad'):
    with fsspec_open_as_path(path, 'r') as fr:
      return ad.read_h5ad(fr)
  elif path.endswith('.csv'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return ad.read_text(fr, delimiter=',')
  elif path.endswith('.tsv'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return ad.read_text(fr, delimiter='\t')
  elif path.endswith('.txt') or path.endswith('.tab') or path.endswith('.data'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return ad.read_text(fr, delimiter=None)
  elif path.endswith('.xlsx'):
    with fsspec_open_as_path(path, 'r') as fr:
      return ad.read_excel(fr)
  elif path.endswith('.gctx'):
    return anndata_from_gctx(path)
  elif path.endswith('.gct'):
    return anndata_from_gct(path)
  elif path.endswith('.h5'):
    with fsspec_open_as_path(path, 'r') as fr:
      return ad.read_hdf(fr)
  else:
    raise NotImplementedError

def np_jsonifyable(x):
  x_ = x.astype('object')
  x_[np.isnan(x)] = 'nan'
  x_[np.isposinf(x)] = 'inf'
  x_[np.isneginf(x)] = '-inf'
  return x_.tolist()

def gene_count_matrix(url):
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
  with upsert_file('.h5ad') as f:
    d.write_h5ad(f.file)
  return gene_count_matrix(f.url)
