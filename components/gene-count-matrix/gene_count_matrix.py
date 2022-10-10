import anndata as ad

def anndata_from_gctx(path):
  import h5py
  f = h5py.File(path, 'r')
  return ad.AnnData(
    X=f['0']['DATA']['0']['MATRIX'],
    obs=f['0']['META']['0']['ROW'],
    var=f['0']['META']['0']['COL'],
  )

def anndata_from_gct(path):
  import fsspec
  import pandas as pd
  with fsspec.open(path, 'r') as fr:
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
    return ad.read_h5ad(path)
  elif path.endswith('.csv'):
    return ad.read_text(path, delimiter=',')
  elif path.endswith('.tsv'):
    return ad.read_text(path, delimiter='\t')
  elif path.endswith('.txt') or path.endswith('.tab') or path.endswith('.data'):
    return ad.read_text(path, delimiter=None)
  elif path.endswith('.xlsx'):
    return ad.read_excel(path)
  elif path.endswith('.gctx'):
    return anndata_from_gctx(path)
  elif path.endswith('.gct'):
    return anndata_from_gct(path)
  elif path.endswith('.h5'):
    return ad.read_hdf(path)
  else:
    raise NotImplementedError

def gene_count_matrix(url):
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  d = anndata_from_path(url)
  return dict(
    url=url,
    shape=d.shape,
    head_columns=d.var_names[:5].tolist(), tail_columns=d.var_names[-5:].tolist(),
    head_index=d.obs_names[:5].tolist(), tail_index=d.obs_names[-5:].tolist(),
    head_values=d.X[:5, :5].tolist() + d.X[:5, -5:].tolist(), tail_values=d.X[-5:, :5].tolist() + d.X[-5:, -5:].tolist(),
  )
