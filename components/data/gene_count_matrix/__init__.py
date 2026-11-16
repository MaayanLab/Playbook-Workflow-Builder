import typing
import numpy as np
import anndata as ad
import anndata.io as anndata_io
import pandas as pd
from components.core.file import File, file_as_seekable, file_as_stream, upsert_file

class GeneCountMatrix(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def anndata_from_gctx(file: File):
  ''' GCTX (the HDF5 GCT): the matrix (`0/DATA/0/matrix`, or `MATRIX`), the samples' ids &
  metadata (`0/META/COL`, or `0/META/0/COL`) & the genes' (`0/META/ROW`, or `0/META/0/ROW`).
  As in GCT, the columns are the samples. The matrix's orientation is that of the ids (cmapPy's
  files store it as samples x genes), strings are stored as bytes.
  '''
  import h5py
  def first(f, *paths):
    for p in paths:
      if p in f: return f[p]
    raise KeyError(f"Not a GCTX file: none of {', '.join(paths)} found")
  def values(dataset):
    v = dataset[()]
    return v.astype(str) if v.dtype.kind == 'S' else v
  def meta(group):
    return pd.DataFrame({k: values(group[k]) for k in group if k != 'id'}, index=pd.Index(values(group['id'])))
  with file_as_seekable(file) as fr, h5py.File(fr, 'r') as f:
    matrix = first(f, '0/DATA/0/matrix', '0/DATA/0/MATRIX')[()]
    samples = meta(first(f, '0/META/COL', '0/META/0/COL'))
    genes = meta(first(f, '0/META/ROW', '0/META/0/ROW'))
  if matrix.shape == (samples.shape[0], genes.shape[0]):
    X = matrix
  elif matrix.shape == (genes.shape[0], samples.shape[0]):
    X = matrix.T
  else:
    raise ValueError(f"The GCTX matrix {matrix.shape} doesn't match its {samples.shape[0]} column & {genes.shape[0]} row ids")
  return ad.AnnData(X=X, obs=samples, var=genes)

def anndata_from_gct(file: File):
  ''' GCT: `#1.2` (or `#1.3`), the dimensions (1.3: also the number of row & column metadata),
  then a table of genes x samples, after its row metadata columns (1.2: Description) & column
  metadata rows (1.3, which are the samples' metadata)
  '''
  import pandas as pd
  with file_as_stream(file, 'r') as fr:
    version = fr.readline()
    dims = list(map(int, fr.readline().split()))
    n_row_meta = dims[2] if len(dims) > 2 else 1
    n_col_meta = dims[3] if len(dims) > 3 else 0
    df = pd.read_csv(fr, sep='\t', index_col=0)
  values = df.iloc[n_col_meta:, n_row_meta:].apply(pd.to_numeric)
  d = ad.AnnData(values.T)
  if n_col_meta:
    d.obs = df.iloc[:n_col_meta, n_row_meta:].T.reindex(d.obs_names)
  return d

def anndata_from_file(file: File):
  ''' Read from a bunch of different formats, get an anndata file
  '''
  if file['filename'].endswith('.h5ad'):
    with file_as_seekable(file) as fr:
      return ad.read_h5ad(fr)
  elif file['filename'].endswith('.csv'):
    with file_as_stream(file, 'r') as fr:
      return anndata_io.read_text(fr, delimiter=',').transpose()
  elif file['filename'].endswith('.csv.gz'):
    import gzip
    with file_as_stream(file, 'rb') as fr_compressed:
      with gzip.open(fr_compressed, 'rt') as fr:
        return anndata_io.read_text(fr, delimiter=',').transpose()
  elif file['filename'].endswith('.tsv'):
    with file_as_stream(file, 'r') as fr:
      return anndata_io.read_text(fr, delimiter='\t').transpose()
  elif file['filename'].endswith('.tsv.gz'):
    import gzip
    with file_as_stream(file, 'rb') as fr_compressed:
      with gzip.open(fr_compressed, 'rt') as fr:
        return anndata_io.read_text(fr, delimiter='\t').transpose()
  elif file['filename'].endswith('.txt') or file['filename'].endswith('.tab') or file['filename'].endswith('.data'):
    with file_as_stream(file, 'r') as fr:
      try:
        return anndata_io.read_text(fr, delimiter=None).transpose()
      except Exception as e:
        raise Exception(f"Failed to read {file['filename']} as a matrix (please use .csv, or .tsv)") from e
  elif file['filename'].endswith('.xlsx'):
    with file_as_seekable(file) as fr:
      # (anndata's read_excel only reads paths, this is a file: the first sheet, first column the index)
      return ad.AnnData(pd.read_excel(fr, sheet_name=0, index_col=0).apply(pd.to_numeric)).transpose()
  elif file['filename'].endswith('.gctx'):
    return anndata_from_gctx(file)
  elif file['filename'].endswith('.gct'):
    return anndata_from_gct(file)
  elif file['filename'].endswith('.h5'):
    with file_as_seekable(file) as fr:
      try:
        return anndata_io.read_hdf(fr).transpose()
      except Exception as e:
        raise Exception(f"Failed to read {file['filename']} as a matrix (please use .gctx, or .h5ad)") from e
  else:
    raise NotImplementedError

def np_jsonifyable(x):
  if x.dtype != 'object':
    x_ = x.astype('object')
    x_[np.isnan(x)] = 'nan'
    x_[np.isposinf(x)] = 'inf'
    x_[np.isneginf(x)] = '-inf'
    return x_.tolist()
  return x.tolist()

def gene_count_matrix(file: File) -> GeneCountMatrix:
  ''' We'll preserve the file url but include various properties useful
  for visualization. If the file is invalid, reading it will fail.
  '''
  d = anndata_from_file(file)
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
  d = anndata_from_file(m)
  d = d.T
  with upsert_file('.h5ad') as f:
    d.write_h5ad(f.file)
  return gene_count_matrix(f)
