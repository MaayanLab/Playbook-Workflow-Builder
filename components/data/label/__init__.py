import numpy as np
import pandas as pd
from components.core.file import File, upsert_file
from components.data.anndata import anndata, anndata_from_file

def jsonifyable_float(x):
  if type(x) != float: return x
  if np.isnan(x): return 'nan'
  if np.isposinf(x): return 'inf'
  if np.isneginf(x): return '-inf'
  return x

def get_metadata_from_anndata(file: File):
  ad = anndata_from_file(file)
  if not ad.obs.columns.shape[0]: ad.obs['index'] = ad.obs.index
  return {
    column: {
      row: jsonifyable_float(value)
      for row, value in record.items()
    }
    for column, record in ad.obs.to_dict().items()
  }

def unjsonifyable_float(x):
  ''' The inverse of jsonifyable_float '''
  if type(x) != str: return x
  return { 'nan': np.nan, 'inf': np.inf, '-inf': -np.inf }.get(x, x)

def is_missing(x):
  return not isinstance(x, str) and pd.isna(x)

def writable_column(column: pd.Series) -> pd.Series:
  ''' h5ad can't store an object column which mixes strings with other values, e.g. NaN from
  jsonifyable_float round-tripped as 'nan' among numbers, or labels among numbers. Numbers (and
  their 'nan'/'inf' encodings) become floats, anything else becomes strings (missing stays missing).
  '''
  if column.dtype != object: return column
  column = column.map(unjsonifyable_float)
  present = column[~column.map(is_missing)]
  if all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_)) for v in present):
    return column.astype(float)
  return column.map(lambda v: v if is_missing(v) else str(v))

def update_anndata_metadata(file: File, data: dict):
  ad = anndata_from_file(file)
  for col, vector in data.items():
    # the values are set on a copy of the column which can hold them, then the column is replaced:
    #  setting them in place keeps the column's type, which fails (pandas 3) for e.g. labels in an
    #  empty (float) metadata column, or values new to a categorical one
    column = ad.obs[col].astype(object) if col in ad.obs.columns else pd.Series(np.nan, index=ad.obs.index, dtype=object)
    # (only the samples in the matrix: a sample which isn't would be added as a row)
    keys = [k for k in vector.keys() if k in column.index]
    if vector and not keys:
      samples = list(map(str, ad.obs.index))
      raise ValueError(
        f"None of the samples given for {col!r} are in the matrix (e.g. {', '.join(map(repr, list(vector.keys())[:3]))}),"
        f" its samples are: {', '.join(samples[:50])}{f' (and {len(samples) - 50} more)' if len(samples) > 50 else ''}"
      )
    column.loc[keys] = [vector[k] for k in keys]
    ad.obs[col] = writable_column(column)
  with upsert_file('.h5ad', description=file.get('description')) as f:
    ad.write_h5ad(f.file)
  return anndata(f)
