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

def update_anndata_metadata(file: File, data: dict):
  ad = anndata_from_file(file)
  for col, vector in data.items():
    # the values are set on a copy of the column which can hold them, then the column is replaced:
    #  setting them in place keeps the column's type, which fails (pandas 3) for e.g. labels in an
    #  empty (float) metadata column, or values new to a categorical one
    column = ad.obs[col].astype(object) if col in ad.obs.columns else pd.Series(np.nan, index=ad.obs.index, dtype=object)
    # (only the samples in the matrix: a sample which isn't would be added as a row)
    keys = [k for k in vector.keys() if k in column.index]
    column.loc[keys] = [vector[k] for k in keys]
    ad.obs[col] = column
  with upsert_file('.h5ad', description=file.get('description')) as f:
    ad.write_h5ad(f.file)
  return anndata(f)
