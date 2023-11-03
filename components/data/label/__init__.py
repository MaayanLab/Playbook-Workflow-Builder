import numpy as np
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
    keys = vector.keys()
    ad.obs.loc[keys, col] = [vector[k] for k in keys]
  with upsert_file('.h5ad', description=file.get('description')) as f:
    ad.write_h5ad(f.file)
  return anndata(f)
