import io
import requests
import pandas as pd
import anndata as ad
from components.core.file import upsert_file
from components.data.gene_count_matrix import gene_count_matrix

def fetch_samples(samples):
  if len(samples) > 100: raise RuntimeError('Too many samples, please select a smaller subset of samples')
  # get archs4 samples from archs4 api
  req = requests.post(
    'https://api.archs4.maayanlab.cloud/data/expression',
    headers={'Accept': 'text/tab-separated-values'},
    json={'geo_accession':samples},
  )
  req.raise_for_status()
  df = pd.read_csv(io.BytesIO(req.content), sep='\t', index_col=0)
  # make missing samples nan
  for sample in set(samples) - set(df.columns):
    df[sample] = float('nan')
  # register anndata
  d = ad.AnnData(df.T)
  with upsert_file('.h5ad') as f:
    d.write_h5ad(f.file)
  return gene_count_matrix(f)
