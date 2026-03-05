import io
import requests
import pandas as pd
import anndata as ad
from components.core.file import upsert_file
from components.data.gene_count_matrix import gene_count_matrix
from components.data.metadata_matrix import metadata_matrix
from time import sleep

def fetch_samples(samples, species):
  if len(samples) > 100: raise RuntimeError('Too many samples, please select a smaller subset of samples')
  with requests.Session() as s:
    s.headers.update({"Accept": "application/json"})
    # get archs4 samples from archs4 api
    req = s.post(
      'https://maayanlab.cloud/sigpy/data/samples',
      json={'gsm_ids':samples, 'species':species}
    )
    req.raise_for_status()
    task_id = req.json()['task_id']
    # wait for resolve
    status = 'PENDING'
    while status in {"PENDING", "PROCESSING"}:
      req = s.get(f'https://maayanlab.cloud/sigpy/data/samples/status/{task_id}')
      req.raise_for_status()
      status = req.json()['status']
      sleep(0.05)
    # get expression download
    req = s.get(
      f'https://maayanlab.cloud/sigpy/data/samples/download/{task_id}',
      headers={"Accept": "text/tab-separated-values"}
    )
    req.raise_for_status()
  df = pd.read_csv(io.BytesIO(req.content), sep='\t', index_col=0, compression='zip')
  # make missing samples nan
  for sample in set(samples) - set(df.columns):
    df[sample] = float('nan')
  # register anndata
  d = ad.AnnData(df.T)
  with upsert_file('.h5ad') as f:
    d.write_h5ad(f.file)
  return gene_count_matrix(f)

def fetch_samples_meta(samples, species):
  if len(samples) > 100: raise RuntimeError('Too many samples, please select a smaller subset of samples')
  # get archs4 samples from archs4 api
  req = requests.post(
    'https://maayanlab.cloud/sigpy/meta/samplemeta',
    headers={'Accept': 'application/json'},
    json={'samples':samples, 'species':species}
  )
  req.raise_for_status()
  df = pd.read_json(io.BytesIO(req.content), dtype='string').T
  control_terms = {'wt', 'wildtype', 'control', 'cntrl', 'ctrl', 'uninfected', 'normal', 'untreated', 'unstimulated', 'shctrl', 'ctl', 'healthy', 'sictrl', 'sicontrol', 'ctr', 'wild', 'dmso', 'sint'}
  condition = df['characteristics'].map(lambda x: 'Control' if any(term in x.lower() for term in control_terms) else 'Perturbation')
  backup = (['Control']*(df.shape[0]//2))+(['Perturbation']*((df.shape[0])-(df.shape[0]//2)))
  df['Suggested Condition'] = condition if condition.nunique()==2 else backup
  # register anndata
  with upsert_file('.tsv') as f:
    df.to_csv(f.file, sep='\t')
  return metadata_matrix(f)
