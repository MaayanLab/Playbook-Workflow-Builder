import requests
import numpy as np
import pandas as pd
import scipy.stats as st
import typing as t
from components.data.anndata import AnnDataFile
from components.data.gene_count_matrix import GeneCountMatrix, anndata_from_file

targetranger_url = 'https://targetranger.maayanlab.cloud'

def targetscreener(m: t.Union[AnnDataFile, GeneCountMatrix], *, bg):
  ''' We prepare the incoming data for target screener by producing summary stats
  and translate the output into Scored[Gene]
  '''
  d = anndata_from_file(m)
  if 'Type: Control or Perturbation' in d.obs.columns:
    d = d[d.obs['Type: Control or Perturbation'] == 'Perturbation', :]
  # compute statistics which targetranger uses for the query
  n = d.X.shape[0]
  mean = d.X.mean(axis=0)
  std = d.X.std(axis=0)
  genes = pd.DataFrame({
    gene: { 'mean': mean[i], 'std': std[i] }
    for i, gene in enumerate(d.var.index)
  })
  # submit to targetranger API
  req = requests.post(
    f"{targetranger_url}/api/query_db_targets",
    json=dict(
      bg=bg,
      inputData=dict(genes=genes.to_dict(), n=n),
    ),
  )
  assert req.status_code == 200
  res = req.json()
  if len(res) == 0: raise Exception(f"No targets found")
  res = pd.DataFrame(res).dropna()
  # Translate response into output
  res = res.sort_values('adj_p', ascending=True)
  res['zscore'] = st.norm.ppf(1-res['adj_p'])
  res.rename({ 'gene': 'term' }, axis=1, inplace=True)
  res['zscore'].replace(np.nan, 'nan', inplace=True)
  res['zscore'].replace(np.inf, 'inf', inplace=True)
  res['zscore'].replace(-np.inf, '-inf', inplace=True)
  return res[['term', 'zscore']].to_dict(orient='records')
