import requests
import numpy as np
import pandas as pd
import scipy.stats as st
from components.data.gene_count_matrix import anndata_from_path

targetranger_url = 'https://targetranger.maayanlab.cloud'

# https://github.com/MaayanLab/TargetRanger/blob/main/pages/api/query_db_targets.js#L18-L27
bg_lookup = {
  'ARCHS4': 0,
  'GTEx_transcriptomics': 1,
  'Tabula_Sapiens': 2,
  'CCLE_transcriptomics': 3,
  'HPM': 4,
  'HPA': 5,
  'GTEx_proteomics': 6,
  'CCLE_proteomics': 7,
}

def targetscreener(url, *, bg):
  ''' We prepare the incoming data for target screener by producing summary stats
  and translate the output into Scored[Gene]
  '''
  d = anndata_from_path(url)
  # compute statistics which targetranger uses for the query
  n = d.X.shape[1]
  mean = d.X.mean(axis=1)
  std = d.X.std(axis=1)
  genes = pd.DataFrame({
    gene: { 'mean': mean[i], 'std': std[i] }
    for i, gene in enumerate(d.obs.index)
  })
  # submit to targetranger API
  req = requests.post(
    f"{targetranger_url}/api/query_db_targets",
    json=dict(
      bg=bg_lookup[bg],
      inputData=dict(genes=genes.to_dict(), n=n),
    ),
  )
  assert req.status_code == 200
  res = req.json()
  if len(res) == 0: raise Exception(f"No targets found")
  res = pd.DataFrame(res)
  # Translate response into output
  res = res.sort_values('adj_p', ascending=True)
  res['zscore'] = st.norm.ppf(1-res['adj_p'])
  res.rename({ 'gene': 'term' }, axis=1, inplace=True)
  res['zscore'].replace(np.nan, 'nan', inplace=True)
  res['zscore'].replace(np.inf, 'inf', inplace=True)
  res['zscore'].replace(-np.inf, '-inf', inplace=True)
  return res[['term', 'zscore']].to_dict(orient='records')
