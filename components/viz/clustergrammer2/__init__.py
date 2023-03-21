import json
import pandas as pd
from components.data.gene_count_matrix import anndata_from_path
from clustergrammer2 import net

def clustergrammer2_from_df(url):
  df = anndata_from_path(url)
  net.load_df(df.to_df())
  net.cluster()
  data = dict(network=json.loads(net.export_net_json()))
  data['network']['order'] = {}
  data['network']['order']['row'] = 'rank'
  data['network']['order']['col'] = 'clust'
  return data

def clustergrammer2_from_gmt(gmt):
  df = pd.DataFrame({
    term: { gene: 1 for gene in record['set'] }
    for term, record in gmt.items()
  }).fillna(0)
  net.load_df(df)
  net.cluster()
  data = dict(network=json.loads(net.export_net_json()))
  data['network']['order'] = {}
  data['network']['order']['row'] = 'rank'
  data['network']['order']['col'] = 'clust'
  return data
