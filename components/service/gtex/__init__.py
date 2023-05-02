import requests
import pandas as pd

def gtex_resolve_genecode_id(geneId):
  req = requests.get(
    'https://gtexportal.org/api/v2/reference/gene',
    params=dict(
      geneId=geneId,
      pageSize=1,
      format='json',
    )
  )
  res = req.json()
  return res['data'][0]['gencodeId']

def gtex_gene_expression(geneSymbol: str, datasetId: str='gtex_v8'):
  gencodeId = gtex_resolve_genecode_id(geneSymbol)
  req = requests.get(
    'https://gtexportal.org/api/v2/expression/medianGeneExpression',
    params=dict(
      format='json',
      gencodeId=gencodeId,
      datasetId=datasetId,
    )
  )
  res = req.json()['data']
  if len(res) == 0: raise Exception(f"No information for gene with identifier {geneSymbol} found in GTEx")
  res = pd.DataFrame(res)
  res['zscore'] = (res['median'] - res['median'].mean()) / res['median'].std()
  res.rename({ 'tissueSiteDetailId': 'term' }, axis=1, inplace=True)
  return res[['term', 'zscore']].sort_values('zscore', ascending=False).to_dict(orient='records')
