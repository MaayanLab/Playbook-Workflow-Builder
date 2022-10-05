import requests
import pandas as pd
import scipy.stats as st
from functools import partial

def combine_pvalues(pvalues, method='fisher', select='pvalue'):
  ''' A helper for accessing this method via pd.agg which expects a scaler result
  '''
  statistic, pvalue = st.combine_pvalues(pvalues, method=method)
  return dict(statistic=statistic, pvalue=pvalue)[select]

def gtex_singleTissueEqtl(geneSymbol: str, datasetId: str='gtex_v8'):
  res = requests.get(
    'https://gtexportal.org/rest/v1/association/singleTissueEqtl',
    params=dict(
      format='json',
      geneSymbol=geneSymbol,
      datasetId=datasetId,
    )
  )
  gtex_results = res.json()['singleTissueEqtl']
  if len(gtex_results) == 0: raise Exception(f"No information for gene with identifier {geneSymbol} found in GTEx")
  gtex_combined_stouffer_statistic = (
    pd.DataFrame(gtex_results).groupby('tissueSiteDetailId')['pValue']
      .agg(partial(combine_pvalues, method='stouffer', select='statistic'))
      .to_frame('combined_stouffer_statistic')
      .reset_index()
      .sort_values('combined_stouffer_statistic', ascending=False)
  )
  gtex_combined_stouffer_statistic['group'] = gtex_combined_stouffer_statistic['tissueSiteDetailId'].apply(lambda name: name.split('_', maxsplit=1)[0])
  return [
    dict(tissue=row['tissueSiteDetailId'], zscore=row['combined_stouffer_statistic'])
    for _, row in gtex_combined_stouffer_statistic.iterrows()
  ]
