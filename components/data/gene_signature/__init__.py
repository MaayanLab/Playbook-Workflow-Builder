import typing
import numpy as np
import pandas as pd
from components.core.file import File, file_as_path, file_as_stream
from components.data.gene_count_matrix import np_jsonifyable

class Signature(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[int, typing.Literal['nan'], typing.Literal['inf'], typing.Literal['-inf']]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def signature_from_file(file: File):
  ''' Read from a bunch of different formats, get a metadata table
  '''
  if file['filename'].endswith('.csv'):
    with file_as_stream(file, 'r') as fr:
      df = pd.read_csv(fr, index_col=0)
  elif file['filename'].endswith('.tsv'):
    with file_as_stream(file, 'r') as fr:
      df = pd.read_csv(fr, sep='\t', index_col=0)
  elif file['filename'].endswith('.txt') or file['filename'].endswith('.tab') or file['filename'].endswith('.data'):
    with file_as_stream(file, 'r') as fr:
      df = pd.read_csv(fr, sep=None, index_col=0, engine='python')
  elif file['filename'].endswith('.xlsx'):
    with file_as_path(file, 'r') as fr:
      df = pd.read_excel(fr, index_col=0)
  else:
    raise NotImplementedError
  assert set(df.columns.tolist()) & { 'Pval', 'AdjPval', 'LogFC' }
  return df

def gene_signature(file: File) -> Signature:
  d = signature_from_file(file)
  if d.shape[0] >= 10:
    top = 5
    bottom = 5
  elif d.shape[0] > 5:
    top = 5
    bottom = d.shape[0] - top
  else:
    top = d.shape[0] - 1
    bottom = 1

  index = np.concatenate([d.index[:top], d.index[-bottom:]]).tolist()
  columns = d.columns.tolist()
  values = np_jsonifyable(np.concatenate([
    d.iloc[:top], d.iloc[-bottom:]
  ]))
  ellipses = [
    top if len(index) != d.shape[0] else None,
    None,
  ]
  return dict(
    file,
    shape=d.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )

def resonable_geneset_threshold(d, up_down=False):
  ''' Maybe this should be manually specified by the user but for now we try several sane choices
  and pick the one which results in a geneset >0 but closest to 250 up/down genes.
  '''
  candidates = [
    (abs(n_up - 250) + abs(n_down - 250), col, thresh)
    for col, thresh in [('Pval', 0.05),('Pval', 0.01),('Pval', 0.001),('AdjPval', 0.05),('AdjPval', 0.01),('AdjPval', 0.001)]
    for n_up, n_down in ((((d['LogFC']>0) & (d[col]<thresh)).sum(), ((d['LogFC']<0) & (d[col]<thresh)).sum()),)
    if ((n_up > 0 and n_down > 0) if up_down else (n_up+n_down)>0)
  ]
  if not candidates:
    return ('Pval', 0.05)
  _, col, thresh = min(candidates)
  return col, thresh

def gmt_from_sig(sig: Signature):
  d = signature_from_file(sig).sort_values('Pval', ascending=True)
  col, thresh = resonable_geneset_threshold(d, up_down=True)
  up = d[(d['LogFC'] > 0) & (d[col] < thresh)].index.tolist()
  down = d[(d['LogFC'] < 0) & (d[col] < thresh)].index.tolist()
  return {
    f"{sig.get('description', 'signature')} Up Genes": {
      'description': f"Significant up genes for the {sig.get('description', 'signature')}",
      'set': up
    },
    f"{sig.get('description', 'signature')} Down Genes": {
      'description': f"Significant down genes for the {sig.get('description', 'signature')}",
      'set': down
    }
  }

def geneset_from_sig(sig: Signature, direction):
  d = signature_from_file(sig).sort_values('Pval', ascending=True)
  col, thresh = resonable_geneset_threshold(d, up_down=True)
  if direction == 'up':
    top = d[(d['LogFC'] > 0) & (d[col] < thresh)].index.tolist()
  else:
    top = d[(d['LogFC'] < 0) & (d[col] < thresh)].index.tolist()
  return {
    'description': f"Significant {direction} genes from the {sig.get('description', 'signature')}",
    'set': top
  }

def jsonifyable_float(x):
  if np.isnan(x): return 'nan'
  if np.isposinf(x): return 'inf'
  if np.isneginf(x): return '-inf'
  return x

def scored_genes_from_sig(sig: Signature):
  from scipy.stats import norm
  d = signature_from_file(sig)
  col, thresh = resonable_geneset_threshold(d)
  d = d[d[col]<thresh]
  zscores = pd.Series(np.sign(d['LogFC']) * norm.ppf(1-d[col]), index=d.index)
  return [
    dict(term=term, zscore=jsonifyable_float(zscore))
    for term, zscore in zscores.to_dict().items()
  ]
