import typing
import numpy as np
import pandas as pd
from components.core.file import File, file_as_path, file_as_stream

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
  assert set(df.columns.tolist()) & { 'Pval', 'LogFC' }
  return df

def np_jsonifyable(x):
  x_ = x.astype('object')
  return x_.tolist()

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

def gmt_from_sig(sig: Signature):
  d = signature_from_file(sig).sort_values('Pval', ascending=True)
  up_250 = d[d['LogFC'] > 0].index[:250].tolist()
  down_250 = d[d['LogFC'] < 0].index[:250].tolist()
  return {
    f"{sig['description']} Up Genes": {
      'description': f"Top 250 up genes for the {sig['description']}",
      'set': up_250
    },
    f"{sig['description']} Down Genes": {
      'description': f"Top 250 down genes for the {sig['description']}",
      'set': down_250
    }
  }

def geneset_from_sig(sig: Signature, direction):
  d = signature_from_file(sig).sort_values('Pval', ascending=True)
  if direction == 'up':
    top_250 = d[d['LogFC'] > 0].index[:250].tolist()
  else:
    top_250 = d[d['LogFC'] < 0].index[:250].tolist()
  return {
    'description': f"{direction.capitalize()} genes from the {sig['description']}",
    'set': top_250
  }

def jsonifyable_float(x):
  if np.isnan(x): return 'nan'
  if np.isposinf(x): return 'inf'
  if np.isneginf(x): return '-inf'
  return x

def scored_genes_from_sig(sig: Signature):
  from scipy.stats import norm
  d = signature_from_file(sig)
  zscores = pd.Series(np.sign(d['LogFC']) * norm.ppf(1-d['Pval']), index=d.index)
  return [
    dict(term=term, zscore=jsonifyable_float(zscore))
    for term, zscore in zscores.to_dict().items()
    if zscore > 3 or zscore < -3
  ]
