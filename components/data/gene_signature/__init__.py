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
      return pd.read_csv(fr, index_col=0)
  elif file['filename'].endswith('.tsv'):
    with file_as_stream(file, 'r') as fr:
      return pd.read_csv(fr, sep='\t', index_col=0)
  elif file['filename'].endswith('.txt') or file['filename'].endswith('.tab') or file['filename'].endswith('.data'):
    with file_as_stream(file, 'r') as fr:
      return pd.read_csv(fr, sep=None, index_col=0, engine='python')
  elif file['filename'].endswith('.xlsx'):
    with file_as_path(file, 'r') as fr:
      return pd.read_excel(fr, index_col=0)
  else:
    raise NotImplementedError

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
  d = signature_from_file(sig)
  col = d.columns[0]
  comparison = col.split('vs.')[0].strip() + ' vs. ' + col.split('vs.')[1].split(':')[0].strip()
  up_250 = d.sort_values(by=col, ascending=False)[:250].index.tolist()
  down_250 = d.sort_values(by=col, ascending=True)[:250].index.tolist()
  return {
    f'{comparison} Up Genes': {
      'description': f'Top 250 up genes for {comparison}',
      'set': up_250
    },
    f'{comparison} Down Genes': {
      'description': f'Top 250 down genes for {comparison}',
      'set': down_250
    }
  }

def geneset_from_sig(sig: Signature, direction):
  d = signature_from_file(sig)
  col = d.columns[0]
  comparison = col.split('vs.')[0].strip() + ' vs. ' + col.split('vs.')[1].split(':')[0].strip()
  if direction == 'up':
    top_250 = d.sort_values(by=col, ascending=False)[:250].index.tolist()
  else:
    top_250 = d.sort_values(by=col, ascending=True)[:250].index.tolist()
  return {
    'description': f'{direction.capitalize()} Genes {comparison}',
    'set': top_250
  }

def scored_genes_from_sig(sig: Signature):
  from scipy.stats import zscore
  d = signature_from_file(sig)
  col = d.columns[0]
  scores = pd.Series(zscore(d[col]), index=d.index)
  return [
    dict(term=term, zscore=zscore)
    for term, zscore in scores.to_dict().items()
    if zscore > 3 or zscore < -3
  ]
