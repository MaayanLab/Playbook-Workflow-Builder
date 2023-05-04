import numpy as np
import pandas as pd
from components.core.file import fsspec_open_as_path, fsspec_open_as_iterator

def signature_from_path(path):
  ''' Read from a bunch of different formats, get a metadata table
  '''
  if path.endswith('.csv'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return pd.read_csv(fr, index_col=0)
  elif path.endswith('.tsv'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return pd.read_csv(fr, sep='\t', index_col=0)
  elif path.endswith('.txt') or path.endswith('.tab') or path.endswith('.data'):
    with fsspec_open_as_iterator(path, 'r') as fr:
      return pd.read_csv(fr, sep=None, index_col=0, engine='python')
  elif path.endswith('.xlsx'):
    with fsspec_open_as_path(path, 'r') as fr:
      return pd.read_excel(fr, index_col=0)
  else:
    raise NotImplementedError

def np_jsonifyable(x):
  x_ = x.astype('object')
  return x_.tolist()

def gene_signature(url):
  d = signature_from_path(url)
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
    url=url,
    shape=d.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )

def gmt_from_sig(sig):
  d = signature_from_path(sig['url'])
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

def geneset_from_sig(sig, direction):
  d = signature_from_path(sig['url'])
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