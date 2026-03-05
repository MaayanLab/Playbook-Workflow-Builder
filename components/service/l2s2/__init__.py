import io
import requests
import pandas as pd
import numpy as np
from components.core.file import File, upsert_file
import typing
from time import sleep
from components.data.gene_count_matrix import np_jsonifyable

class L2S2EnrichmentResults(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[str, int, float, bool, None]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def l2s2_enrichment_results(df: pd.DataFrame) -> L2S2EnrichmentResults:
  index = df.index.tolist()
  columns = df.columns.tolist()
  values = np_jsonifyable(df.fillna('').to_numpy())
  ellipses = [None,None]
  with upsert_file('.tsv') as f:
    df.to_csv(f.file, sep='\t')
  return dict(
    f,
    shape=df.shape,
    index=index,
    columns=columns,
    values=values,
    ellipses=ellipses,
  )

def fetch_enrichment_results(up_id: str, down_id: str) -> dict:
    req = requests.get(
        f'https://l2s2.maayanlab.cloud/enrichpair/download?datasetup={up_id}&datasetdown={down_id}',
        headers={'Accept': 'text/tab-separated-values'},
    )
    req.raise_for_status()
    df = pd.read_csv(io.BytesIO(req.content), sep='\t').drop(columns=['signatureCount'])
    return l2s2_enrichment_results(df)


def fetch_geneset_enrichment_results(id: str) -> dict:
    req = requests.get(
        f'https://l2s2.maayanlab.cloud/enrich/download?dataset={id}',
        headers={'Accept': 'text/tab-separated-values'},
    )
    req.raise_for_status()
    df = pd.read_csv(io.BytesIO(req.content), sep='\t').drop(columns=['geneSetHash'])
    return l2s2_enrichment_results(df)