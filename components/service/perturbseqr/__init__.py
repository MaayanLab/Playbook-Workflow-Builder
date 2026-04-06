import io
import requests
import pandas as pd
import numpy as np
from components.core.file import File, upsert_file
import typing
from components.data.gene_count_matrix import np_jsonifyable
from time import sleep

class PerturbSeqrEnrichmentResults(File, typing.TypedDict):
  shape: typing.Tuple[int, int]
  index: typing.List[str]
  columns: typing.List[str]
  values: typing.List[typing.List[typing.Union[str, int, float, bool, None]]]
  ellipses: typing.Tuple[typing.Union[int, None], typing.Union[int, None]]

def perturbseqr_enrichment_results(df: pd.DataFrame, concat: int = None) -> PerturbSeqrEnrichmentResults:
  index = df.index.astype(str).tolist()
  columns = df.columns.tolist()
  values = np_jsonifyable(df.fillna('').to_numpy())
  ellipses = (concat,None)
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
    mimic_req = requests.get(
        f'https://perturbseqr.maayanlab.cloud/enrichpair/download?datasetup={up_id}&datasetdown={down_id}&sort=pvalue_mimic&maxTotal=1000',
        headers={'Accept': 'text/tab-separated-values'},
    )
    mimic_req.raise_for_status()
    mimic_df = pd.read_csv(io.BytesIO(mimic_req.content), sep='\t')
    reverse_req = requests.get(
        f'https://perturbseqr.maayanlab.cloud/enrichpair/download?datasetup={up_id}&datasetdown={down_id}&sort=pvalue_reverse&maxTotal=1000',
        headers={'Accept': 'text/tab-separated-values'},
    )
    reverse_req.raise_for_status()
    reverse_df = pd.read_csv(io.BytesIO(reverse_req.content), sep='\t')
    df = pd.concat([mimic_df,reverse_df]).drop_duplicates().drop(columns=['perturbationId', 'signatureCount']).reset_index(drop=True)
    ellipses = None
    if mimic_df.shape[0] == 1000:
      ellipses = 1000
    elif reverse_df.shape[0] == 1000:
      ellipses = df.shape[0]-reverse_df.shape[0]
    return perturbseqr_enrichment_results(df, ellipses)


def fetch_geneset_enrichment_results(id: str) -> dict:
    req = requests.get(
        f'https://perturbseqr.maayanlab.cloud/enrich/download?dataset={id}&maxTotal=1000',
        headers={'Accept': 'text/tab-separated-values'},
    )
    req.raise_for_status()
    df = pd.DataFrame(io.BytesIO(req.content), sep='\t').drop(columns=['geneSetHash']).reset_index(drop=True)
    return perturbseqr_enrichment_results(df)


def extract_perturbation_set(perturbseqr_results: PerturbSeqrEnrichmentResults, perturbation_type: typing.Literal['Drug','Gene'], direction: typing.Literal['Mimickers','Reversers']) -> dict[str,list]:
  perturbseqr_datasets = {
    'Drug': {"LINCS L1000 CP","Tahoe-100M","Microarrays CMap","NIBR DRUG-seq","SciPlex","DeepCover MoA","CREEDS Chem","RummaGEO Chem","Ginkgo Bioworks"},
    'Gene': {"LINCS L1000 XPR","Perturb Atlas Human","Perturb Atlas Mouse","CREEDS Gene","RummaGEO Gene","Replogle et al.","CM4AI"}
  }
  perturbseqr_columns = {
    'Mimickers': 'pvalueMimic',
    'Reversers': 'pvalueReverse',
  }
  df = pd.DataFrame(perturbseqr_results['values'],columns=perturbseqr_results['columns'])
  perturbation_df = df[df['dataset'].isin(perturbseqr_datasets[perturbation_type])]
  col = perturbseqr_columns[direction]
  perturbation_df[col] = pd.to_numeric(perturbation_df[col], errors='coerce')
  perturbation_set = perturbation_df[perturbation_df[col] < 0.05]['perturbation'].replace("",pd.NA).dropna().unique().tolist()
  if not perturbation_set:
    raise ValueError()

  return { 'description':'', 'set': perturbation_set }
