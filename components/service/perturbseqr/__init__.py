import io
import requests
import pandas as pd
import numpy as np
import typing

def extract_perturbation_set(perturbseqr_signature: dict[str,str], perturbation_type: typing.Literal['Drug','Gene'], direction: typing.Literal['Mimickers','Reversers']) -> dict[str,list]:
  perturbseqr_datasets = {
    'Drug': {"LINCS L1000 CP","Tahoe-100M","Microarrays CMap","NIBR DRUG-seq","SciPlex","DeepCover MoA","CREEDS Chem","RummaGEO Chem","Ginkgo Bioworks"},
    'Gene': {"LINCS L1000 XPR","Perturb Atlas Human","Perturb Atlas Mouse","CREEDS Gene","RummaGEO Gene","Replogle et al.","CM4AI"}
  }
  perturbseqr_columns = {
    'Mimickers': {
      'sort':'pvalue_mimic',
      'column':'pvalueMimic'
    },'Reversers': {
      'sort':'pvalue_reverse',
      'column':'pvalueReverse'
    }
  }
  up_id = perturbseqr_signature['up_id']
  down_id = perturbseqr_signature['down_id']
  sort = perturbseqr_columns[direction]['sort']
  req = requests.get(
        f'https://perturbseqr.maayanlab.cloud/enrichpair/download?datasetup={up_id}&datasetdown={down_id}&sort={sort}&maxTotal=1000',
        headers={'Accept': 'text/tab-separated-values'},
    )
  req.raise_for_status()
  df = pd.read_csv(io.BytesIO(req.content), sep='\t')
  if df.empty:
    raise ValueError('No perturbation signatures.')
  perturbation_df = df[df['dataset'].isin(perturbseqr_datasets[perturbation_type])]
  col = perturbseqr_columns[direction]['column']
  perturbation_df[col] = pd.to_numeric(perturbation_df[col], errors='coerce')
  perturbation_set = perturbation_df[perturbation_df[col] < 0.05]['perturbation'].replace("",pd.NA).dropna().unique().tolist()
  if not perturbation_set:
    raise ValueError('No perturbations.')

  return { 'description':'', 'set': perturbation_set }
