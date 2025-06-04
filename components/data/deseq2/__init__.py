from components.data.gene_count_matrix import anndata_from_file
from components.data.gene_signature import gene_signature
from components.core.file import upsert_file
from maayanlab_bioinformatics.dge import deseq2_differential_expression
import pandas as pd

# Function for computing signatures with characteristic direction
def deseq2(anndata):
  if 'Type: Control or Perturbation' in anndata.obs.columns:
    col = 'Type: Control or Perturbation'
    grp_ids = ['Control', 'Perturbation']
  else:
    col = anndata.obs.columns[-1]
    grp_ids = anndata.obs[col].unique()
    assert len(grp_ids) == 2, "Expected two possible groups"

  ctrl_ids = anndata.obs[anndata.obs[col] == grp_ids[0]].index.tolist()
  ctrl_mask = [x in ctrl_ids for x in anndata.obs_names]
  case_ids = anndata.obs[anndata.obs[col] == grp_ids[1]].index.tolist()
  case_mask = [x in case_ids for x in anndata.obs_names]

  ctrl_df = pd.DataFrame(
    data = anndata.X[ctrl_mask, :],
    index = anndata.obs_names[ctrl_mask],
    columns = anndata.var_names
  ).T

  case_df = pd.DataFrame(
    data = anndata.X[case_mask, :],
    index = anndata.obs_names[case_mask],
    columns = anndata.var_names
  ).T

  signature = deseq2_differential_expression(
    ctrl_df,
    case_df
  ).rename({
    'stat': 'Statistic',
    'pvalue': 'Pval',
    'padj': 'AdjPval',
    'log2FoldChange': 'LogFC',
  }, axis=1)[['Statistic', 'Pval', 'AdjPval', 'LogFC']]
  return signature

def deseq2_from_matrix(file):
  anndata = anndata_from_file(file)
  gene_sig = deseq2(anndata)
  with upsert_file('.tsv', description=f"Gene signature computed by the DESeq2 analysis from the {file.get('description') or 'file'}") as f:
    gene_sig.to_csv(f.file, sep='\t')
  return gene_signature(f)

