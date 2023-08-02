from components.data.gene_count_matrix import anndata_from_file
from components.data.metadata_matrix import metadata_from_file
from components.data.gene_signature import gene_signature
from components.core.file import upsert_file
from maayanlab_bioinformatics.dge import characteristic_direction, logfc_differential_expression
import pandas as pd
import scipy.stats

# Function for computing signatures with characteristic direction
def cd_signature(anndata):
  col = anndata.obs.columns[0]
  grp_ids = anndata.obs[col].unique()

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

  signature = pd.concat([
    characteristic_direction(
      ctrl_df,
      case_df,
    ),
    logfc_differential_expression(
      ctrl_df,
      case_df,
    ),
  ], axis=1)
  signature['Pval'] = 1-scipy.stats.norm.sf(abs(scipy.stats.zscore(signature['CD-coefficient'])))*2
  signature.sort_values('Pval', ascending=True, inplace=True)
  signature['AdjPval'] = scipy.stats.false_discovery_control(signature['Pval'], method='bh')
  return signature

def cd_from_matrix(file):
  anndata = anndata_from_file(file)

  # cd
  gene_sig = cd_signature(anndata)

  with upsert_file('.tsv', description=f"Gene signature computed by the Characteristic Direction method of the {file.get('description') or 'file'}") as f:
    gene_sig.to_csv(f.file, sep='\t')

  return gene_signature(f)
