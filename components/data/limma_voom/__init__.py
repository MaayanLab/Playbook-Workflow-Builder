from components.data.gene_count_matrix import anndata_from_file
from components.data.metadata_matrix import metadata_from_file
from components.data.gene_signature import gene_signature
from components.core.file import upsert_file
from maayanlab_bioinformatics.dge import characteristic_direction , limma_voom_differential_expression
import pandas as pd

# Function for computing signatures with characteristic direction
def limma_voom(anndata):
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

  signature = limma_voom_differential_expression(
    ctrl_df,
    case_df
  )
  #signature = signature.sort_values("CD-coefficient", ascending=False)
  #signature = signature.rename(columns={'CD-coefficient': f"{grp_ids[0]} vs. {grp_ids[1]}:CD-coefficient"})
  return signature

def limma_voom_from_matrix(anndata):
  anndata = anndata_from_file(anndata)

  # cd
  gene_sig = limma_voom(anndata)

  with upsert_file('.tsv') as f:
    gene_sig.to_csv(f.file, sep='\t')

  return gene_signature(f)

