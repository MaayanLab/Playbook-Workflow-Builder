from components.data.gene_count_matrix import anndata_from_file
from components.data.metadata_matrix import metadata_from_file
from components.data.gene_signature import gene_signature
from components.core.file import upsert_file
from maayanlab_bioinformatics.dge import characteristic_direction
import pandas as pd

# Function for computing signatures with characteristic direction
def cd_signature(data_mat, meta_mat):
  col = meta_mat.columns[0]
  grp_ids = meta_mat[col].unique()

  ctrl_ids = meta_mat[meta_mat[col] == grp_ids[0]].index.tolist()
  ctrl_mask = [x in ctrl_ids for x in data_mat.var_names]
  case_ids = meta_mat[meta_mat[col] == grp_ids[1]].index.tolist()
  case_mask = [x in case_ids for x in data_mat.var_names]

  ctrl_df = pd.DataFrame(
    data = data_mat.X[:, ctrl_mask],
    columns = data_mat.var_names[ctrl_mask],
    index = data_mat.obs_names
  )

  case_df = pd.DataFrame(
    data = data_mat.X[:, case_mask],
    columns = data_mat.var_names[case_mask],
    index = data_mat.obs_names
  )

  signature = characteristic_direction(
    ctrl_df,
    case_df
  )
  signature = signature.sort_values("CD-coefficient", ascending=False)
  signature = signature.rename(columns={'CD-coefficient': f"{grp_ids[0]} vs. {grp_ids[1]}:CD-coefficient"})
  return signature

def cd_from_matrix(data, meta):
  data_df = anndata_from_file(data)
  meta_df = metadata_from_file(meta)

  # cd
  gene_sig = cd_signature(data_df, meta_df)

  with upsert_file('.tsv') as f:
    gene_sig.to_csv(f.file, sep='\t')

  return gene_signature(f)
