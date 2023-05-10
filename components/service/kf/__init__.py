import pandas as pd
import requests


def main(ensembl_id: str) -> dict:
    """
    Retrieves gene expression data for the specified Ensembl gene ID from the Open Pediatric Cancer Atlas API.
    
    Args:
        ensembl_id (str): The Ensembl gene ID to retrieve data for.
        
    Returns:
        dict: A dictionary containing the gene expression data as key-value pairs.
            The dictionary is sorted by TPM_sd in descending order.
    """

    api_url = 'https://openpedcan-api.d3b.io/tpm/gene-all-cancer/json'
    query_params = {'ensemblId':ensembl_id,'includeTumorDesc':'primaryOnly'}

    with requests.Session() as session:
        response = session.get(api_url,params=query_params)
        response.raise_for_status()
        

    request_df = pd.read_json(path_or_buf=response.text,orient='records')

    gene_expression_df = request_df[['TPM_mean','TPM_sd','TPM_median',
                                     'Disease','Gene_symbol','Gene_Ensembl_ID',
                                     'Dataset']].copy(True)

    gene_expression_df.sort_values(by='TPM_sd',ascending=False,inplace=True)

    return gene_expression_df.to_dict(orient='records')