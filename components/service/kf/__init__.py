import pandas as pd
import requests as rq


def main(ensembl_id):
    api_url = 'https://openpedcan-api.d3b.io/tpm/gene-all-cancer/json'
    query_params = {'ensemblId':ensembl_id,'includeTumorDesc':'primaryOnly'}

    get_request = rq.get(api_url,query_params)
    request_df = pd.read_json(get_request.text,'records')

    gene_expression_df = request_df[['TPM_mean','TPM_sd','TPM_median',
                                     'Disease','Gene_symbol','Gene_Ensembl_ID',
                                     'Dataset']].copy(True)

    gene_expression_df.sort_values(by='TPM_sd',ascending=False,inplace=True)

    return gene_expression_df.to_dict(orient='records')
