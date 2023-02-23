import os
import pandas as pd
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient

__all__ = ['PedcbioIngestor']

class PedcbioIngestor:

    url = 'https://pedcbioportal.kidsfirstdrc.org/api/api-docs'
    token = os.environ['PEDCBIO_TOKEN']

    def __init__(self):
        self.cbio = PedcbioIngestor.get_request_client()

    @classmethod
    def get_request_client(cls) -> SwaggerClient:
        client = RequestsClient()
        client.session.headers = {'Authorization': f'Bearer {PedcbioIngestor.token}'}

        cbioportal = SwaggerClient.from_url(
            PedcbioIngestor.url,
            http_client=client,
            config={"validate_swagger_spec": False}
        )

        return cbioportal

    def get_mutations(self):
        mutations = self.cbio.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
            molecularProfileId="pbta_all_mutations", # {study_id}_mutations gives default mutations profile for study 
            sampleListId="pbta_all_all", # {study_id}_all includes all samples
            projection="SUMMARY" # include gene info
        ).result()


        mdf = pd.DataFrame.from_dict([
            # python magic that combines two dictionaries:
            dict(
                {k:getattr(m,k) for k in dir(m)})
            # create one item in the list for each mutation
            for m in mutations
        ])

        return mdf    


    def get_molecular_data(self, target_gene='7157'):

        rna = self.cbio.Molecular_Data.getAllMolecularDataInMolecularProfileUsingGET(
            entrezGeneId=target_gene, #7157
            sampleListId="pbta_all_all",
            molecularProfileId="pbta_all_rna_seq_v2_mrna",
            projection="SUMMARY"
        ).result()


        rdf = pd.DataFrame.from_dict([
            dict(
                {k: getattr(item,k) for k in dir(item)},
                **{k:getattr(item.gene,k) for k in dir(item.gene)})
            # create one item in the list for each mutation
            for item in rna
        ])

        return rdf


    def get_patients(self):
        patients = self.cbio.Patients.getAllPatientsInStudyUsingGET(
            studyId='pbta_all',
            projection='DETAILED'
        ).result()


        pdf = pd.DataFrame.from_dict([
            dict(
                {k: getattr(patient,k) for k in dir(patient)},
                **{k:getattr(patient.cancerStudy,k) for k in dir(patient.cancerStudy)})
            # create one item in the list for each mutation
            for patient in patients
        ])

        pdf.drop(columns=['cancerStudy'],inplace=True)

        return pdf 