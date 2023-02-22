import os
import pandas as pd
from bravado.client import SwaggerClient
from bravado.requests_client import RequestsClient


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


class TumorsOfGene:

    def __init__(self):
        self.rna_df: pd.DataFrame = None
        self.mutations_df: pd.DataFrame = None
        self.patients_df: pd.DataFrame = None


    def build_tumor_table(self) -> pd.DataFrame:
        self.rna_df = self.rna_df[['entrezGeneId','patientId','sampleId']]

        self.mutations_df = self.mutations_df[['entrezGeneId',
                                               'mutationType',
                                               'mutationStatus',
                                               'startPosition',
                                               'endPosition']]

        self.patients_df = self.patients_df[['patientId',
                                              'cancerTypeId']]

        rna_mut_df = self.rna_df.merge(self.mutations_df,
                                 how='inner',
                                 on='entrezGeneId')

        return rna_mut_df.merge(self.patients_df,
                                how='inner',
                                on='patientId')

def main(gene):
    ingestor = PedcbioIngestor()

    tumors_of_gene = TumorsOfGene()

    tumors_of_gene.mutations_df = ingestor.get_mutations()
    tumors_of_gene.rna_df = ingestor.get_molecular_data()
    tumors_of_gene.patients_df = ingestor.get_patients()

    tumors_of_gene_df = tumors_of_gene.build_tumor_table()[['entrezGeneId','patientId','sampleId',
                                     'mutationType','mutationStatus','startPosition',
                                     'endPosition','cancerTypeId']]
    
    tumors_of_gene_df.drop_duplicates(inplace=True)

    return tumors_of_gene_df.to_dict(orient='records')