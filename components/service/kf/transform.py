import pandas as pd

__all__ = ['TumorsOfGene']

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
                                on='patientId')[['entrezGeneId','patientId','sampleId',
                                     'mutationType','mutationStatus','startPosition',
                                     'endPosition','cancerTypeId']]