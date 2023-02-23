import pandas as pd
from . ingestor import PedcbioIngestor
from . transform import TumorsOfGene


def main(gene):
    ingestor = PedcbioIngestor()

    tumors_of_gene = TumorsOfGene()

    tumors_of_gene.mutations_df = ingestor.get_mutations()
    tumors_of_gene.rna_df = ingestor.get_molecular_data(gene)
    tumors_of_gene.patients_df = ingestor.get_patients()

    tumors_of_gene_df: pd.DataFrame = tumors_of_gene.build_tumor_table()
    
    tumors_of_gene_df.drop_duplicates(inplace=True)

    return tumors_of_gene_df.to_dict(orient='records')