import typing
from components.core.file import File, upsert_file
from components.data.anndata import anndata, AnnDataFile
from components.data.gene_count_matrix import anndata_from_file
from maayanlab_bioinformatics.harmonization import ncbi_genes_lookup
from tempfile import gettempdir
import anndata as ad
import pandas as pd

biotype = typing.Literal["unknown", "tRNA", "rRNA", "snRNA", "scRNA", "snoRNA", "protein-coding", "pseudo", "transposon", "miscRNA", "ncRNA", "biological-region", "other"]

def create_lookup(species:typing.Literal['human','mouse']='human', biotypes:typing.Set[biotype]=None) -> typing.Dict[str, str]:
    organism = 'Mammalia/Mus_musculus' if species=='mouse' else 'Mammalia/Homo_sapiens'
    return ncbi_genes_lookup(
        organism=organism,
        filters=lambda ncbi: ncbi['type_of_gene'].isin(biotypes),
        cache_dir=gettempdir()
    )

def gene_biotype_filter(anndata_file: AnnDataFile, species:typing.Literal['human','mouse']='human', biotypes:typing.List[biotype]=None) -> AnnDataFile:
    if not biotypes:
        raise ValueError("No biotypes selected.")
    lookup = create_lookup(species, set(biotypes))
    input_anndata = anndata_from_file(anndata_file)
    gene_map = input_anndata.var.index.map(lookup)
    masked_anndata = input_anndata[:,(~gene_map.isna())].copy()
    output_anndata = ad.AnnData(
        X = pd.DataFrame(masked_anndata.X).T.groupby(masked_anndata.var_names).sum().T.to_numpy(),
        obs = input_anndata.obs,
        var=masked_anndata.var.groupby(masked_anndata.var_names).first()
    )
    if output_anndata.var.shape[0] == 0:
        raise ValueError(f"No genes remaining after filtering for biotypes: {', '.join(biotypes)}.")
    with upsert_file('.h5ad', description=anndata_file.get('description')) as f:
        output_anndata.write_h5ad(f.file)
    return anndata(f)
