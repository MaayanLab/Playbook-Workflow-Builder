import gzip
import json
import requests

def load_gmt(url):
    with requests.get(url) as gmt_res:
        gmt_res.raise_for_status()

        if url.endswith(".gz"):
            gmt_text = gzip.decompress(gmt_res.content).decode('utf-8')
        else:
            gmt_text = gmt_res.text

    gmt = {}
    genes = set()
    gene_set_length = 0
    for line in gmt_text.split('\n'):
        if len(line.strip().split('\t'))==1: continue
        term, _, *geneset = line.strip().split('\t')
        if geneset:
            geneset = list(geneset)
            gmt[term] = {"set": geneset}
            genes.update(geneset)
            gene_set_length += len(geneset)

    return gmt

def load_gmts(datasets):
    gmts = []
    for (key,dataset) in datasets:
        dataset["gmt"] = load_gmt(dataset["url"])
        gmts.append({"key":key, "dataset":dataset})

    return gmts
