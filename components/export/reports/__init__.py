from components.data.gene_count_matrix import anndata_from_file
from components.core.file import upsert_file

import html
import io
import numpy as np
import os
import pandas as pd
import re
import requests
import sys
from textwrap import dedent

import asyncio
from openai import AsyncOpenAI
from xml.etree import ElementTree as ET

from adjustText import adjust_text
import base64
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def log(message: str):
    print(message, file=sys.stderr, flush=True)


SYSTEM_PROMPT = """You are a scientific writing assistant generating sections of a bioinformatics re-analysis report.

The re-analysis workflow retrieves RNA-seq expression data and metadata from ARCHS4, performs 
differential expression analysis using limma-voom to create a gene signature, and submits the 
signature to Enrichr and Perturb-Seqr for enrichment analysis.

Output rules that apply to every response:
- Use only plain text prose. No markdown, no bullet points, no numbered lists, no bold, no italics, no headers.
- When an instruction explicitly asks for an outline, use a dash to begin each point. In all other cases, write in paragraph prose only.
- Cite sources inline using this exact format (no other brackets are acceptable): [[[key1, key2]]]. Do not invent citation keys, do not reformat them. Citations you include MUST come only from input publications, results, or explicit instructions.
- If you are explicitly instructed to exclude citations for a section, you must do so.
- Do not speculate or introduce claims beyond what the provided source material explicitly supports.
- When drawing on PMC articles, use only the title, abstract, and introduction sections unless the instruction specifies otherwise. Do not draw from the methods, results, or discussion of the source articles."""


def format_citations(input_text:str):
    def repl(match):
        cites = re.findall(r'\[\[\[([^\]]+)\]\]\]', match.group(0))
        return f"[[[{','.join(cites)}]]]"

    double_pattern = r'(?:(?<!\[)\[\[([^\]]+)\]\](?!\]))+'
    triple_citations = re.sub(double_pattern, repl, input_text)
    
    pattern = r'(?:\[\[\[([^\]]+)\]\]\])+'
    merged_citations = re.sub(pattern, repl, triple_citations)

    return re.sub(r"\[\[\[(.+?)\]\]\]", r"~\\cite{\1}", merged_citations)


async def generate_section(client: AsyncOpenAI, model:str, prompt:str, data:str):
    
    history = [
        {"role":"system","content":SYSTEM_PROMPT},
        {"role":"user","content":prompt},
        {"role": "user", "content": data}
    ]

    response = await client.responses.create(
        model=model,
        input=history
    )

    history.extend(
        (
            {"role": "assistant", "content": response.output_text},
            {
                "role": "user",
                "content": dedent('''
                    Review the paragraph above against the instructions given earlier in this conversation.
                    Check and correct each of the following:
                        1. All citations use the exact format [[[key]]] with no variations. You must use triple-brackets for formatting purposes.
                        2. No citation keys have been added that were not present in the source material or allowed reference lists.
                        3. The output is plain prose only with no markdown, bullet points, or headings.
                        4. No claims are made that go beyond the scope of the provided data or source material.
                    Return only the corrected paragraph. Do not explain your changes.''').strip(),
            },
        )
    )

    response = await client.responses.create(
        model=model,
        input=history

    )

    return response.output_text.encode('ascii', errors='ignore').decode('ascii')
                            
    

async def write_report_abstract(client: AsyncOpenAI, model:str, geo_accession:str, introduction:dict[str,str], methods:str, results:str, discussion:dict[str,str]):
    log("Writing abstract...")
    prompt = dedent('''
    Write a single plain text paragraph abstract of approximately 200 words for a re-analysis report.

    Cover the following points in order, briefly:
    1. The biological or clinical problem the original study addressed.
    2. The motivation for re-analyzing the data from that study.
    3. What the re-analysis workflow did (include tools used but not specific methods)
    4. One or two specific findings from the discussion that best represent the re-analysis results.

    Do not copy sentences verbatim from the provided sections.
    You MUST NOT include citations or reference keys of any kind.''').strip()
    data = dedent(f'''
    GEO Accession:{geo_accession}
    Introduction:{introduction}
    Methods: {methods}
    Results: {results}
    Discussion: {discussion}''').strip()

    abstract = await generate_section(client, model, prompt, data)
    abstract = re.sub(r"\[\[\[(.+?)\]\]\]", "", abstract)
    abstract = re.sub(r"~\\cite\{(?:.+?)\}","", abstract).strip(',')

    log("Abstract complete.")
    return abstract


async def write_report_introduction(client: AsyncOpenAI, model:str, geo_accession:str, pmc_articles:list, labelled_samples:pd.DataFrame):
    log("Writing introduction...")
    pmc_articles = str(pmc_articles)

    problem_prompt = dedent('''
    Write a plain text paragraph introducing the biological or clinical problem  investigated in the GEO study.
    Use only the titles, abstracts, and introductions of the provided PMC articles as your source.
    Describe the problem being studied and why it is relevant. Do not describe the methods used in the original study.
    Do not mention the re-analysis. 
    Include inline citations for any claims you make using references that appear in the articles.''').strip()

    background_prompt = dedent('''
    Write a plain text paragraph providing biological background that contextualizes the problem introduced in the GEO study.
    Use only the titles, abstracts, and introductions of the provided PMC articles as your source.
    Focus on established biology, prior work, or relevant context that motivates the study.
    Expand on the previous problem paragraph, but avoid repeating information or abbreviations that have already been introduced.
    Do not describe any analysis methods from the original study or the re-analysis.
    Do not repeat the specific problem statement - assume that has already been introduced.
    Include inline citations for any claims you make using references that appear in the articles.''').strip()

    signature_prompt = dedent('''
    Create a short name for the signature being re-analyzed by identifying the condition(s)
    that are being used as perturbations in the labelled samples. Return only the name of the signature. For example:
    "metformin signature" or "melanoma cell line signature". Avoid acronyms and be as concise as possible.''').strip()


    introduction_problem, SIGNATURE_NAME= await asyncio.gather(
        generate_section(client, model,problem_prompt, pmc_articles),
        generate_section(client, model, signature_prompt, labelled_samples)
    )
    introduction_background = await generate_section(client, model, background_prompt, f"Articles:{pmc_articles}\nIntroduction Problem:{introduction_problem}")

    introduction_motivation = dedent(f'''
    In order to further investigate any underlying mechanisms or regulatory activity, we performed a re-analysis of samples from {geo_accession} ~\\cite{{{geo_accession}}}
    to create a workflow analyzing the {SIGNATURE_NAME} signature utillizing bioinformatics tools. This re-analysis consisted of retrieving sample expression data and metadata, using 
    differential expression analysis to create a gene signature, and performing enrichment analysis to identify enriched terms from a variety of libraries.''').strip().replace('\n', ' ')
    log("Introduction complete.")

    return {
        "problem":introduction_problem, 
        "background":introduction_background,
        "motivation":introduction_motivation
    }


def write_report_methods(geo_accession:str, labelld_samples:pd.DataFrame):
    return dedent(f'''The workflow starts with selecting {geo_accession} as the search term. GEO studies were identified matching {geo_accession} using 
    ARCHS4 ~\\cite{{ARCHS4}} term search. The GEO study accession was used to fetch the linked publication accession from PMC. Gene expression counts and sample 
    metadata for published samples were obtained from ARCHS4 ~\\cite{{ARCHS4}}. An AnnData file was prepared from the input data and metadata ~\\cite{{AnnData}}. Genes from the 
    anndata matrix were filtered to include protein-coding genes. The samples were then labeled as either control or perturbation to allow for 
    further analysis. The AnnData file was then visualized as a bar plot representing library sizes. Dimensionality reduction of the data was performed 
    using PCA with the normalization set to log-counts-per-million (logCPM). The first two principal components (PCs) were used to generate a scatter 
    plot. The AnnData file was then analyzed using differential expression by Limma-Voom ~\\cite{{limma, voom}} to create a gene signature using the selected conditions. 
    The data in the differential expression table was then visualized as a volcano plot. The up-regulated genes were extracted from the gene signature 
    computed by the Limma-Voom analysis from the file. The gene sets containing significant up and down genes were extracted from the gene signatureand submitted to 
    Enrichr ~\\cite{{Enrichr}}. The gene sets were enriched against the GO Biological Process 2023 ~\\cite{{GO}}, KEGG 2021 Human ~\\cite{{KEGG}}, ChEA 2022 ~\\cite{{ChEA}}, and KOMP2 Mouse Phenotypes 
    2022 ~\\cite{{KOMP2}} libraries to identify statistically significant enriched biological processes, pathways, transcription factors and phenotypes. 
    Significant genes were extracted from the gene signature and submitted to Perturb-Seqr ~\\cite{{Perturb-Seqr}} to identify small molecules and single
    gene perturbations producing gene expression profiles similar or opposite to the signature.''').strip().replace('\n',' ')


def write_report_results(geo_accession:str, labelld_samples:pd.DataFrame, signature:dict[str,str|float], enrichr_results:dict[str,dict[str,pd.DataFrame]], perturbseqr_results:dict[str,pd.DataFrame]):
    def clean_term(term: str):
        return term.split('(')[0].strip()

    def natural_join(items):
        """
        Oxford comma join:
        A
        A and B
        A, B, and C
        """
        if len(items) == 0:
            return ""
        if len(items) == 1:
            return items[0]
        if len(items) == 2:
            return f"{items[0]} and {items[1]}"

        return f"{', '.join(items[:-1])}, and {items[-1]}"

    def summarize_library(library: str, terms, n_terms: int = 2):
        cleaned_terms = []
        seen = set()

        for term in terms:
            cleaned = clean_term(term)
            normalized = cleaned.lower()

            if normalized not in seen:
                seen.add(normalized)
                cleaned_terms.append(cleaned)

            if len(cleaned_terms) >= n_terms:
                break

        joined_terms = natural_join(cleaned_terms)
        lib_lower = library.lower()

        standard_rules = [
            ("go", lambda: f"GO Biological Process terms related to {joined_terms.lower()} ~\\cite{{GO}}"),
            ("kegg", lambda: f"KEGG pathways involving {joined_terms.lower()} ~\\cite{{KEGG}}"),
            ("chea", lambda: f"ChEA transcription factors including {joined_terms} ~\\cite{{ChEA}}"),
            ("komp", lambda: f"KOMP2 mouse phenotypes associated with {joined_terms.lower()} ~\\cite{{KOMP2}}"),
        ]

        perturbseqr_citations = {
            "lincs l1000 xpr": "LINCS",
            "lincs l1000 cp": "LINCS",
            "perturb atlas human": "PerturbAtlas",
            "perturb atlas mouse": "PerturbAtlas",
            "creeds gene": "CREEDS",
            "creeds chem": "CREEDS",
            "rummageo gene": "RummaGEO",
            "rummageo chem": "RummaGEO",
            "replogle et al.": "Replogle",
            "cm4ai": "CM4AI",
            "tahoe-100m": "Tahoe",
            "microarrays cmap": "CMap",
            "nibr drug-seq": "NIBR",
            "sciplex": "SciPlex",
            "deepcover moa": "DeepCoverMoA",
            "ginkgo bioworks": "Ginkgo",
        }

        for pattern, formatter in standard_rules:
            if pattern in lib_lower:
                return formatter()

        if lib_lower in perturbseqr_citations:
            return f"{joined_terms} from {library} ~\\cite{{{perturbseqr_citations[lib_lower]}}}"

        return f"{joined_terms} from {library}"


    def format_enrichr_section(results, n_terms=2):
        library_summaries = [
            summarize_library(lib.rsplit("_", 1)[0], terms_df["term"], n_terms=n_terms)
            for lib, terms_df in results.items()
        ]

        return natural_join(library_summaries)
    

    def format_perturbseqr_section(df, n_terms_per_dataset: int = 2):
        summaries = [
            summarize_library(dataset_name,dataset_df["Perturbation"],n_terms=n_terms_per_dataset)
            for dataset_name, dataset_df in df.groupby("Dataset")
        ]

        return natural_join(summaries)

    up_genes,down_genes = [],[]
    for gene in signature:
        symbol = gene['term']
        score = float(gene['zscore'])
        if score and score > 0:
            up_genes.append(symbol)
        elif score and score < 0:
            down_genes.append(symbol)
    down_genes = down_genes[::-1]

    return dedent(f'''
    Library size distributions across samples were visualized to assess sequencing depth and sample consistency (Figure \\ref{{fig:librarySizes}}).
    PCA of normalized expression profiles revealed the relationship between control and perturbation samples, while also displaying additional study samples not included in the differential expression analysis (Figure \\ref{{fig:PCAScatter}}).
    Differential expression analysis identified {len(up_genes)} significantly up-regulated genes and {len(down_genes)} significantly down-regulated genes after logFC and adjusted p-value filtering (Figure \\ref{{fig:volcanoScatter}}).
    Among the most strongly up-regulated genes were {', '.join(up_genes[:4])}, and {up_genes[4]}, whereas prominent down-regulated genes included {', '.join(down_genes[:4])}, and {down_genes[4]}.
    Functional enrichment analysis of the up-regulated gene set identified associations with {format_enrichr_section(enrichr_results['enrichr_up'])} (Figure \\ref{{fig:upEnrichrBars}}). 
    Analysis of the down-regulated gene set identified enrichment for {format_enrichr_section(enrichr_results['enrichr_down'])} (Figure \\ref{{fig:downEnrichrBars}}).
    Signature search revealed gene perturbations including {format_perturbseqr_section(perturbseqr_results['mimic_gene_signatures'])} (Table \\ref{{table:perturbseqrGeneMimickers}}) and drug perturbations including {format_perturbseqr_section(perturbseqr_results['mimic_drug_signatures'])} (Table \\ref{{table:perturbseqrDrugMimickers}}) as mimickers.
    Top reversers includeded genes perturbations {format_perturbseqr_section(perturbseqr_results['reverse_gene_signatures'])} (Table \\ref{{table:perturbseqrGeneReversers}}) and drugs including {format_perturbseqr_section(perturbseqr_results['reverse_drug_signatures'])} (Table \\ref{{table:perturbseqrDrugReversers}}).
    ''').strip().replace('\n',' ')


async def write_report_discussion(client: AsyncOpenAI, model:str, geo_accession:str, pmc_articles:list, labelled_samples:pd.DataFrame, enrichr_results:dict[str,dict[str,pd.DataFrame]], perturbseqr_results:dict[str,pd.DataFrame]):
    log("Writing discussion...")
    pmc_articles = str(pmc_articles)

    enrichr_prompt = dedent('''
    Write a plain text paragraph analyzing the Enrichr enrichment results provided for the up-regulated and down-regulated gene sets from this re-analysis.
    Do not list the terms to introduce them, assume this has already been done.
    Focus on: terms that appear or are consistent across multiple libraries, and any complementary or contrasting patterns between the up and down gene sets.
    Only discuss terms that are present in the results provided.
    Do not introduce terms or biological processes not present in the data.
    Cite only the libraries whose results you are directly referencing, choosing from: 
    [[[GO]]], [[[KEGG]]], [[[ChEA]]], [[[KOMP2]]].''').strip()

    perturbseqr_prompt = dedent('''
    Write a plain text paragraph analyzing the Perturb-Seqr results provided for the gene signature from this re-analysis.
    Focus on: small molecules or genetic perturbations that appear in each mimicker and reverser result table and what those patterns suggest about the biology of the signature.
    Do not list the terms to introduce them, assume this has already been done.
    Highlight similarities in perturbations within each result table.
    Only discuss entries that are present in the results provided.
    Do not introduce perturbations or mechanisms not present in the data.
    When you reference a mimicker or reverser signature, you MUST cite the associated resource (found in Dataset field), choosing from:
    [[[CMap]]], [[[CM4AI]]], [[[CREEDS]]], [[[DeepCoverMoA]]], [[[Ginkgo]]], [[[LINCS]]], [[[NIBR]]], 
    [[[PerturbAtlas]]], [[[RummaGEO]]], [[[Replogle]]], [[[SciPlex]]], [[[Tahoe]]].''').strip()

    discussion_enrichr, discussion_perturbseqr = await asyncio.gather(
        generate_section(client, model, enrichr_prompt, str(enrichr_results)),
        generate_section(client, model, perturbseqr_prompt, str(perturbseqr_results)),
    )

    conclusion_prompt = dedent('''
    Write a plain text concluding paragraph for a re-analysis report.
    Summarize the most notable findings from the Enrichr and Perturb-Seqr analyses provided.
    Focus on findings that are consistent with or directly relevant to the biology of the samples.
    Where a finding is unexpected relative to the sample context, note it briefly.
    Do not propose mechanisms unless they are directly supported by the enrichment results provided.
    Connect the findings back to the research question described in the Enrichr and Perturb-Seqr paragraphs.
    Do not include citations.''').strip()

    conclusion_data = f'''Samples: {labelled_samples}\nEnrichr Results:{discussion_enrichr}\nPerturb-Seqr Results:{discussion_perturbseqr}'''
    discussion_conclusion = await generate_section(client, model, conclusion_prompt, conclusion_data)

    log("Discussion complete.")

    return {
        "enrichr":discussion_enrichr, 
        "perturbseqr":discussion_perturbseqr,
        "conclusion":discussion_conclusion
    }    


# Data extraction utility functions
def retrieve_pmc_articles(pmc_set: set[str]):
    ids = ",".join(pid.replace("PMC", "") for pid in pmc_set)
    resp = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={"id": ids, "db": "pmc"},
    )
    resp.raise_for_status()
    return resp.text

def _page_range(el, fpage_tag="fpage", lpage_tag="lpage"):
    fpage = el.findtext(fpage_tag, "").strip()
    lpage = el.findtext(lpage_tag, "").strip()
    return f"{fpage}--{lpage}" if fpage else ""

def _names_to_author_list(name_els) -> list[str]:
    authors = []
    for name in name_els:
        surname = name.findtext("surname", "").strip()
        given = name.findtext("given-names", "").strip()
        if surname:
            authors.append(f"{surname}, {given}" if given else surname)
    return authors

def _extract_citation_author_list(cit) -> list[str]:
    person_group = cit.find("person-group[@person-group-type='author']")
    if person_group is not None:
        name_els = person_group.findall("name") + person_group.findall("string-name")
        if name_els:
            return _names_to_author_list(name_els)
    name_els = cit.findall("name") + cit.findall("string-name")
    return _names_to_author_list(name_els)

def _ref_dict(entry_type: str, key: str, fields: dict) -> dict:
    clean = {"id": key, "type": entry_type}
    for k, v in fields.items():
        if isinstance(v, list):
            if v:
                clean[k] = v
        else:
            v = str(v).strip()
            if v and v != "--":
                clean[k] = html.unescape(v)
    return clean

def extract_text_with_citations(elem):
    parts = []
    if elem.text:
        parts.append(elem.text)
    for child in elem:
        if child.tag == "xref" and child.attrib.get("ref-type") == "bibr":
            parts.append(f"[{child.get('rid')}]")
        else:
            parts.append(extract_text_with_citations(child))
        if child.tail:
            parts.append(child.tail)
    return "".join(parts)

def parse_section(sec):
    title_el = sec.find("title")
    title = "".join(title_el.itertext()).strip() if title_el is not None else ""
    paragraphs = [
        extract_text_with_citations(p).strip()
        for p in sec.findall("p")
    ]
    text = "\n\n".join(p for p in paragraphs if p)
    blocks = [{"title": title, "text": text}]
    for subsec in sec.findall("sec"):
        blocks.extend(parse_section(subsec))
    return blocks

def _parse_article_bib(geo_accession: str, front) -> dict:
    doi = front.findtext(".//article-id[@pub-id-type='doi']", "").strip()
    return _ref_dict("article", geo_accession, {
        "title":   front.findtext(".//article-title", "").strip(),
        "authors": _names_to_author_list(front.findall(".//contrib-group//name")),
        "journal": front.findtext(".//journal-id[@journal-id-type='nlm-ta']", ""),
        "year":    front.findtext(".//year", ""),
        "volume":  front.findtext(".//volume", ""),
        "pages":   _page_range(front),
        "doi":     doi,
        "pmid":    front.findtext(".//article-id[@pub-id-type='pmid']", ""),
        **({"url": f"https://doi.org/{doi}"} if doi else {}),
    })

def _parse_reference_bib(ref) -> dict | None:
    ref_id = ref.attrib.get("id", "ref_unknown")
    cit = ref.find("element-citation") or ref.find("mixed-citation")
    if cit is None:
        return None

    pub_type = cit.attrib.get("publication-type", "")
    entry_type = "article" if pub_type == "journal" else "misc"
    doi = cit.findtext(".//pub-id[@pub-id-type='doi']", "").strip()

    return _ref_dict(entry_type, ref_id, {
        "title":   (cit.findtext("article-title") or "").strip(),
        "authors": _extract_citation_author_list(cit),
        "journal": (cit.findtext("source") or "").strip(),
        "year":    cit.findtext("year", ""),
        "volume":  cit.findtext("volume", ""),
        "pages":   _page_range(cit),
        "doi":     doi,
        "pmid":    cit.findtext(".//pub-id[@pub-id-type='pmid']", ""),
        **({"url": f"https://doi.org/{doi}"} if doi else {}),
    })

def parse_pmc_xml(geo_accession: str, pmc_set: set[str]):
    root = ET.fromstring(retrieve_pmc_articles(pmc_set))

    articles = []
    for article in root:
        front = article.find(".//front")
        abstract_el = article.find(".//abstract")

        abstract = " ".join(
            "".join(p.itertext()).strip()
            for p in abstract_el.findall(".//p")
        ) if abstract_el is not None else ""

        body = article.find(".//body")
        sections = []
        if body is not None:
            for child in body:
                if child.tag == "p":
                    if text := extract_text_with_citations(child).strip():
                        sections.append({"title": "", "text": text})
                elif child.tag == "sec":
                    sections.extend(parse_section(child))

        references = [_parse_article_bib(geo_accession, front)]
        for ref in article.findall(".//ref"):
            if entry := _parse_reference_bib(ref):
                references.append(entry)

        articles.append({
            "title":      article.findtext(".//article-title", "").strip(),
            "abstract":   abstract,
            "body":       sections,
            "references": references,
        })

    return articles


def extract_labelled_samples(anndata_file):
    meta_df = anndata_from_file(anndata_file).obs
    meta_df['Type: Control or Perturbation'] = meta_df['Type: Control or Perturbation'].astype(str).replace("nan", pd.NA).copy()
    return meta_df.dropna()


def extract_enrichr_results(up_enrichment, down_enrichment):
    up_id = up_enrichment.pop("enrichr_up_id")["shortId"]
    up_result = {
        library:pd.DataFrame(results["scored"]).head(10) for library,results in up_enrichment.items()
    }
    down_id = down_enrichment.pop("enrichr_down_id")["shortId"]
    down_result = {
        library:pd.DataFrame(results["scored"]).head(10) for library,results in down_enrichment.items()
    }
    return {"enrichrUp":up_id,"enrichrDown":down_id},{"enrichr_up":up_result, "enrichr_down":down_result}


def extract_perturbseqr_results(perturbseqr_ids:dict[str,str], n:int=20):
    def clean_df(df:pd.DataFrame, dir:str, pert:str):
        if dir=="mimic":
            or_col="oddsRatioMimic"
            pval_col = "pvalueMimic"
            apval_col = "adjPvalueMimic"
            drop_cols = ["signatureCount", "nReverseOverlap", "oddsRatioReverse", "pvalueReverse", "adjPvalueReverse"]
        else:
            or_col="oddsRatioReverse"
            pval_col = "pvalueReverse"
            apval_col = "adjPvalueReverse"
            drop_cols = ["signatureCount", "nReverseOverlap", "oddsRatioMimic", "pvalueMimic", "adjPvalueMimic"]
        colnames = ["Dataset", "Perturbation", "Perturbation ID", "Cell Line", "Timepoint", "Concentration", "MoA", "FDA Approved", "Gene Set Size Up", "Gene Set Size Down", "n Overlap", "Odds Ratio", "p-value", "Adjusted p-value"]
        df = df.head(n).drop(columns=drop_cols).rename_axis("Rank", inplace=False)
        df[or_col] = df[or_col].astype(float).map(lambda x: np.format_float_positional(x, precision=5))
        df[pval_col] = df[pval_col].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
        df[apval_col] = df[apval_col].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
        df.index = df.index.astype(int) + 1
        df.columns=colnames
        if pert=='gene':
            df = df.drop(columns=['Concentration','MoA'])
        return df

    up_id = perturbseqr_ids['up_id']
    down_id = perturbseqr_ids['down_id']
    gene_libraries = 'LINCS L1000 XPR,Perturb Atlas Human,Perturb Atlas Mouse,CREEDS Gene,RummaGEO Gene,Replogle et al.,CM4AI'
    drug_libraries = 'LINCS L1000 CP,Tahoe-100M,Microarrays CMap,NIBR DRUG-seq,SciPlex,DeepCover MoA,CREEDS Chem,RummaGEO Chem,Ginkgo Bioworks'
    perturbseqr_url = 'https://perturbseqr.maayanlab.cloud/enrichpair/download'
    s = requests.Session()
    s.headers.update({'Accept': 'text/tab-separated-values'})
    
    mimic_gene_req = s.get(url=perturbseqr_url, params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_mimic","maxTotal":n, "libraries":gene_libraries})
    mimic_gene_req.raise_for_status()
    mimic_gene_df = clean_df(pd.read_csv(io.BytesIO(mimic_gene_req.content), sep='\t'),'mimic', 'gene')

    mimic_drug_req = s.get(url=perturbseqr_url, params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_mimic","maxTotal":n, "libraries":drug_libraries})
    mimic_drug_req.raise_for_status()
    mimic_drug_df = clean_df(pd.read_csv(io.BytesIO(mimic_drug_req.content), sep='\t'),'mimic', 'drug')

    reverse_gene_req = s.get(url=perturbseqr_url, params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_reverse","maxTotal":n, "libraries":gene_libraries})
    reverse_gene_req.raise_for_status()
    reverse_gene_df = clean_df(pd.read_csv(io.BytesIO(reverse_gene_req.content), sep='\t'),'reverse', 'gene')

    reverse_drug_req = s.get(url=perturbseqr_url, params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_reverse","maxTotal":n, "libraries":drug_libraries})
    reverse_drug_req.raise_for_status()
    reverse_drug_df = clean_df(pd.read_csv(io.BytesIO(reverse_drug_req.content), sep='\t'),'reverse', 'drug')

    return {"mimic_gene_signatures":mimic_gene_df, "mimic_drug_signatures":mimic_drug_df, "reverse_gene_signatures":reverse_gene_df, "reverse_drug_signatures":reverse_drug_df}


# Figure utility functions
def truncate(text:str, length:int=78):
    if len(text) >= length:
        text = f'{text[:length-3]}...'
    return text


def decode_bdata(bdata, dtype=np.float64):
    raw = base64.b64decode(bdata)
    return np.frombuffer(raw, dtype=dtype)


def extract_gene_names(hovertext):
    genes = []
    for h in hovertext:
        m = re.search(r"<b>(.*?)</b>", h)
        genes.append(m[1] if m else "NA")
    return np.array(genes)


def extract_pvals(hovertext):
    pvals = []
    for h in hovertext:
        if m := re.search(
            r"p = ([0-9]+(?:\.[0-9]+)?(?:e[+-]?[0-9]+)?)", h, re.IGNORECASE
        ):
            pvals.append(float(m[1]))
        else:
            raise ValueError(f"Could not parse p-value from: {h}")
    return np.array(pvals)


def make_library_size_barplot(library_size_plotly):
    plotly_data = library_size_plotly['data'][0]

    x = decode_bdata(plotly_data['x']['bdata'], dtype=np.int32)
    y = plotly_data['y']

    x = x[::-1]
    y = y[::-1]

    fig, ax = plt.subplots(figsize=(12, 6), dpi=300)

    bars = ax.barh(
        y,
        x,
        height=0.7,
        edgecolor='none'
    )

    ax.set_xlabel("Library Size", fontsize=18)

    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.xaxis.get_offset_text().set_fontsize(14)

    ax.grid(axis='x', alpha=0.25)
    ax.grid(axis='y', visible=False)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()
    with upsert_file('.pdf') as f, upsert_file('.png') as f1:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f1.file, format='png', dpi=300, bbox_inches='tight')
    return f,f1


def make_pca_scatter(pca_scatter_plotly):
    plotly_data = pca_scatter_plotly['data'][0]
    plotly_scene = pca_scatter_plotly['layout']

    x = decode_bdata(plotly_data['x']['bdata'])
    y = decode_bdata(plotly_data['y']['bdata'])
    colors = np.array(plotly_data['marker']['color'])

    color_map = {
        'blue': '#1f77b4',
        'red': '#d62728',
        'grey': 'lightgray'
    }
    mpl_colors = np.array([color_map.get(c, 'black') for c in colors])


    fig = plt.figure(figsize=(12, 12), dpi=300)
    ax = fig.add_subplot(111)

    ax.scatter(
        x, y,
        c=mpl_colors,
        s=120,
        alpha=0.85,
    )

    ax.set_box_aspect(1)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Control',
            markerfacecolor='#1f77b4', markersize=14),
        Line2D([0], [0], marker='o', color='w', label='Perturbation',
            markerfacecolor='#d62728', markersize=14),
        Line2D([0], [0], marker='o', color='w', label='Other',
            markerfacecolor='lightgray', markersize=14)
    ]

    ax.set_xlabel(plotly_scene['xaxis']['title']['text'], fontsize=24, labelpad=10)
    ax.set_ylabel(plotly_scene['yaxis']['title']['text'], fontsize=24, labelpad=10)
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 0.5), loc='upper left', fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    with upsert_file('.pdf') as f, upsert_file('.png') as f1:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f1.file, format='png', dpi=300, bbox_inches='tight')
    return f,f1


def make_volcano_scatter(volcano_scatter_plotly):
    plotly_data = volcano_scatter_plotly['data'][0]

    x = decode_bdata(plotly_data['x']['bdata'], dtype=np.float64)
    colors = plotly_data['marker']['color']
    hovertext = plotly_data['hovertext']

    genes = extract_gene_names(hovertext)
    pvals = extract_pvals(hovertext)
    y = -np.log10(pvals)

    mpl_colors = []
    for c in colors:
        if c == "blue":
            mpl_colors.append("blue")
        elif c == "red":
            mpl_colors.append("red")
        else:
            mpl_colors.append("lightgray")
    mpl_colors = np.array(mpl_colors)

    p_thresh = 0.05
    sig = (pvals < p_thresh) & (np.abs(x)>1.5)

    dist = np.sqrt(x**2 + y**2)

    up = (x > 0) & sig
    down = (x < 0) & sig

    top_up_idx = np.argsort(dist[up])[-10:]
    top_down_idx = np.argsort(dist[down])[-10:]

    up_indices = np.where(up)[0][top_up_idx]
    down_indices = np.where(down)[0][top_down_idx]

    label_indices = np.concatenate([up_indices, down_indices])

    nonsig_idx = np.where(~sig)[0]
    keep_fraction = 0.5
    sampled_nonsig_idx = np.random.choice(nonsig_idx,size=int(len(nonsig_idx) * keep_fraction),replace=False)
    keep_idx = np.concatenate([np.where(sig)[0],sampled_nonsig_idx])

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

    ax.scatter(x[keep_idx], y[keep_idx], c=mpl_colors[keep_idx], s=20, alpha=0.7, linewidths=0)

    ax.axvline(-1.5, linestyle="--", linewidth=1, alpha=0.4)
    ax.axvline(1.5, linestyle="--", linewidth=1, alpha=0.4)
    ax.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1, alpha=0.4)

    xpad = (x.max() - x.min()) * 0.05
    ypad = (y.max() - y.min()) * 0.08
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.set_ylim(y.min() - ypad, y.max() + ypad)

    texts = [
        ax.text(
            x[i],
            y[i],
            genes[i],
            fontsize=13,
            fontweight='bold',
            ha='center',
            va='center',
        )
        for i in label_indices
    ]
    # Repel labels
    adjust_text(
        texts,
        x=x[label_indices],
        y=y[label_indices],
        ax=ax,
        expand_text=(1.2, 1.4),
        expand_points=(1.2, 1.4),
        force_text=(0.5, 1.0),
        force_points=(0.2, 0.5),
        arrowprops=dict(arrowstyle="-", lw=0.5, color="#585858"),
        only_move={'points':'y', 'text':'xy'},
        lim=300
    )

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Upregulated',
            markerfacecolor='red', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Downregulated',
            markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='Not significant',
            markerfacecolor='lightgray', markersize=8)
    ]

    ax.legend(handles=legend_elements, loc='lower right', frameon=False)

    ax.set_xlabel("log2FC", fontsize=12)
    ax.set_ylabel("-log10(P)", fontsize=12)

    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    with upsert_file('.pdf') as f, upsert_file('.png') as f1:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f1.file, format='png', dpi=300, bbox_inches='tight')
    return f,f1


def make_enrichr_barplot(scored_enrichr_libraries: dict[str,pd.DataFrame]):
    fig, axes = plt.subplots(2, 2, figsize=(18,10), dpi=300)
    dfs = scored_enrichr_libraries.values()
    titles = ["GO Biological Process", "KEGG Pathways", "ChEA Transcription Factors", "KOMP2 Mouse Phenotypes"]
    columns = ["Scored Pathway or Biological Processes", "Scored Pathway or Biological Processes", "Scored Genes", "Scored Phenotypes"]
    colors = ["skyblue","lightgreen","salmon","plum"]
    panel_labels = ["A", "B", "C", "D"]

    for ax, i, df, column, title, color, label in zip(axes.flatten(), range(4), dfs, columns, titles, colors, panel_labels):
        top = df.copy()[::-1]
        top["term"] = top["term"].map(truncate)
        cap = top['zscore'].min()*10
        top['zscore'] = top['zscore'].clip(upper=cap)

        bars = ax.barh(range(len(top)), top['zscore'], color=color)

        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.set_title(title, fontsize=16)
        ax.set_xlabel("ZScore", fontsize=12)

        for i, (term, z) in enumerate(zip(top['term'], top['zscore'])):
            ax.text(
                z * 0.02,
                i,
                term,
                va='center',
                ha='left',
                fontsize=15,
                color='black',
                clip_on=True
            )

        ax.text(-0.02, 1.05, label, transform=ax.transAxes,fontsize=16, fontweight='bold', va='top', ha='right')

    plt.tight_layout()
    with upsert_file('.pdf') as f, upsert_file('.png') as f1:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f1.file, format='png', dpi=300, bbox_inches='tight')
    return f,f1


def make_figures(plots, enrichr_results:dict[str,dict[str,pd.DataFrame]], supplement:dict[str,str]):
    library_pdf, librar_png = make_library_size_barplot(plots["library_sizes_plot"])
    pca_pdf, pca_png = make_pca_scatter(plots["pca_plot"])
    volcano_pdf, volcano_png = make_volcano_scatter(plots["volcano_plot"])
    enrichr_up_pdf, enrichr_up_png = make_enrichr_barplot(enrichr_results["enrichr_up"])
    enrichr_down_pdf, enrichr_down_png = make_enrichr_barplot(enrichr_results["enrichr_down"])
    enrichr_up_id = supplement["enrichrUp"]
    enrichr_down_id = supplement["enrichrDown"]
    return {
        "librarySizes": {
            "pdf":library_pdf,
            "png":librar_png,
            "caption":"Library sizes for each sample in the dataset, shown as total mapped read counts per sample. Samples are labeled by their GEO accession identifier (GSM). Consistent library sizes across samples indicate uniform sequencing depth suitable for differential expression analysis."
        },
        "PCAScatter": {
            "pdf":pca_pdf,
            "png":pca_png,
            "caption":"Principal component analysis (PCA) of normalized gene expression profiles across all samples. Each point represents one sample, colored by experimental condition: control (blue), perturbation (red), and additional study samples not included in the differential expression analysis (gray). Axes indicate the percentage of total variance explained by each principal component. Normalization was performed using log-counts-per-million (logCPM)."
        },
        "volcanoScatter": {
            "pdf":volcano_pdf,
            "png":volcano_png,
            "caption":"Volcano plot of differential expression results comparing perturbation to control samples. Each point represents one protein-coding gene; the x-axis shows the log2 fold-change and the y-axis shows statistical significance as -log10(adjusted p-value). Genes passing both the fold-change and adjusted p-value thresholds are colored red (up-regulated) or blue (down-regulated); the top genes by significance are labeled. Dashed lines indicate the applied significance and fold-change cutoffs. Insignificant points are randomly downsampled by a factor of 0.5.",
        },
        "upEnrichrBars": {
            "pdf":enrichr_up_pdf,
            "png":enrichr_up_png,
            "caption":f"Enrichment analysis of the up-regulated gene signature using Enrichr~\\cite{{Enrichr}}. Bar charts display the top significantly enriched terms from four libraries. A. GO Biological Process 2023~\\cite{{GO}}, B. KEGG 2021 Human~\\cite{{KEGG}}, C. ChEA 2022~\\cite{{ChEA}}, and D. KOMP2 Mouse Phenotypes 2022~\\cite{{KOMP2}}. Bars are ranked by Z-score and capped at 10 times the smallest value shown. Color indicates the source library. The full enrichment results are available to view at \\href{{https://maayanlab.cloud/enrichr/enrich?dataset={enrichr_up_id}}}{{Enrichr}}.",
        },
        "downEnrichrBars": {
            "pdf":enrichr_down_pdf,
            "png":enrichr_down_png,
            "caption":f"Enrichment analysis of the down-regulated gene signature using Enrichr~\\cite{{Enrichr}}. Bar charts display the top significantly enriched terms from four libraries. A. GO Biological Process 2023~\\cite{{GO}}, B. KEGG 2021 Human~\\cite{{KEGG}}, C. ChEA 2022~\\cite{{ChEA}}, and D. KOMP2 Mouse Phenotypes 2022~\\cite{{KOMP2}}. Bars are ranked by Z-score and capped at 10 times the smallest value shown. Color indicates the source library. The full enrichment results are available to view at \\href{{https://maayanlab.cloud/enrichr/enrich?dataset={enrichr_down_id}}}{{Enrichr}}."
        }
    }


# Table utility functions
def make_table(table_data:pd.DataFrame):
    with upsert_file('.tsv') as f:
        table_data.to_csv(f.file, sep='\t')
    return f


def make_tables(perturbseqr_results, supplement:dict[str,str]):
    gene_libraries = 'LINCS L1000 XPR,Perturb Atlas Human,Perturb Atlas Mouse,CREEDS Gene,RummaGEO Gene,Replogle et al.,CM4AI'
    drug_libraries = 'LINCS L1000 CP,Tahoe-100M,Microarrays CMap,NIBR DRUG-seq,SciPlex,DeepCover MoA,CREEDS Chem,RummaGEO Chem,Ginkgo Bioworks'
    perturbseqr_url = f'https://perturbseqr.maayanlab.cloud/enrichpair?dataset={supplement["perturbseqrUpGenes"]}&dataset={supplement["perturbseqrDownGenes"]}'
    return {
        "perturbseqrGeneMimickers": {
            "file":make_table(perturbseqr_results["mimic_gene_signatures"]),
            "caption":f"Single gene perturbations whose transcriptional profiles most closely resemble the query signature, as identified by Perturb-seqr~\\cite{{Perturb-Seqr}}. Each row lists the source dataset, perturbed gene, cell line, timepoint, and statistical measures including gene set overlap size, odds ratio, and adjusted p-value. The full signature search results are available to view at \\href{{{perturbseqr_url}&view=table&dir=up&sort=pvalue_mimic&libraries={gene_libraries}}}{{Perturb-Seqr}}."
        },
        "perturbseqrDrugMimickers": {
            "file":make_table(perturbseqr_results["mimic_drug_signatures"]),
            "caption":f"Small molecule perturbations whose transcriptional profiles most closely resemble the query signature, as identified by Perturb-seqr~\\cite{{Perturb-Seqr}}. Each row lists the source dataset, compound, cell line, timepoint, concentration, mechanism of action (MoA), FDA approval status, and statistical measures including gene set overlap size, odds ratio, and adjusted p-value. The full signature search results are available to view at \\href{{{perturbseqr_url}&view=table&dir=up&sort=pvalue_mimic&libraries={drug_libraries}}}{{Perturb-Seqr}}."
        },
        "perturbseqrGeneReversers": {
            "file":make_table(perturbseqr_results["reverse_gene_signatures"]),
            "caption":f"Single gene perturbations whose transcriptional profiles are most opposite to the query signature, as identified by Perturb-seqr~\\cite{{Perturb-Seqr}}. Each row lists the source dataset, perturbed gene, cell line, timepoint, and statistical measures including gene set overlap size, odds ratio, and adjusted p-value. The full signature search results are available to view at \\href{{{perturbseqr_url}&view=table&dir=down&sort=pvalue_reverse&libraries={gene_libraries}}}{{Perturb-Seqr}}."
        },
        "perturbseqrDrugReversers": {
            "file":make_table(perturbseqr_results["reverse_drug_signatures"]),
            "caption":f"Small molecule perturbations whose transcriptional profiles are most opposite to the query signature, as identified by Perturb-seqr~\\cite{{Perturb-Seqr}}. Each row lists the source dataset, compound, cell line, timepoint, concentration, mechanism of action (MoA), FDA approval status, and statistical measures including gene set overlap size, odds ratio, and adjusted p-value. The full signature search results are available to view at \\href{{{perturbseqr_url}&view=table&dir=down&sort=pvalue_reverse&libraries={drug_libraries}}}{{Perturb-Seqr}}."
        },
    }


def make_references(pmc_articles):
    refs = [
        {"id":"AnnData","type":"article","title":"anndata: Access and store annotated data matrices","authors":["Virshup, Isaac","Rybakov, Sergei","Theis, Fabian J.","Angerer, Philipp","Wolf, F. Alexander"],"year":"2024","journal":"Journal of Open Source Software","volume":"9","pages":"4371","doi":"10.1101/2021.12.16.473007","url":"https://doi.org/10.1101/2021.12.16.473007"},
        {"id":"ARCHS4","type":"article","title":"Massive mining of publicly available RNA-seq data from human and mouse","authors":["Lachmann, Alexander","Torre, Denis","Keenan, Alexandra B.","Jagodnik, Kathleen M.","Lee, Hoyjin J.","Wang, Lily","Silverstein, Moshe C.","Ma'ayan, Avi"],"year":"2018","journal":"Nature Communications","volume":"9","pages":"1366","doi":"10.1038/s41467-018-03751-6","url":"https://doi.org/10.1038/s41467-018-03751-6"},
        {"id":"voom","type":"article","title":"voom: precision weights unlock linear model analysis tools for RNA-seq read counts","authors":["Law, Charity W","Chen, Yunshun","Shi, Wei","Smyth, Gordon K"], "year":"2014","journal":"Genome Biology","volume":"15","pages":"R29","doi":"10.1186/gb-2014-15-2-r29","url":"https://doi.org/10.1186/gb-2014-15-2-r29"},
        {"id":"limma","type":"article","title":"limma powers differential expression analyses for RNA-sequencing and microarray studies","authors":["Ritchie, Matthew E.","Phipson, Belinda","Wu, Di","Hu, Yifang","Law, Charity W.","Shi, Wei","Smyth, Gordon K."], "year":"2015","journal":"Nucleic Acids Research","volume":"43","pages":"e47","doi":"10.1093/nar/gkv007","url":"https://doi.org/10.1093/nar/gkv007"},
        {"id":"Enrichr","type":"article","title":"Gene Set Knowledge Discovery with Enrichr","authors":["Xie, Zhuorui","Bailey, Allison","Kuleshov, Maxim V.","Clarke, Daniel J. B.","Evangelista, John E.","Jenkins, Sherry L.","Lachmann, Alexander","Wojciechowicz, Megan L.","Kropiwnicki, Eryk","Jagodnik, Kathleen M.","Jeon, Minji","Ma'ayan, Avi"],"year":"2021","journal":"Current Protocols","volume":"1","pages":"e90","doi":"10.1002/cpz1.90","url":"https://doi.org/10.1002/cpz1.90"},
        {"id":"GO","type":"article","title":"Gene Ontology: tool for the unification of biology","authors":["Ashburner, Michael","Ball, Catherine A.","Blake, Judith A.","Botstein, David","Butler, Heather","Cherry, J. Michael","Davis, Allan P.","Dolinski, Kara","Dwight, Selina S.","Eppig, Janan T.","Harris, Midori A.","Hill, David P.","Issel-Tarver, Laurie","Kasarskis, Andrew","Lewis, Suzanna","Matese, John C.","Richardson, Joel E.","Ringwald, Martin","Rubin, Gerald M.","Sherlock, Gavin"],"year":"2000","journal":"Nature Genetics","volume":"25","pages":"25--29","doi":"10.1038/75556","url":"https://doi.org/10.1038/75556"},
        {"id":"KEGG","type":"article","title":"KEGG for taxonomy-based analysis of pathways and genomes","authors":["Kanehisa, Minoru","Furumichi, Miho","Sato, Yoko","Kawashima, Masayuki","Ishiguro-Watanabe, Mari"],"year":"2022","journal":"Nucleic Acids Research","volume":"51","pages":"D587--D592","doi":"10.1093/nar/gkac963","url":"https://doi.org/10.1093/nar/gkac963"},
        {"id":"ChEA","type":"article","title":"ChEA3: transcription factor enrichment analysis by orthogonal omics integration","authors":["Keenan, Alexandra B","Torre, Denis","Lachmann, Alexander","Leong, Ariel K","Wojciechowicz, Megan L","Utti, Vivian","Jagodnik, Kathleen M","Kropiwnicki, Eryk","Wang, Zichen","Ma'ayan, Avi"],"year":"2019","journal":"Nucleic Acids Research","volume":"47","pages":"W212--W224","doi":"10.1093/nar/gkz446","url":"https://doi.org/10.1093/nar/gkz446"},
        {"id":"KOMP2","type":"article","title":"The International Mouse Phenotyping Consortium: comprehensive knockout phenotyping underpinning the study of human disease","authors":["Groza, Tudor","Gomez, Federico Lopez","Mashhadi, Hamed Haseli","Muñoz-Fuentes, Violeta","Gunes, Osman","Wilson, Robert","Cacheiro, Pilar","Frost, Anthony","Keskivali-Bond, Piia","Vardal, Bora","McCoy, Aaron","Cheng, Tsz Kwan","Santos, Luis","Wells, Sara","Smedley, Damian","Mallon, Ann-Marie","Parkinson, Helen"],"year":"2022","journal":"Nucleic Acids Research","volume":"51","pages":"D1038--D1045","doi":"10.1093/nar/gkac972","url":"https://doi.org/10.1093/nar/gkac972"},
        {"id":"CMap","type":"article","title":"The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease","authors":["Lamb, Justin","Crawford, Emily D.","Peck, David","Modell, Joshua W.","Blat, Irene C.","Wrobel, Matthew J.","Lerner, Jim","Brunet, Jean-Philippe","Subramanian, Aravind","Ross, Kenneth N.","Reich, Michael","Hieronymus, Haley","Wei, Guo","Armstrong, Scott A.","Haggarty, Stephen J.","Clemons, Paul A.","Wei, Ru","Carr, Steven A.","Lander, Eric S.","Golub, Todd R."],"year":"2006","journal":"Science","volume":"313","pages":"1929--1935","doi":"10.1126/science.1132939","url":"https://doi.org/10.1126/science.1132939"},
        {"id":"CM4AI","type":"article","title":"A Perturbation Cell Atlas of Human Induced Pluripotent Stem Cells","authors":["Nourreddine, Sami","Doctor, Yesh","Dailamy, Amir","Forget, Antoine","Lee, Yi-Hung","Chinn, Becky","Khaliq, Hammza","Polacco, Benjamin","Muralidharan, Monita","Pan, Emily","Zhang, Yifan","Sigaeva, Alina","Hansen, Jan Niklas","Gao, Jiahao","Parker, Jillian A.","Obernier, Kirsten","Clark, Timothy","Chen, Jake Y.","Metallo, Christian","Lundberg, Emma","Ideker, Trey","Krogan, Nevan","Mali, Prashant"],"year":"2024","journal":"bioRxiv","doi":"10.1101/2024.11.03.621734","url":"https://doi.org/10.1101/2024.11.03.621734"},
        {"id":"CREEDS","type":"article","title":"Extraction and analysis of signatures from the Gene Expression Omnibus by the crowd","authors":["Wang, Zichen","Monteiro, Caroline D.","Jagodnik, Kathleen M.","Fernandez, Nicolas F.","Gundersen, Gregory W.","Rouillard, Andrew D.","Jenkins, Sherry L.","Feldmann, Axel S.","Hu, Kevin S.","McDermott, Michael G.","Duan, Qiaonan","Clark, Neil R.","Jones, Matthew R.","Kou, Yan","Goff, Troy","Woodland, Holly","Amaral, Fabio M R.","Szeto, Gregory L.","Fuchs, Oliver","Schüssler-Fiorenza Rose, Sophia M.","Sharma, Shvetank","Schwartz, Uwe","Bausela, Xabier Bengoetxea","Szymkiewicz, Maciej","Maroulis, Vasileios","Salykin, Anton","Barra, Carolina M.","Kruth, Candice D.","Bongio, Nicholas J.","Mathur, Vaibhav","Todoric, Radmila D","Rubin, Udi E.","Malatras, Apostolos","Fulp, Carl T.","Galindo, John A.","Motiejunaite, Ruta","Jüschke, Christoph","Dishuck, Philip C.","Lahl, Katharina","Jafari, Mohieddin","Aibar, Sara","Zaravinos, Apostolos","Steenhuizen, Linda H.","Allison, Lindsey R.","Gamallo, Pablo","de Andres Segura, Fernando","Dae Devlin, Tyler","Pérez-García, Vicente","Ma'ayan, Avi"],"year":"2016","journal":"Nature Communications","volume":"7","pages":"12846","doi":"10.1038/ncomms12846","url":"https://doi.org/10.1038/ncomms12846"},
        {"id":"DeepCoverMOA","type":"article","title":"A proteome-wide atlas of drug mechanism of action","authors":["Mitchell, Dylan C.","Kuljanin, Miljan","Li, Jiaming","Van Vranken, Jonathan G.","Bulloch, Nathan","Schweppe, Devin K.","Huttlin, Edward L.","Gygi, Steven P."],"year":"2023","journal":"Nature Biotechnology","volume":"41","pages":"845--857","doi":"10.1038/s41587-022-01539-0","url":"https://doi.org/10.1038/s41587-022-01539-0"},
        {"id":"Ginkgo","type":"online","title":"Ginkgo Bioworks","url":"https://www.ginkgo.bio/"},
        {"id":"LINCS","type":"article","title":"SigCom LINCS: data and metadata search engine for a million gene expression signatures","authors":["Evangelista, John Erol","Clarke, Daniel J B","Xie, Zhuorui","Lachmann, Alexander","Jeon, Minji","Chen, Kerwin","Jagodnik, Kathleen M","Jenkins, Sherry L","Kuleshov, Maxim V","Wojciechowicz, Megan L","Schürer, Stephan C","Medvedovic, Mario","Ma'ayan, Avi"],"year":"2022","journal":"Nucleic Acids Research","volume":"50","pages":"W697--W709","doi":"10.1093/nar/gkac328","url":"https://doi.org/10.1093/nar/gkac328"},
        {"id":"NIBR","type":"article","title":"DRUG-seq Provides Unbiased Biological Activity Readouts for Neuroscience Drug Discovery","authors":["Li, Jingyao","Ho, Daniel J.","Henault, Martin","Yang, Chian","Neri, Marilisa","Ge, Robin","Renner, Steffen","Mansur, Leandra","Lindeman, Alicia","Kelly, Brian","Tumkaya, Tayfun","Ke, Xiaoling","Soler-Llavina, Gilberto","Shanker, Gopi","Russ, Carsten","Hild, Marc","Gubser Keller, Caroline","Jenkins, Jeremy L.","Worringer, Kathleen A.","Sigoillot, Frederic D.","Ihry, Robert J."],"year":"2022","journal":"ACS Chemical Biology","volume":"17","pages":"1401--1414","doi":"10.1021/acschembio.1c00920","url":"https://doi.org/10.1021/acschembio.1c00920"},
        {"id":"PerturbAtlas","type":"article","title":"PerturbAtlas: a comprehensive atlas of public genetic perturbation bulk RNA-seq datasets","authors":["Zhang, Yiming","Zhang, Ting","Yang, Gaoxia","Pan, Zhenzhong","Tang, Min","Wen, Yue","He, Ping","Wang, Yuan","Zhou, Ran"],"year":"2024","journal":"Nucleic Acids Research","volume":"53","pages":"D1112--D1119","doi":"10.1093/nar/gkae851","url":"https://doi.org/10.1093/nar/gkae851"},
        {"id":"Perturb-Seqr","type":"online","title":"Perturb-Seqr","url":"https://perturbseqr.maayanlab.cloud"},
        {"id":"Replogle","type":"article","title":"Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq","authors":["Replogle, Joseph M.","Saunders, Reuben A.","Pogson, Angela N.","Hussmann, Jeffrey A.","Lenail, Alexander","Guna, Alina","Mascibroda, Lauren","Wagner, Eric J.","Adelman, Karen","Lithwick-Yanai, Gila","Iremadze, Nika","Oberstrass, Florian","Lipson, Doron","Bonnar, Jessica L.","Jost, Marco","Norman, Thomas M.","Weissman, Jonathan S."],"year":"2022","journal":"Cell","volume":"185","pages":"2559--2575.e28","doi":"10.1016/j.cell.2022.05.013","url":"https://doi.org/10.1016/j.cell.2022.05.013"},
        {"id":"RummaGEO","type":"article","title":"RummaGEO: Automatic mining of human and mouse gene sets from GEO","authors":["Marino, Giacomo B.","Clarke, Daniel J.B.","Lachmann, Alexander","Deng, Eden Z.","Ma'ayan, Avi"],"year":"2024","journal":"Patterns","volume":"5","pages":"101072","doi":"10.1016/j.patter.2024.101072","url":"https://doi.org/10.1016/j.patter.2024.101072"},
        {"id":"SciPlex","type":"article","title":"Massively multiplex chemical transcriptomics at single-cell resolution","authors":["Srivatsan, Sanjay R.","McFaline-Figueroa, José L.","Ramani, Vijay","Saunders, Lauren","Cao, Junyue","Packer, Jonathan","Pliner, Hannah A.","Jackson, Dana L.","Daza, Riza M.","Christiansen, Lena","Zhang, Fan","Steemers, Frank","Shendure, Jay","Trapnell, Cole"],"year":"2020","journal":"Science","volume":"367","pages":"45--51","doi":"10.1126/science.aax6234","url":"https://doi.org/10.1126/science.aax6234"},
        {"id":"Tahoe","type":"article","title":"Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling","authors":["Zhang, Jesse","Ubas, Airol A","de Borja, Richard","Svensson, Valentine","Thomas, Nicole","Thakar, Neha","Lai, Ian","Winters, Aidan","Khan, Umair","Jones, Matthew G.","Thompson, John D.","Tran, Vuong","Pangallo, Joseph","Papalexi, Efthymia","Sapre, Ajay","Nguyen, Hoai","Sanderson, Oliver","Nigos, Maria","Kaplan, Olivia","Schroeder, Sarah","Hariadi, Bryan","Marrujo, Simone","Salvino, Crina Curca Alec","Gallareta Olivares, Guillermo","Koehler, Ryan","Geiss, Gary","Rosenberg, Alexander","Roco, Charles","Merico, Daniele","Alidoust, Nima","Goodarzi, Hani","Yu, Johnny"],"year":"2025","journal":"bioRxiv","doi":"10.1101/2025.02.20.639398","url":"https://doi.org/10.1101/2025.02.20.639398"}
    ]

    for article in pmc_articles:
        refs.extend(article['references'])

    return refs


async def stage_sections(client: AsyncOpenAI, model:str, geo_accession:str, pmc_articles:list, labelled_samples:pd.DataFrame, enrichr_results:dict[str,dict[str,pd.DataFrame]], perturbseqr_results:dict[str,pd.DataFrame], methods, results):
    pmc_arrticles = [f'{article["title"]}{article["abstract"]}{article["body"]}' for article in pmc_articles]

    title_prompt = '''Create a title for the signature being re-analyzed by identifying the condition(s)
    that are being used as perturbations in the labelled samples. Return only the name of the signature. For example:
    "Re-Analysis of GSE299362 Human iPSC-Derived Hepatocyte-Like Cells Treated with DHT and Insulin". Fill in the
    GEO study accession and adjust the condition to match the samples used in signature generation. Selec the most
    relevant terms to inform a reader what conditions were used to create a signature. Use proper title capitalization
     for important words, maintaining acronyms, if used. Do not include quantities or special symbols.'''

    REPORT_TITLE,introduction, discussion = await asyncio.gather(
        generate_section(client, model, title_prompt, f'GEO study:{geo_accession}\nSamples:{labelled_samples}'),
        write_report_introduction(client, model, geo_accession, pmc_articles, labelled_samples),
        write_report_discussion(client, model, geo_accession, pmc_articles, labelled_samples, enrichr_results, perturbseqr_results)
    )

    abstract = await write_report_abstract(client, model, geo_accession, introduction, methods, results, discussion)

    return REPORT_TITLE,abstract,introduction,discussion


# LaTex cleanup helper functions
def _repair_unbalanced_brackets(text: str):
    """
        [[[key]]   → [[[key]]]
        [[key]]]   → [[[key]]]
        [[[key]]]  → [[[key]]]
    """
    result = []
    i = 0
    n = len(text)

    while i < n:
        if text[i] != '[':
            result.append(text[i])
            i += 1
            continue

        open_start = i
        open_count = 0
        while i < n and text[i] == '[':
            open_count += 1
            i += 1

        content_start = i
        while i < n and text[i] not in '[]':
            i += 1
        content = text[content_start:i]

        close_count = 0
        while i < n and text[i] == ']':
            close_count += 1
            i += 1

        if close_count == 0:
            result.extend((text[open_start:content_start], content))
            continue

        # Balance both sides to whichever count is larger, capped at 3.
        balanced = min(max(open_count, close_count), 3)
        result.append('[' * balanced + content + ']' * balanced)

    return ''.join(result)


def _normalise_brackets(text: str):
    text = re.sub(r'\[{4,}([^\[\]]+)\]{4,}', r'[[[\1]]]', text)
    text = re.sub(r'(?<!\[)\[\[([^\[\]]+)\]\](?!\])', r'[[[\1]]]', text)
    text = re.sub(r'(?<!\[)\[([A-Za-z][A-Za-z0-9.\-_]{1,40})\](?!\])', r'[[[\1]]]', text)

    return text

_CITE_TOKEN = re.compile(r'\[\[\[([^\[\]]+)\]\]\]')

def _emit_cite(keys: list[str]):
    seen: set[str] = set()
    deduped = [k for k in keys if not (k in seen or seen.add(k))]
    return f'[[[{",".join(deduped)}]]]'

def _merge_adjacent_citations(text: str):
    """[[[a]]][[[b]]] → [[[a,b]]]"""
    tokens = _CITE_TOKEN.split(text)
    out_parts: list[str] = []
    pending_keys: list[str] = []
    pending_gap: str = ''

    i = 0
    while i < len(tokens):
        if i % 2 == 0:
            segment = tokens[i]
            if pending_keys:
                if segment.strip() == '':
                    pending_gap += segment
                else:
                    out_parts.append(_emit_cite(pending_keys))
                    pending_keys = []
                    out_parts.append(pending_gap)
                    pending_gap = ''
                    out_parts.append(segment)
            else:
                out_parts.append(segment)
        else:
            new_keys = [k.strip() for k in tokens[i].split(',') if k.strip()]
            if pending_keys:
                pending_keys.extend(new_keys)
            else:
                pending_keys = new_keys
                pending_gap = ''
        i += 1

    if pending_keys:
        out_parts.append(_emit_cite(pending_keys))

    return ''.join(out_parts)

def _split_comma_citations(text: str):
    """Strip stray spaces around keys inside [[[...]]]."""
    def _clean_keys(match: re.Match) -> str:
        keys = [k.strip() for k in match[1].split(',') if k.strip()]
        return f'[[[{",".join(keys)}]]]' if keys else ''

    return _CITE_TOKEN.sub(_clean_keys, text)

def _extract_reference_keys(references: list[dict]):
    return {ref["id"] for ref in references if "id" in ref}

def _find_unresolved_bare_brackets(text: str, valid_keys: set[str]):
    """Find [key] patterns that survived normalisation and match a valid key."""
    found: list[str] = []
    pat = re.compile(r'(?<!\[)\[([A-Za-z][A-Za-z0-9.\-_]{1,40})\](?!\])')
    for match in pat.finditer(text):
        key = match.group(1)
        if key in valid_keys:
            log(f'Bare bracket citation found for known key: {key}')
            found.append(key)
    return found

def _validate_and_repair_citations(
    text: str,
    valid_keys: set[str],
    section_name: str = '',
    strip_unknown: bool = False,
):
    """Warn on unknown keys; optionally remove them."""
    def _check(match: re.Match) -> str:
        keys = [k.strip() for k in match[1].split(',')]
        good = [k for k in keys if k in valid_keys]
        bad  = [k for k in keys if k not in valid_keys]

        if bad:
            log(f'Unknown citation key(s) in {section_name or "unknown section"}: {"", "".join(bad)}')

        kept = good if strip_unknown else keys
        return f'[[[{",".join(kept)}]]]' if kept else ''

    text = _CITE_TOKEN.sub(_check, text)
    text = re.sub(r'~\s*$', '', text, flags=re.MULTILINE)
    return text


def _escape_latex_segment(segment: str) -> str:
    LATEX_ESCAPES = str.maketrans({
        '%': r'\%',
        '$': r'\$',
        '&': r'\&',
    })

    SUPER_SUB = re.compile(r'([_^])([^{\\])')
    segment = segment.translate(LATEX_ESCAPES)
    segment = SUPER_SUB.sub(r'\1{\2}', segment)
    return segment.encode('ascii', errors='ignore').decode('ascii')

def _to_latex_cite_and_escape(text: str):
    """
    Single-pass conversion: split on [[[key]]] tokens, escape the literal
    segments, and replace citation tokens with ~\\cite{keys}.

    Doing both in one pass guarantees that BibTeX keys are never fed through
    the LaTeX escaper.
    """
    parts = _CITE_TOKEN.split(text)
    out: list[str] = []
    for i, part in enumerate(parts):
        if i % 2 == 0:
            out.append(_escape_latex_segment(part))
        elif keys := [k.strip() for k in part.split(',') if k.strip()]:
            out.append(f'~\\cite{{{",".join(keys)}}}')
    return ''.join(out)

def repair_and_validate_report(
    report: dict,
    references: list[str],
    strip_unknown: bool = False,
):
    valid_keys = _extract_reference_keys(references)

    _TEXT_PATHS: list[tuple[str, ...]] = [
        ('abstract',),
        ('introduction', 'problem'),
        ('introduction', 'background'),
        ('introduction', 'motivation'),
        ('discussion', 'enrichr'),
        ('discussion', 'perturbseqr'),
        ('discussion', 'conclusion'),
    ]


    def _get(d: dict, path: tuple) -> str:
        for key in path:
            d = d[key]
        return d  # type: ignore[return-value]


    def _set(d: dict, path: tuple, value: str) -> None:
        for key in path[:-1]:
            d = d[key]
        d[path[-1]] = value

    for path in _TEXT_PATHS:
        raw: str = _get(report, path)
        section_name = '.'.join(path)

        # 1. Balance and normalise bracket variants → [[[key]]]
        fixed = _repair_unbalanced_brackets(raw)
        fixed = _normalise_brackets(fixed)

        # 2. Merge adjacent citations and clean key whitespace
        fixed = _merge_adjacent_citations(fixed)
        fixed = _split_comma_citations(fixed)

        # 3. Validate keys (and optionally strip unknown ones)
        fixed = _validate_and_repair_citations(fixed, valid_keys, section_name, strip_unknown)

        # 4. Catch any bare [key] that survived normalisation
        bare = _find_unresolved_bare_brackets(fixed, valid_keys)
        for key in bare:
            fixed = re.sub(
                r'(?<!\[)\[' + re.escape(key) + r'\](?!\])',
                f'[[[{key}]]]',
                fixed,
            )
            fixed = _validate_and_repair_citations(fixed, valid_keys, section_name, strip_unknown)

        # 5. Emit \cite{} and escape LaTeX specials in ONE pass (prevents escaper from mangling BibTeX keys like Smith_2020)
        fixed = _to_latex_cite_and_escape(fixed)

        _set(report, path, fixed)

    return report['abstract'], report['introduction'], report['discussion']

    client = AsyncOpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    model = os.getenv("OPENAI_MODEL", "gpt-5-nano")

def construct_georeanalysis_report(geo_accession, pmc_set, labelled_samples_anndata, signature, plots, enrichr_up, enrichr_down, perturbseqr):
    labelled_samples = extract_labelled_samples(labelled_samples_anndata).to_json()
    supplement,enrichr_results = extract_enrichr_results(enrichr_up,enrichr_down)
    perturbseqr_results = extract_perturbseqr_results(perturbseqr)
    pmc_articles = parse_pmc_xml(geo_accession,pmc_set)
    supplement['perturbseqrUpGenes'] = perturbseqr['up_id']
    supplement['perturbseqrDownGenes'] = perturbseqr['down_id']
    log("PMC articles retrieved.")

    log("Generating sections...")
    methods = write_report_methods(geo_accession, labelled_samples)
    results = write_report_results(geo_accession, labelled_samples, signature, enrichr_results, perturbseqr_results)
    title,abstract,introduction,discussion = asyncio.run(stage_sections(client, model, geo_accession, pmc_articles, labelled_samples, enrichr_results, perturbseqr_results, methods, results))
    log("Sections complete.")

    log("Collecting references...")
    references = make_references(pmc_articles)
    log("References complete.")

    abstract, introduction, discussion = repair_and_validate_report(
        dict(
            abstract=abstract,
            introduction=introduction,
            discussion=discussion
        ),
        references,
        strip_unknown=False
    )

    log("Formatting validated.")
    figures = make_figures(plots, enrichr_results, supplement)
    log("Figures complete.")
    tables = make_tables(perturbseqr_results, supplement)

    return {
        "geo_accession":geo_accession, 
        "title":title,
        "abstract":abstract,
        "introduction":introduction,
        "methods":methods,
        "results":results,
        "discussion":discussion,
        "figures":figures,
        "tables":tables,
        "references":references,
        "model":model
    }
