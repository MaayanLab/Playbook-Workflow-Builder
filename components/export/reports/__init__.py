from components.data.gene_count_matrix import anndata_from_file
from components.core.file import upsert_file,file_as_path

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
from aiohttp import ClientSession, ClientError
from openai import AsyncOpenAI
from xml.etree import ElementTree as ET

from adjustText import adjust_text
import base64
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from tempfile import gettempdir,TemporaryDirectory
import subprocess
import shutil
import jinja2
import pypdf
import urllib.request
from io import BytesIO
import json

def log(message: str):
    print(message, file=sys.stderr, flush=True)


GEO_SYSTEM_PROMPT = """You are a scientific writing assistant generating sections of a bioinformatics re-analysis report.

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

CROSSING_SYSTEM_PROMPT = """You are a scientific writing assistant generating sections of a bioinformatics report on the intersection of gene sets derived from NIH Common Fund program datasets.

    The workflow identifies statistically significant crossings between two to five gene sets sourced 
    from Common Fund program libraries, computes overlap statistics (Fisher's exact test p-value, 
    Jaccard index, and the intersecting gene list), submits the overlapping gene set to Enrichr for 
    enrichment analysis across biological process, pathway, transcription factor, and phenotype 
    libraries, and synthesizes per-gene literature deepdives drawn from the top 50 most-cited PubMed 
    abstracts for each gene in the intersection.

    Input provided for each section will include some combination of: the names and sizes of the 
    crossing gene sets, crossing statistics, Enrichr enrichment results, dataset and Common Fund 
    program metadata, and per-gene deepdive summaries compiled from PubMed abstracts. 
    Write only from what is explicitly provided in the input for that section.

    Output rules that apply to every response:
    - Use only plain text prose. No markdown, no bullet points, no numbered lists, no bold, no italics, no headers. Special characters (including em-dashes) will be removed and negatively affect formatting. so do not use them.
    - When an instruction explicitly asks for an outline, use a dash to begin each point. In all other cases, write in paragraph prose only.
    - Cite sources inline using this exact format (no other brackets are acceptable): [[[key1, key2]]]. Do not invent citation keys, do not reformat them. Citations you include MUST come only from input publications, results, or explicit instructions.
    - If you are explicitly instructed to exclude citations for a section, you must do so.
    - Do not speculate or introduce claims beyond what the provided source material explicitly supports.
    - Do not infer biological meaning, mechanistic relationships, or experimental context from gene set names, term labels, or dataset identifiers alone. A term name is an identifier, not evidence for a biological claim.
    - Do not introduce prior knowledge about specific genes, pathways, or datasets beyond what the provided input explicitly states. This applies in particular to deepdive synthesis: draw only from the provided per-gene abstract summaries and do not supplement them with external claims about those genes.
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


async def generate_section(client: AsyncOpenAI, model:str, system_prompt:str, prompt:str, data:str):
    
    history = [
        {"role":"system","content":system_prompt},
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

    abstract = await generate_section(client, model, GEO_SYSTEM_PROMPT, prompt, data)
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
        generate_section(client, model, GEO_SYSTEM_PROMPT, problem_prompt, pmc_articles),
        generate_section(client, model, GEO_SYSTEM_PROMPT, signature_prompt, labelled_samples)
    )
    introduction_background = await generate_section(client, model, GEO_SYSTEM_PROMPT, background_prompt, f"Articles:{pmc_articles}\nIntroduction Problem:{introduction_problem}")

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
        generate_section(client, model, GEO_SYSTEM_PROMPT, enrichr_prompt, str(enrichr_results)),
        generate_section(client, model, GEO_SYSTEM_PROMPT, perturbseqr_prompt, str(perturbseqr_results)),
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
    discussion_conclusion = await generate_section(client, model, GEO_SYSTEM_PROMPT, conclusion_prompt, conclusion_data)

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

    ext_link = cit.find(".//ext-link")
    ext_url = ""
    if ext_link is not None:
        ext_url = (
            ext_link.get("{http://www.w3.org/1999/xlink}href")
            or ext_link.get("xlink:href")
            or ""
        ).strip()

    url = ext_url or (f"https://doi.org/{doi}" if doi else "")

    title = (cit.findtext("article-title") or "").strip()
    if not title and pub_type not in ("journal",):
        text_parts = ([cit.text] if cit.text else []) + [
            child.tail for child in cit if child.tail
        ]
        title = " ".join(p.strip() for p in text_parts if p.strip())
        title = " ".join(title.split())

    return _ref_dict(entry_type, ref_id, {
        "title":   title,
        "authors": _extract_citation_author_list(cit),
        "journal": (cit.findtext("source") or "").strip(),
        "year":    cit.findtext("year", ""),
        "volume":  cit.findtext("volume", ""),
        "pages":   _page_range(cit),
        "doi":     doi,
        "pmid":    cit.findtext(".//pub-id[@pub-id-type='pmid']", ""),
        **({"url": url} if url else {}),
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
    with upsert_file('.pdf') as f:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
    return f


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
    with upsert_file('.pdf') as f:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
    return f


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
    with upsert_file('.pdf') as f:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
    return f


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
    with upsert_file('.pdf') as f:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
    return f


def make_figures(plots, enrichr_results:dict[str,dict[str,pd.DataFrame]], supplement:dict[str,str]):
    library_pdf = make_library_size_barplot(plots["library_sizes_plot"])
    pca_pdf = make_pca_scatter(plots["pca_plot"])
    volcano_pdf = make_volcano_scatter(plots["volcano_plot"])
    enrichr_up_pdf = make_enrichr_barplot(enrichr_results["enrichr_up"])
    enrichr_down_pdf = make_enrichr_barplot(enrichr_results["enrichr_down"])
    enrichr_up_id = supplement["enrichrUp"]
    enrichr_down_id = supplement["enrichrDown"]
    return {
        "librarySizes": {
            "file":library_pdf,
            "caption":"Library sizes for each sample in the dataset, shown as total mapped read counts per sample. Samples are labeled by their GEO accession identifier (GSM). Consistent library sizes across samples indicate uniform sequencing depth suitable for differential expression analysis."
        },
        "PCAScatter": {
            "file":pca_pdf,
            "caption":"Principal component analysis (PCA) of normalized gene expression profiles across all samples. Each point represents one sample, colored by experimental condition: control (blue), perturbation (red), and additional study samples not included in the differential expression analysis (gray). Axes indicate the percentage of total variance explained by each principal component. Normalization was performed using log-counts-per-million (logCPM)."
        },
        "volcanoScatter": {
            "file":volcano_pdf,
            "caption":"Volcano plot of differential expression results comparing perturbation to control samples. Each point represents one protein-coding gene; the x-axis shows the log2 fold-change and the y-axis shows statistical significance as -log10(adjusted p-value). Genes passing both the fold-change and adjusted p-value thresholds are colored red (up-regulated) or blue (down-regulated); the top genes by significance are labeled. Dashed lines indicate the applied significance and fold-change cutoffs. Insignificant points are randomly downsampled by a factor of 0.5.",
        },
        "upEnrichrBars": {
            "file":enrichr_up_pdf,
            "caption":f"Enrichment analysis of the up-regulated gene signature using Enrichr~\\cite{{Enrichr}}. Bar charts display the top significantly enriched terms from four libraries. A. GO Biological Process 2023~\\cite{{GO}}, B. KEGG 2021 Human~\\cite{{KEGG}}, C. ChEA 2022~\\cite{{ChEA}}, and D. KOMP2 Mouse Phenotypes 2022~\\cite{{KOMP2}}. Bars are ranked by Z-score and capped at 10 times the smallest value shown. Color indicates the source library. The full enrichment results are available to view at \\href{{https://maayanlab.cloud/enrichr/enrich?dataset={enrichr_up_id}}}{{Enrichr}}.",
        },
        "downEnrichrBars": {
            "file":enrichr_down_pdf,
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
        generate_section(client, model, GEO_SYSTEM_PROMPT, title_prompt, f'GEO study:{geo_accession}\nSamples:{labelled_samples}'),
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

    _UNICODE_TO_LATEX: list[tuple[str, str]] = [
        # Dashes
        ('\u2014', '-'),
        ('\u2013', '-'),
        ('\u2012', '--'),
        ('\u2015', '---'),
        # Spaces
        ('\u00a0', ' '),
        ('\u202f', ' '),
        ('\u2009', ' '),
        ('\u2003', ' '),
        ('\u2002', ' '),
    ]

    for char, replacement in _UNICODE_TO_LATEX:
        if char in segment:
            segment = segment.replace(char, replacement)

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

    _TEXT_PATHS: list[tuple[str, ...]] = [('abstract',)]+ \
    [('introduction', key) for key in report["introduction"].keys()]+ \
    [('discussion', key) for key in report["discussion"].keys()]


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
            if path==('abstract',):
                fixed = re.sub(r'\[+\s*[A-Za-z0-9_][^\[\]]*?\s*\]+', '', fixed)
            fixed = _validate_and_repair_citations(fixed, valid_keys, section_name, strip_unknown)

        # 5. Emit \cite{} and escape LaTeX specials in ONE pass (prevents escaper from mangling BibTeX keys like Smith_2020)
        fixed = _to_latex_cite_and_escape(fixed)

        _set(report, path, fixed)

    return report['abstract'], report['introduction'], report['discussion']

def construct_georeanalysis_report(geo_accession, pmc_set, labelled_samples_anndata, signature, plots, enrichr_up, enrichr_down, perturbseqr):
    client = AsyncOpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    model = os.getenv("OPENAI_MODEL", "gpt-5-nano")
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
        "supplement":supplement,
        "references":references,
        "model":model
    }

def resolve_drs_url(url:str) -> str:

    return url


def construct_geo_report_references(references: list[dict]) -> str:
    entries = []
    for ref in references:
        fields = [
            f"  title     = {{{ref['title']}}}" if ref.get('title') else None,
            f"  author    = {{{' and '.join(ref['authors'])}}}" if ref.get('authors') else None,
            f"  year      = {{{ref['year']}}}"    if ref.get('year')    else None,
            f"  journal   = {{{ref['journal']}}}" if ref.get('journal') else None,
            f"  volume    = {{{ref['volume']}}}"  if ref.get('volume')  else None,
            f"  pages     = {{{ref['pages']}}}"   if ref.get('pages')   else None,
            f"  doi       = {{{ref['doi']}}}"     if ref.get('doi')     else None,
            f"  url       = {{{ref['url']}}}"     if ref.get('url')     else None,
        ]
        fields = [f for f in fields if f is not None]
        entries.append(f"@{ref['type']}{{{ref['id']}, {', '.join(fields)} }}")

    return '\n\n'.join(entries)

def render_georeanalysis_report(report):
    log('Preparing LaTeX bundle...')
    geo_accession = report["geo_accession"]
    with TemporaryDirectory() as texdir:
        public_url = os.getenv("PUBLIC_URL")
        urllib.request.urlretrieve(f"{public_url}/geo-report/geo_report.tex", f"{texdir}/geo_report.tex")
        urllib.request.urlretrieve(f"{public_url}/geo-report/ExcelAtFIT.cls", f"{texdir}/ExcelAtFIT.cls")
        urllib.request.urlretrieve(f"{public_url}/geo-report/maayan-lab-logo.png", f"{texdir}/maayan-lab-logo.png")
        os.makedirs(f"{texdir}/figures", exist_ok=True)
        for figure,figurefile in report["figures"].items():
            with file_as_path(figurefile["file"]) as path:
                shutil.copy(path, f"{texdir}/figures/{figure}.pdf")
        os.makedirs(f"{texdir}/tables", exist_ok=True)
        for table,tablefile in report["tables"].items():
            with file_as_path(tablefile["file"]) as path:
                shutil.copy(path, f"{texdir}/tables/{table}.tsv")
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(texdir),
            variable_start_string='(((',
            variable_end_string=')))',
        )
        template = env.get_template("geo_report.tex").render(
            title=report['title'],
            abstract=report['abstract'],
            introduction = report['introduction'],
            methods = report['methods'],
            results = report['results'],
            discussion = report['discussion'],
            figures = report['figures'],
            tables = report['tables'],
            model = report['model'],
        )
        log('Rendering LaTeX template...')
        with open(f"{texdir}/{geo_accession}.tex", "w") as texfile:
            texfile.write(template)
        with open(f"{texdir}/references.bib", "w") as bibfile:
            bibfile.write(construct_geo_report_references(report["references"]))
        
        with upsert_file('.zip') as zipfile:
            base = zipfile.file.rsplit('.', 1)[0]
            shutil.make_archive(base, "zip", texdir)
        log('Bundle complete.')
        log('Beginning PDF compilation...')
            
        with subprocess.Popen(
            ['latexmk', '-pdf', '-lualatex', '-cd', '-f', '-interaction=nonstopmode', f"{texdir}/{geo_accession}.tex"],
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True) as proc:
            for line in proc.stdout:
                if line.startswith("Run"):
                    log(line.strip())
        
        with upsert_file(".pdf") as pdffile:
            reader = pypdf.PdfReader(f'{texdir}/{geo_accession}.pdf')
            writer = pypdf.PdfWriter()
            for page in reader.pages:
                writer.add_page(page)
            writer.write(pdffile.file)
        
    zipfile["filename"] = f"{geo_accession}.zip"
    pdffile["filename"] = f"{geo_accession}.pdf"

    return {"bundle":zipfile, "pdf":pdffile}

def make_venn_diagram(venn_diagram_plotly):
    img_source = venn_diagram_plotly['layout']['images'][0]['source']
    b64_str = img_source.split(',', 1)[1]
    img_bytes = base64.b64decode(b64_str)

    img = plt.imread(BytesIO(img_bytes))

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    ax.imshow(img)
    ax.axis('off')

    plt.tight_layout()
    with upsert_file('.pdf') as f:
        plt.savefig(f.file, format='pdf', dpi=300, bbox_inches='tight')
    plt.close()
    return f

def get_gpt4o_gene_summary(gene_symbol: str) -> str | None:
    """Fetch just the GPT-4o gene summary from GSFM via tRPC."""
    url = "https://gsfm.maayanlab.cloud/api/trpc/gene_info"
    params = {
        "input": json.dumps(gene_symbol)
    }
    
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    
    data = resp.json()
    gene_info = data["result"]["data"]
    
    raw = gene_info.get("deepdive_gpt4o_description")
    if not raw:
        return None
    
    return json.loads(raw)

def format_deepdive(symbol, deepdive):
    def md(props):
        fn = COMPONENTS.get(props.get("type"))
        return fn(props) if fn else ""

    def component_recurse(props):
        return "".join(md(child) for child in props.get("children", []))

    COMPONENTS = {
        "root": lambda props: component_recurse(props).strip(),

        "p": lambda props: component_recurse(props).strip() + "\n\n",

        "t": lambda props: props.get("text", ""),

        # flatten string
        "b": lambda props: component_recurse(props),
        "i": lambda props: component_recurse(props),

        # keep link text
        "a": lambda props: component_recurse(props),

        # single reference
        "f": lambda props: f"[{props.get('ref')}]",

        # grouped references
        "fg": lambda props: (
            "[" +
            ",".join(md(child) for child in props.get("children", [])) +
            "]"
        ),

        # individual ref inside fg
        "fg_f": lambda props: props.get("ref"),

        # reference range
        "fg_fs": lambda props: (
            f"{props.get('start_ref')}-{props.get('end_ref')}"
        ),

        # references group
        "rg": lambda props: (
            "\nReferences\n\n"
            #+ "-" * 10 + "\n"
            + component_recurse(props).strip()
        ),

        # reference
        "r": lambda props: (
            f"[{props.get('ref')}] "
            f"{component_recurse(props).strip()}\n"
        ),
    }

    rendered = md(deepdive).rstrip().split('\n\n\n', 1)[0]

    ref_to_doi = {}

    def extract_dois(node):
        if not isinstance(node, dict):
            return

        node_type = node.get("type")

        if node_type == "r":
            ref_number = node.get("ref")
            doi = None

            children = node.get("children", [])
            for i, child in enumerate(children):
                # look for "DOI: " text node
                if (
                    child.get("type") == "t"
                    and "DOI" in child.get("text", "")
                ):
                    if i + 1 < len(children):
                        next_child = children[i + 1]
                        if next_child.get("type") == "a":
                            doi = "".join(
                                c.get("text", "")
                                for c in next_child.get("children", [])
                                if c.get("type") == "t"
                            )
                    break

            if ref_number is not None and doi:
                ref_to_doi[f'{symbol}_{int(ref_number)}'] = doi

        for child in node.get("children", []):
            extract_dois(child)

    extract_dois(deepdive)

    return rendered,ref_to_doi

class CrossingReport():
    '''
    Organize and generate the sections of a crossing report, analyzing a selecgted crossing by integrating
    information about the crossed GMTs, involved terms, enrichment analysis, and research DeepDives of the
    genes from the intersecting set.
    '''
    def __init__(self, client: AsyncOpenAI, model:str, datasets, gmts, crossing, venn_diagram, intersecting_genes, enrichr):
        self.client = client
        self.model = model
        self.system_prompt = CROSSING_SYSTEM_PROMPT
        self.datasets = datasets
        self.gmts = gmts
        self.crossing = crossing
        self.venn_diagram = venn_diagram
        self.intersecting_genes = intersecting_genes
        self.enrichr = enrichr
        self.title = None
        self.abstract = None
        self.introduction = None
        self.methods = None
        self.results = None
        self.discussion = None
        self.references = []
        log("Initialized")


    async def write_report_abstract(self):
        log("Writing abstract...")
        prompt = dedent('''
        Write a single plain text paragraph abstract of approximately 200 words for a gene set crossing 
        analysis report. You should introduce the problem of discrete datasets, summarize why crossing
        Common Fund datasets is useful, what was done, and what was found. You must make sure to include
        the name of each resource and why each included gene set is useful based on the results, using
        the discussion section as a guide to justify your reasoning.Your abstract should serve as an overall
        summary of the report. Do not include statistics, a header, citations, or any symbols. Highlight that the
        report explores one crossing, that the intersectiung genes are in all gene sets and that the crossing
        is being investigated because of their maximal inclusion among all possible overlaps, and highlight
        why the relevant drug may be useful. Your response should closely resemble this example format:
        Biomedical omics resources are often analyzed as discrete datasets, which fragments biological
        insight and limits discovery of shared mechanisms across conditions and species. Crossing Common
        Fund datasets can reveal molecular programs, increase robustness, assess reproducibility, and
        connect molecular changes to physiology and intervention. We intersected four complementary
        gene sets spanning intervention, aging, tissue state, and metabolism: an imeglimin response in
        adipose from RummaGEO, a heart aging signature from GTEx, a brown adipose exercise consensus from
        MoTrPAC, and an metabolite centered associations-derived NAD+ enzyme gene set from Metabolomics
        Workbench. Gene sets were harmonized and crossed, and the intersection was profiled for pathway
        and regulator enrichment. The intersection converged on a coherent mitochondrial oxidative
        metabolism program that linked imeglimin exposure, thermogenic biology, and age related expression
        shifts. Shared genes included pyruvate dehydrogenase complex subunits PDHA1 PDHB DLAT DLD and DLST,
        tricarboxylic acid cycle enzymes IDH3B MDH1 MDH2, respiratory chain components NDUFS1 NDUFS2, and LDHB.
        Enriched processes highlighted aerobic respiration, pyruvate metabolism, acetyl CoA biosynthesis, the
        malate aspartate shuttle, NADH oxidation, and lipoic acid metabolism, consistent with tighter coupling
        of pyruvate entry to citrate cycle throughput and complex one activity under NAD plus centered redox control.
        The overlap with the heart aging signature suggests partial restoration of youthful mitochondrial
        programs, while the brown adipose consensus supports increased thermogenic capacity. Together, these
        findings nominate imeglimin as a potential modulator of aging, exercise-like adaptations, and metabolic
        flexibility.

        Do not copy sentences verbatim from the provided sections.
        You MUST NOT include citations or reference keys of any kind.''').strip()
        data = dedent(f'''
        Introduction:{self.introduction}
        Methods: {self.methods}
        Results: {self.results}
        Discussion: {self.discussion}''').strip()

        abstract = await generate_section(self.client, self.model, self.system_prompt, prompt, data)
        abstract = re.sub(r"\[\[\[(.+?)\]\]\]", "", abstract)
        abstract = re.sub(r"~\\cite\{(?:.+?)\}","", abstract).strip(',')

        log("Abstract complete.")
        self.abstract = abstract
        return abstract


    async def write_report_introduction(self):
        log("Writing introduction...")
        terms = [val for (key,val) in self.crossing.items() if ("term" in key and "Length" not in key)]
        dataset_info = [(dataset["dataset"]["name"],dataset["dataset"]["resource"]) for dataset in self.datasets]

        static_introduction = dedent('''
        Common Fund programs have generated many diverse and robust multiomics datasets in the past 20 years.
        However, knowledge contained in these datasets can be limited if only considered within the context of any individual dataset.
        To address this, we have converted Common Fund datasets from the Common Fund Data Ecosystem (CFDE) [[[CFDE}]]] into gene set libraries  to combine multiple datasets at a time and investigate interesting overlapping features.
        Pairwise gene set library crossing has been explored previously in GeneSetCart [[[GeneSetCart]]], and Harmonizome 3.0 [[[Harmonizome]]], where gene sets from two libraries where significant overlaps were found using Fisher's exact test.
        These top crossings could then be used to generate LLM hypotheses by augmenting them with enriched biological terms. Here we extend the strategy by combining more libraries at once.
        Each combination of three gene sets (one from each library) is analyzed to identify signficant pairwise intersections and a robust 3-way intersecting gene set.
        ''')

        datasets = ", and ".join(", ".join([f'{name} [[[{resource}]]]' for (name, resource) in dataset_info]).rsplit(", ",1))
        dynamic_introduction = dedent(f'''
        For this report, we combined gene sets from the {datasets} datasets.
        Among the top-ranked crossings, we identified an overlap with rank {self.crossing["rank"]} and a p-value of {self.crossing["pvalue"]:.5e} warranting further investigation.
        ''')

        terms_prompt = dedent('''
        Introduce the gene sets found in this top-ranked crossing between 3 datasets.
        You must make sure to include each term name in the provided top-ranked crossing to provide context.
        Include a simple explanation of what each gene set represents using only the provided inforatmion.
        Your response will be inserted into a template, so do not inclue an introduction or conclusion.
        Aovid repeating the names of each geene set and using the term "signature" instead of "gene set".
        Make sure to mention that the GlyGen set represents proteins phosphorylated by the glycan.
        Return LaTeX compliant plaintext paragraph.
        Do not use any LaTeX or Markdown expressions.''').strip()

        introduction_terms = await generate_section(
            self.client,
            self.model,
            self.system_prompt,
            terms_prompt,
            f"Terms:{terms}\nIntroduction:{static_introduction}\n{dynamic_introduction}\nDatasets:{dataset_info}",
        )

        introduction = {
            "cfde":static_introduction, 
            "dataset":dynamic_introduction,
            "terms":introduction_terms
        }

        log("Introduction complete.")
        self.introduction = introduction
        return introduction


    def write_report_methods(self):
        dataset_info = [
            (
                dataset["dataset"]["name"],dataset["dataset"]["resource"],dataset["dataset"]["processing"]
            ) for dataset in self.datasets
        ]
        dataset_names = ", and ".join(", ".join([name for (name,_,_) in dataset_info]).rsplit(",", 1))
        dataset_strings = "\n".join([dataset_processing for (_,_,dataset_processing) in dataset_info])
        methods = dedent(f'''
            {dataset_strings}
            The workflow starts with selecting the {dataset_names} datasets. Gene matrix transpose (GMT) files were retrieved from the CFDE
            Workbench ~\\cite{{CFDE}}. GMTs were crossed to analyze all possible combinations of one gene set from each library. For each combination,
            the p-value was computed by computing the hypergeometric test of each gene set against the intersection of the remaining sets, using the
            most conservative p-value as the value for the crossing. The top 5,000 combinations were kept, ranked by ascending (most conservative) p-value
            and descending Jaccard index. The crossing with rank {self.crossing["rank"]} was selected for further investigation, and the intersecting gene
            set was extracted. The terms involved in the crossing were extracted from their respective gene set library and used to create a {len(dataset_info)}-way
            Venn diagram. The intersecting gene set was submitted to Enrichr ~\\cite{{Enrichr}}. The gene sets were enriched against the GO Biological
            Process 2023 ~\\cite{{GO}}, KEGG 2021 Human ~\\cite{{KEGG}}, ChEA 2022 ~\\cite{{ChEA}}, and GWAS Catalog 2019 ~\\cite{{GWAS}} libraries to
            identify statistically significant enriched biological processes, pathways, transcription factors and phenotypes.'''
        ).strip().replace('\n', ' ')
        
        self.methods = methods
        return methods


    def write_report_results(self,  deepdive_results):
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
                ("gwas", lambda: f"GWAS Catalog phenotypes associated with {joined_terms.lower()} ~\\cite{{GWAS}}"),
            ]

            for pattern, formatter in standard_rules:
                if pattern in lib_lower:
                    return formatter()

            return f"{joined_terms} from {library}"

        def format_enrichr_section(results, n_terms=2):
            library_summaries = [
                summarize_library(lib.rsplit("_", 1)[0], terms_df["term"], n_terms=n_terms)
                for lib, terms_df in results.items()
            ]

            return natural_join(library_summaries)

        results = dedent(f'''
        The crossing had a rank of {self.crossing["rank"]}, p-value of {self.crossing["pvalue"]}, and Jaccard index of {self.crossing["jaccard"]}.
        The intersecting gene set identified {self.crossing["overlap"]} genes appearing in all involved sets (Figure \\ref{{fig:venn}}).
        The intersecting gne set consists of {", and".join(self.crossing["genes"].rsplit(" ", 1))}.
        Functional enrichment analysis of the gene set identified associations with {format_enrichr_section(self.enrichr)} (Figure \\ref{{fig:enrichment}}). 
        ''').strip().replace('\n',' ')

        self.results = results
        return results


    async def write_report_discussion(self, deepdive_results):
        log("Writing discussion...")

        deepdive_prompt = dedent('''
        Create between 1 and 3 paragraphs summarizing the shared biological roles of this gene set.

        Instructions:
        1. Reason over the similarities to produce a bullet point outline to guide you.
        2. Each paragraph must begin with an introductory sentence describing the main theme.
        3. Each paragraph should include 2-4 points describing supporting biological mechanisms, pathways, or phenotypes.

        Content constraints:
        - Focus on biological processes related to aging and glycosylation.
        - Only use information explicitly present in the provided summaries.
        - Do NOT introduce outside knowledge or speculation.

        Citation rules (VERY IMPORTANT):
        - Every sub-bullet MUST end with the citations that support the statement.
        - Citations must be copied EXACTLY from the summaries in the format:
        [[[GENE_REFNUM]]]
        Examples: [[[A1BG_1, STAT3_3]]]; [[[ABCA2_1, ABCA2_2, TMEM127B_4]]
        - The gene symbol & reference number format are immutable identifiers.
        - Numeric citations like [1] or [1,3] are INVALID.
        - Do NOT renumber, merge, or alter citation tokens.
        - Do NOT group citations by gene; simply list all relevant citations.

        Style:
        - Avoid repeating the same mechanism across sentences.
        - Combine similar gene functions where possible.
        ''').strip()

        enrichr_prompt = dedent('''
        Write a plain text paragraph analyzing the Enrichr enrichment results provided for the up-regulated and down-regulated gene sets from this re-analysis.
        Do not list the terms to introduce them, assume this has already been done.
        Focus on: terms that appear or are consistent across multiple libraries, and any complementary or contrasting patterns between the up and down gene sets.
        Only discuss terms that are present in the results provided.
        Do not introduce terms or biological processes not present in the data.
        Cite only the libraries whose results you are directly referencing, choosing from: 
        [[[GO]]], [[[KEGG]]], [[[ChEA]]], [[[GWAS]]].''').strip()

        discussion_deepdive, discussion_enrichr,  = await asyncio.gather(
            generate_section(self.client, self.model, self.system_prompt, deepdive_prompt, str(deepdive_results)),
            generate_section(self.client, self.model, self.system_prompt, enrichr_prompt, str(self.enrichr)),
        )

        conclusion_prompt = dedent('''
        Write a plain text concluding paragraph for a re-analysis report.
        Summarize the most notable findings from the crossing, gene deepdives, and Enrichr analyses provided.
        Focus on findings that are consistent with or directly relevant to the biology of the terms.
        Where a finding is unexpected relative to the biological context, note it briefly.
        Do not propose mechanisms unless they are directly supported by the enrichment results provided.
        Connect the findings back to the research question described in the DeepDive and Enrichr paragraphs.
        Do not include citations.''').strip()

        conclusion_data = f'''Crossing: {self.crossing}\nEnrichr Results:{discussion_enrichr}\nPerturb-Seqr Results:{discussion_deepdive}'''
        discussion_conclusion = await generate_section(self.client, self.model, self.system_prompt, conclusion_prompt, conclusion_data)

        log("Discussion complete.")

        discussion =  {
            "deepdive":discussion_deepdive,
            "enrichr":discussion_enrichr, 
            "conclusion":discussion_conclusion
        }
        self.discussion = discussion

        return discussion

    async def make_references(self, deepdive_results):
        def extract(ref: str, field: str, default=None):
            m = re.search(rf'{field}\s*=\s*\{{([^}}]*)\}}', ref)
            if not m:
                return [] if field == "author" else default
            value = m.group(1).strip()
            if field == "author":
                return [a.strip() for a in value.split(" and ")]
            return value

        semaphore = asyncio.Semaphore(2)

        async def doi_to_ref(session: ClientSession, key: str, doi: str):
            async with semaphore:
                try:
                    async with session.get(
                        url=f"https://api.crossref.org/works/{doi}/transform/application/x-bibtex?mailto=danieljbclarkemssm@gmail.edu",
                    ) as res:
                        if not res.ok:
                            return {"id": key, "type": "misc", "doi": doi, "url": f"https://doi.org/{doi}"}

                        ref = await res.text(encoding="utf-8")

                        if not ref.startswith("@"):
                            return {"id": key, "type": "article", "doi": doi, "url": f"https://doi.org/{doi}"}

                        return {
                            "id": key,
                            "type": "article",
                            "title": extract(ref, "title"),
                            "year": extract(ref, "year"),
                            "authors": extract(ref, "author"),
                            "journal": extract(ref, "journal"),
                            "volume": extract(ref, "volume"),
                            "pages": extract(ref, "pages"),
                            "doi": doi,
                            "url": f"https://doi.org/{doi}",
                        }

                except ClientError:
                    return {"id": key, "type": "misc", "doi": doi, "url": f"https://doi.org/{doi}"}

                finally:
                    await asyncio.sleep(0.5)


        refs = [
            {"id": "CFDE", "type": "article", "title":"The CFDE Workbench: Integrating Metadata and Processed Data from Common Fund Programs","year":"2026", "authors": ["Evangelista, John Erol", "Clarke, Daniel J.B.", "Byrd, Anna I.", "Srinivasan, Shivaramakrishna", "Srinivasan, Sumana", "Maurya, Mano R.", "Jenkins, Sherry L.", "Diamant, Ido", "Sanchez, Ethan", "Xie, Zhuorui", "Olaiya, Stephanie", "Kim, Heesu", "Marino, Giacomo B.", "Ahmed, Nasheath", "Ramachandran, Srinivasan", "Subramaniam, Shankar", "Ma'ayan, Avi"], "journal":"Journal of Molecular Biology", "pages":"169631", "doi": "10.1016/j.jmb.2026.169631", "url": "https://doi.org/10.1016/j.jmb.2026.169631"},
            {"id": "Harmonizome", "type": "article", "title":"Harmonizome 3.0: integrated knowledge about genes and proteins from diverse multi-omics resources","year":"2024", "authors": ["Diamant, Ido", "Clarke, Daniel\u00a0J B", "Evangelista, John\u00a0Erol", "Lingam, Nathania", "Ma'ayan, Avi"], "journal":"Nucleic Acids Research", "volume":"53", "pages":"D1016--D1028", "doi": "10.1093/nar/gkae1080", "url": "https://doi.org/10.1093/nar/gkae1080"},
            {"id": "GeneSetCart", "type": "article", "title":"GeneSetCart: assembling, augmenting, combining, visualizing, and analyzing gene sets","year":"2025", "authors": ["Marino, Giacomo B", "Olaiya, Stephanie", "Evangelista, John Erol", "Clarke, Daniel J B", "Ma'ayan, Avi"], "journal":"GigaScience", "volume":"14", "doi": "10.1093/gigascience/giaf025", "url": "https://doi.org/10.1093/gigascience/giaf025"},
            {"id": "GlyGen", "type": "article", "title":"GlyGen: Computational and Informatics Resources for Glycoscience","year":"2019", "authors": ["York, William S", "Mazumder, Raja", "Ranzinger, Rene", "Edwards, Nathan", "Kahsay, Robel", "Aoki-Kinoshita, Kiyoko F", "Campbell, Matthew P", "Cummings, Richard D", "Feizi, Ten", "Martin, Maria", "Natale, Darren A", "Packer, Nicolle H", "Woods, Robert J", "Agarwal, Gaurav", "Arpinar, Sena", "Bhat, Sanath", "Blake, Judith", "Castro, Leyla Jael Garcia", "Fochtman, Brian", "Gildersleeve, Jeffrey", "Goldman, Radoslav", "Holmes, Xavier", "Jain, Vinamra", "Kulkarni, Sujeet", "Mahadik, Rupali", "Mehta, Akul", "Mousavi, Reza", "Nakarakommula, Sandeep", "Navelkar, Rahi", "Pattabiraman, Nagarajan", "Pierce, Michael J", "Ross, Karen", "Vasudev, Preethi", "Vora, Jeet", "Williamson, Tatiana", "Zhang, Wenjin"], "journal":"Glycobiology", "volume":"30", "pages":"72--73", "doi": "10.1093/glycob/cwz080", "url": "https://doi.org/10.1093/glycob/cwz080"},
            {"id": "GTEx", "type": "article", "title":"The GTEx Consortium atlas of genetic regulatory effects across human tissues","year":"2020", "authors": ["", "Aguet, Fran\u00e7ois", "Anand, Shankara", "Ardlie, Kristin G.", "Gabriel, Stacey", "Getz, Gad A.", "Graubert, Aaron", "Hadley, Kane", "Handsaker, Robert E.", "Huang, Katherine H.", "Kashin, Seva", "Li, Xiao", "MacArthur, Daniel G.", "Meier, Samuel R.", "Nedzel, Jared L.", "Nguyen, Duyen T.", "Segr\u00e8, Ayellet V.", "Todres, Ellen", "Balliu, Brunilda", "Barbeira, Alvaro N.", "Battle, Alexis", "Bonazzola, Rodrigo", "Brown, Andrew", "Brown, Christopher D.", "Castel, Stephane E.", "Conrad, Donald F.", "Cotter, Daniel J.", "Cox, Nancy", "Das, Sayantan", "de Goede, Olivia M.", "Dermitzakis, Emmanouil T.", "Einson, Jonah", "Engelhardt, Barbara E.", "Eskin, Eleazar", "Eulalio, Tiffany Y.", "Ferraro, Nicole M.", "Flynn, Elise D.", "Fresard, Laure", "Gamazon, Eric R.", "Garrido-Mart\u00edn, Diego", "Gay, Nicole R.", "Gloudemans, Michael J.", "Guig\u00f3, Roderic", "Hame, Andrew R.", "He, Yuan", "Hoffman, Paul J.", "Hormozdiari, Farhad", "Hou, Lei", "Im, Hae Kyung", "Jo, Brian", "Kasela, Silva", "Kellis, Manolis", "Kim-Hellmuth, Sarah", "Kwong, Alan", "Lappalainen, Tuuli", "Li, Xin", "Liang, Yanyu", "Mangul, Serghei", "Mohammadi, Pejman", "Montgomery, Stephen B.", "Mu\u00f1oz-Aguirre, Manuel", "Nachun, Daniel C.", "Nobel, Andrew B.", "Oliva, Meritxell", "Park, YoSon", "Park, Yongjin", "Parsana, Princy", "Rao, Abhiram S.", "Reverter, Ferran", "Rouhana, John M.", "Sabatti, Chiara", "Saha, Ashis", "Stephens, Matthew", "Stranger, Barbara E.", "Strober, Benjamin J.", "Teran, Nicole A.", "Vi\u00f1uela, Ana", "Wang, Gao", "Wen, Xiaoquan", "Wright, Fred", "Wucher, Valentin", "Zou, Yuxin", "Ferreira, Pedro G.", "Li, Gen", "Mel\u00e9, Marta", "Yeger-Lotem, Esti", "Barcus, Mary E.", "Bradbury, Debra", "Krubit, Tanya", "McLean, Jeffrey A.", "Qi, Liqun", "Robinson, Karna", "Roche, Nancy V.", "Smith, Anna M.", "Sobin, Leslie", "Tabor, David E.", "Undale, Anita", "Bridge, Jason", "Brigham, Lori E.", "Foster, Barbara A.", "Gillard, Bryan M.", "Hasz, Richard", "Hunter, Marcus", "Johns, Christopher", "Johnson, Mark", "Karasik, Ellen", "Kopen, Gene", "Leinweber, William F.", "McDonald, Alisa", "Moser, Michael T.", "Myer, Kevin", "Ramsey, Kimberley D.", "Roe, Brian", "Shad, Saboor", "Thomas, Jeffrey A.", "Walters, Gary", "Washington, Michael", "Wheeler, Joseph", "Jewell, Scott D.", "Rohrer, Daniel C.", "Valley, Dana R.", "Davis, David A.", "Mash, Deborah C.", "Branton, Philip A.", "Barker, Laura K.", "Gardiner, Heather M.", "Mosavel, Maghboeba", "Siminoff, Laura A.", "Flicek, Paul", "Haeussler, Maximilian", "Juettemann, Thomas", "Kent, W. James", "Lee, Christopher M.", "Powell, Conner C.", "Rosenbloom, Kate R.", "Ruffier, Magali", "Sheppard, Dan", "Taylor, Kieron", "Trevanion, Stephen J.", "Zerbino, Daniel R.", "Abell, Nathan S.", "Akey, Joshua", "Chen, Lin", "Demanelis, Kathryn", "Doherty, Jennifer A.", "Feinberg, Andrew P.", "Hansen, Kasper D.", "Hickey, Peter F.", "Jasmine, Farzana", "Jiang, Lihua", "Kaul, Rajinder", "Kibriya, Muhammad G.", "Li, Jin Billy", "Li, Qin", "Lin, Shin", "Linder, Sandra E.", "Pierce, Brandon L.", "Rizzardi, Lindsay F.", "Skol, Andrew D.", "Smith, Kevin S.", "Snyder, Michael", "Stamatoyannopoulos, John", "Tang, Hua", "Wang, Meng", "Carithers, Latarsha J.", "Guan, Ping", "Koester, Susan E.", "Little, A. Roger", "Moore, Helen M.", "Nierras, Concepcion R.", "Rao, Abhi K.", "Vaught, Jimmie B.", "Volpi, Simona"], "journal":"Science", "volume":"369", "pages":"1318--1330", "doi": "10.1126/science.aaz1776", "url": "https://doi.org/10.1126/science.aaz1776"},
            {"id": "HuBMAP", "type": "article", "title":"Human BioMolecular Atlas Program (HuBMAP): 3D Human Reference Atlas construction and usage","year":"2025", "authors": ["B\u00f6rner, Katy", "Blood, Philip D.", "Silverstein, Jonathan C.", "Ruffalo, Matthew", "Satija, Rahul", "Teichmann, Sarah A.", "Pryhuber, Gloria J.", "Misra, Ravi S.", "Purkerson, Jeffrey M.", "Fan, Jean", "Hickey, John W.", "Molla, Gesmira", "Xu, Chuan", "Zhang, Yun", "Weber, Griffin M.", "Jain, Yashvardhan", "Qaurooni, Danial", "Kong, Yongxin", "Abramson, Jakub", "Anderson, David", "Ardlie, Kristin", "Arends, Mark J.", "Aronow, Bruce J.", "Bajema, Rachel", "Baldock, Richard A.", "Barnowski, Ross", "Barwinska, Daria", "Bernard, Amy", "Betancur, David", "Bidanta, Supriya", "Bj\u00f6rklund, Frida", "Bolin, Axel", "Boppana, Avinash", "Boulter, Luke", "Browne, Kristen", "Brusko, Maigan A.", "Burger, Albert", "Campbell-Thompson, Martha", "Cao-Berg, Ivan", "Caron, Anita R.", "Carroll, Megan", "Chadwick, Chrystal", "Chen, Haoran", "Chen, Lu", "de Bono, Bernard", "Deutsch, Gail", "Ding, Song-Lin", "Donahue, Sean", "El-Achkar, Tarek M.", "Eskaros, Adel", "Falo, Louis", "Farrow, Melissa", "Ferkowicz, Michael J.", "Fisher, Stephen A.", "Gee, James C.", "Germain, Ronald N.", "Ginda, Michael", "Ginty, Fiona", "Gitomer, Sarah A.", "Goldstone, Melanie B.", "Gustilo, Katherine S.", "Hagood, James S.", "Halushka, Marc K.", "Haniffa, Muzlifah A.", "Hanna, Peter", "Hardi, Josef", "He, Yongqun Oliver", "Honick, Brendan John", "Houghton, Derek", "Itkin, Maxim", "Jain, Sanjay", "Jardine, Laura", "Jiang, Z. Gordon", "Ju, Yingnan", "Karunamurthy, Arivarasan", "Kelleher, Neil L.", "Kendall, Timothy J.", "Kruse, Angela R. S.", "Laronda, Monica M.", "Laurent, Louise C.", "Laurenti, Elisa", "Lee, Sujin", "Lein, Ed", "Li, Chenran", "Li, Zhuoyan", "Lin, Shin", "Lin, Yiing", "Lindsay, Scott A.", "Longacre, Teri A.", "Lundberg, Emma", "Maier, Libby", "Malhotra, Rajeev", "Martinez Casals, Anna", "Masci, Anna Maria", "Mathews, Clayton E.", "McDonough, Elizabeth", "McLaughlin, James A.", "Menon, Rajasree", "Menon, Vilas", "Miller, Jeremy A.", "Morgan, Richard", "M\u00fcller, Werner", "Murphy, Robert F.", "Musen, Mark A.", "Nakshatri, Harikrishna", "Nawijn, Martijn C.", "Neumann, Elizabeth K.", "Nigra, Debra J.", "O'Neill, Kathleen", "Parast, Mana M.", "Patel, Ushma", "Pei, Liming", "Phatnani, Hemali", "Phillips, Gesina A.", "Pouch, Alison M.", "Powers, Alvin C.", "Puerto, Juan F.", "Puig-Barbe, Aleix", "Quardokus, Ellen M.", "Radtke, Andrea J.", "Rajbhandari, Presha", "Record, Elizabeth G.", "Roberts, Drucilla J.", "Ropelewski, Alexander J.", "Rowe, David", "Ruschman, Nancy L.", "Saunders, Diane C.", "Scheuermann, Richard H.", "Schey, Kevin L.", "Schilling, Birgit", "Schlehlein, Heidi", "Schwenk, Melissa", "Scibek, Robin", "Seifert, Robert P.", "Shirey, Bill", "Shivkumar, Kalyanam", "Siletti, Kimberly", "Simmons, J. Alan", "Singhal, Dhruv", "Snyder, Michael", "Spraggins, Jeffrey M.", "Stanley, Valentina", "Strand, Douglas W.", "Sunshine, Joel C.", "Surrette, Christine", "Suzuki, Ayako", "Tata, Purushothama Rao", "Taylor, Deanne M.", "Theriault, Todd", "Theriault, Tracey", "Thomas, Jerin Easo", "Tsui, Elizabeth L.", "Uranic, Jackie", "Valerius, M. Todd", "Van Valen, David", "Vezina, Chad M.", "Vlachos, Ioannis S.", "Wang, Fusheng", "Wang, Xuefei \u2018Julie'", "Wasserfall, Clive H.", "Welling, Joel S.", "Werlein, Christopher", "Winfree, Seth", "Wright, Devin M.", "Yao, Li", "Yuan, Zhou", "Zhang, Ted", "Bueckle, Andreas", "Herr, Bruce W."], "journal":"Nature Methods", "volume":"22", "pages":"845--860", "doi": "10.1038/s41592-024-02563-5", "url": "https://doi.org/10.1038/s41592-024-02563-5"},
            {"id": "IDG", "type": "article", "title":"DrugCentral 2023 extends human clinical data and integrates veterinary drugs","year":"2022", "authors": ["Avram, Sorin", "Wilson, Thomas B", "Curpan, Ramona", "Halip, Liliana", "Borota, Ana", "Bora, Alina", "Bologa, Cristian\u00a0G", "Holmes, Jayme", "Knockel, Jeffrey", "Yang, Jeremy\u00a0J", "Oprea, Tudor\u00a0I"], "journal":"Nucleic Acids Research", "volume":"51", "pages":"D1276--D1287", "doi": "10.1093/nar/gkac1085", "url": "https://doi.org/10.1093/nar/gkac1085"},
            {"id": "KOMP2", "type": "article", "title":"The International Mouse Phenotyping Consortium: comprehensive knockout phenotyping underpinning the study of human disease","year":"2022", "authors": ["Groza, Tudor", "Gomez, Federico Lopez", "Mashhadi, Hamed Haseli", "Mu\u00f1oz-Fuentes, Violeta", "Gunes, Osman", "Wilson, Robert", "Cacheiro, Pilar", "Frost, Anthony", "Keskivali-Bond, Piia", "Vardal, Bora", "McCoy, Aaron", "Cheng, Tsz Kwan", "Santos, Luis", "Wells, Sara", "Smedley, Damian", "Mallon, Ann-Marie", "Parkinson, Helen"], "journal":"Nucleic Acids Research", "volume":"51", "pages":"D1038--D1045", "doi": "10.1093/nar/gkac972", "url": "https://doi.org/10.1093/nar/gkac972"},
            {"id": "LINCS", "type": "article", "title":"SigCom LINCS: data and metadata search engine for a million gene expression signatures","year":"2022", "authors": ["Evangelista, John Erol", "Clarke, Daniel J B", "Xie, Zhuorui", "Lachmann, Alexander", "Jeon, Minji", "Chen, Kerwin", "Jagodnik, Kathleen\u00a0M", "Jenkins, Sherry L", "Kuleshov, Maxim\u00a0V", "Wojciechowicz, Megan\u00a0L", "Sch\u00fcrer, Stephan\u00a0C", "Medvedovic, Mario", "Ma'ayan, Avi"], "journal":"Nucleic Acids Research", "volume":"50", "pages":"W697--W709", "doi": "10.1093/nar/gkac328", "url": "https://doi.org/10.1093/nar/gkac328"},
            {"id": "MoTrPAC", "type": "article", "title":"Molecular Transducers of Physical Activity Consortium (MoTrPAC): Mapping the Dynamic Responses to Exercise","year":"2020", "authors": ["Sanford, James A.", "Nogiec, Christopher D.", "Lindholm, Malene E.", "Adkins, Joshua N.", "Amar, David", "Dasari, Surendra", "Drugan, Jonelle K.", "Fern\u00e1ndez, Facundo M.", "Radom-Aizik, Shlomit", "Schenk, Simon", "Snyder, Michael P.", "Tracy, Russell P.", "Vanderboom, Patrick", "Trappe, Scott", "Walsh, Martin J.", "Adkins, Joshua N.", "Amar, David", "Dasari, Surendra", "Drugan, Jonelle K.", "Evans, Charles R.", "Fernandez, Facundo M.", "Li, Yafeng", "Lindholm, Malene E.", "Nogiec, Christopher D.", "Radom-Aizik, Shlomit", "Sanford, James A.", "Schenk, Simon", "Snyder, Michael P.", "Tomlinson, Lyl", "Tracy, Russell P.", "Trappe, Scott", "Vanderboom, Patrick", "Walsh, Martin J.", "Lee Alekel, D.", "Bekirov, Iddil", "Boyce, Amanda T.", "Boyington, Josephine", "Fleg, Jerome L.", "Joseph, Lyndon J.O.", "Laughlin, Maren R.", "Maruvada, Padma", "Morris, Stephanie A.", "McGowan, Joan A.", "Nierras, Concepcion", "Pai, Vinay", "Peterson, Charlotte", "Ramos, Ed", "Roary, Mary C.", "Williams, John P.", "Xia, Ashley", "Cornell, Elaine", "Rooney, Jessica", "Miller, Michael E.", "Ambrosius, Walter T.", "Rushing, Scott", "Stowe, Cynthia L.", "Jack Rejeski, W.", "Nicklas, Barbara J.", "Pahor, Marco", "Lu, Ching-ju", "Trappe, Todd", "Chambers, Toby", "Raue, Ulrika", "Lester, Bridget", "Bergman, Bryan C.", "Bessesen, David H.", "Jankowski, Catherine M.", "Kohrt, Wendy M.", "Melanson, Edward L.", "Moreau, Kerrie L.", "Schauer, Irene E.", "Schwartz, Robert S.", "Kraus, William E.", "Slentz, Cris A.", "Huffman, Kim M.", "Johnson, Johanna L.", "Willis, Leslie H.", "Kelly, Leslie", "Houmard, Joseph A.", "Dubis, Gabriel", "Broskey, Nick", "Goodpaster, Bret H.", "Sparks, Lauren M.", "Coen, Paul M.", "Cooper, Dan M.", "Haddad, Fadia", "Rankinen, Tuomo", "Ravussin, Eric", "Johannsen, Neil", "Harris, Melissa", "Jakicic, John M.", "Newman, Anne B.", "Forman, Daniel D.", "Kershaw, Erin", "Rogers, Renee J.", "Nindl, Bradley C.", "Page, Lindsay C.", "Stefanovic-Racic, Maja", "Barr, Susan L.", "Rasmussen, Blake B.", "Moro, Tatiana", "Paddon-Jones, Doug", "Volpi, Elena", "Spratt, Heidi", "Musi, Nicolas", "Espinoza, Sara", "Patel, Darpan", "Serra, Monica", "Gelfond, Jonathan", "Burns, Aisling", "Bamman, Marcas M.", "Buford, Thomas W.", "Cutter, Gary R.", "Bodine, Sue C.", "Esser, Karyn", "Farrar, Rodger P.", "Goodyear, Laurie J.", "Hirshman, Michael F.", "Albertson, Brent G.", "Qian, Wei-Jun", "Piehowski, Paul", "Gritsenko, Marina A.", "Monore, Matthew E.", "Petyuk, Vladislav A.", "McDermott, Jason E.", "Hansen, Joshua N.", "Hutchison, Chelsea", "Moore, Samuel", "Gaul, David A.", "Clish, Clary B.", "Avila-Pacheco, Julian", "Dennis, Courtney", "Kellis, Manolis", "Carr, Steve", "Jean-Beltran, Pierre M.", "Keshishian, Hasmik", "Mani, D.R.", "Clauser, Karl", "Krug, Karsten", "Mundorff, Charlie", "Pearce, Cadence", "Ivanova, Anna A.", "Ortlund, Eric A.", "Maner-Smith, Kristal", "Uppal, Karan", "Zhang, Tiantian", "Sealfon, Stuart C.", "Zaslavsky, Elena", "Nair, Venugopalan", "Li, SiDe", "Jain, Nimisha", "Ge, YongChao", "Sun, Yifei", "Nudelman, German", "Ruf-zamojski, Frederique", "Smith, Gregory", "Pincas, Nhanna", "Rubenstein, Aliza", "Anne Amper, Mary", "Seenarine, Nitish", "Lappalainen, Tuuli", "Lanza, Ian R.", "Sreekumaran Nair, K.", "Klaus, Katherine", "Montgomery, Stephen B.", "Smith, Kevin S.", "Gay, Nicole R.", "Zhao, Bingqing", "Hung, Chia-Jiu", "Zebarjadi, Navid", "Balliu, Brunilda", "Fresard, Laure", "Burant, Charles F.", "Li, Jun Z.", "Kachman, Maureen", "Soni, Tanu", "Raskind, Alexander B.", "Gerszten, Robert", "Robbins, Jeremy", "Ilkayeva, Olga", "Muehlbauer, Michael J.", "Newgard, Christopher B.", "Ashley, Euan A.", "Wheeler, Matthew T.", "Jimenez-Morales, David", "Raja, Archana", "Dalton, Karen P.", "Zhen, Jimmy", "Suk Kim, Young", "Christle, Jeffrey W.", "Marwaha, Shruti", "Chin, Elizabeth T.", "Hershman, Steven G.", "Hastie, Trevor", "Tibshirani, Robert", "Rivas, Manuel A."], "journal":"Cell", "volume":"181", "pages":"1464--1474", "doi": "10.1016/j.cell.2020.06.004", "url": "https://doi.org/10.1016/j.cell.2020.06.004"},
            {"id": "Metabolomics", "type": "article", "title":"Metabolomics Workbench: An international repository for metabolomics data and metadata, metabolite standards, protocols, tutorials and training, and analysis tools","year":"2015", "authors": ["Sud, Manish", "Fahy, Eoin", "Cotter, Dawn", "Azam, Kenan", "Vadivelu, Ilango", "Burant, Charles", "Edison, Arthur", "Fiehn, Oliver", "Higashi, Richard", "Nair, K. Sreekumaran", "Sumner, Susan", "Subramaniam, Shankar"], "journal":"Nucleic Acids Research", "volume":"44", "pages":"D463--D470", "doi": "10.1093/nar/gkv1042", "url": "https://doi.org/10.1093/nar/gkv1042"},
            {"id": "RummaGEO", "type": "article", "title":"RummaGEO: Automatic mining of human and mouse gene sets from GEO","year":"2024", "authors": ["Marino, Giacomo B.", "Clarke, Daniel J.B.", "Lachmann, Alexander", "Deng, Eden Z.", "Ma'ayan, Avi"], "journal":"Patterns", "volume":"5", "pages":"101072", "doi": "10.1016/j.patter.2024.101072", "url": "https://doi.org/10.1016/j.patter.2024.101072"},
            {"id":"Enrichr","type":"article","title":"Gene Set Knowledge Discovery with Enrichr","authors":["Xie, Zhuorui","Bailey, Allison","Kuleshov, Maxim V.","Clarke, Daniel J. B.","Evangelista, John E.","Jenkins, Sherry L.","Lachmann, Alexander","Wojciechowicz, Megan L.","Kropiwnicki, Eryk","Jagodnik, Kathleen M.","Jeon, Minji","Ma'ayan, Avi"],"year":"2021","journal":"Current Protocols","volume":"1","pages":"e90","doi":"10.1002/cpz1.90","url":"https://doi.org/10.1002/cpz1.90"},
            {"id":"GO","type":"article","title":"Gene Ontology: tool for the unification of biology","authors":["Ashburner, Michael","Ball, Catherine A.","Blake, Judith A.","Botstein, David","Butler, Heather","Cherry, J. Michael","Davis, Allan P.","Dolinski, Kara","Dwight, Selina S.","Eppig, Janan T.","Harris, Midori A.","Hill, David P.","Issel-Tarver, Laurie","Kasarskis, Andrew","Lewis, Suzanna","Matese, John C.","Richardson, Joel E.","Ringwald, Martin","Rubin, Gerald M.","Sherlock, Gavin"],"year":"2000","journal":"Nature Genetics","volume":"25","pages":"25--29","doi":"10.1038/75556","url":"https://doi.org/10.1038/75556"},
            {"id":"KEGG","type":"article","title":"KEGG for taxonomy-based analysis of pathways and genomes","authors":["Kanehisa, Minoru","Furumichi, Miho","Sato, Yoko","Kawashima, Masayuki","Ishiguro-Watanabe, Mari"],"year":"2022","journal":"Nucleic Acids Research","volume":"51","pages":"D587--D592","doi":"10.1093/nar/gkac963","url":"https://doi.org/10.1093/nar/gkac963"},
            {"id":"ChEA","type":"article","title":"ChEA3: transcription factor enrichment analysis by orthogonal omics integration","authors":["Keenan, Alexandra B","Torre, Denis","Lachmann, Alexander","Leong, Ariel K","Wojciechowicz, Megan L","Utti, Vivian","Jagodnik, Kathleen M","Kropiwnicki, Eryk","Wang, Zichen","Ma'ayan, Avi"],"year":"2019","journal":"Nucleic Acids Research","volume":"47","pages":"W212--W224","doi":"10.1093/nar/gkz446","url":"https://doi.org/10.1093/nar/gkz446"},
            {"id": "GWAS", "type": "article", "title":"The NHGRI-EBI GWAS Catalog: standards for reusability, sustainability and diversity","year":"2024", "authors": ["Cerezo, Maria", "Sollis, Elliot", "Ji, Yue", "Lewis, Elizabeth", "Abid, Ala", "Bircan, Karatu\u011f\u00a0Ozan", "Hall, Peggy", "Hayhurst, James", "John, Sajo", "Mosaku, Abayomi", "Ramachandran, Santhi", "Foreman, Amy", "Ibrahim, Arwa", "McLaughlin, James", "Pendlington, Zo\u00eb", "Stefancsik, Ray", "Lambert, Samuel A", "McMahon, Aoife", "Morales, Joannella", "Keane, Thomas", "Inouye, Michael", "Parkinson, Helen", "Harris, Laura W"], "journal":"Nucleic Acids Research", "volume":"53", "pages":"D998--D1005", "doi": "10.1093/nar/gkae1070", "url": "https://doi.org/10.1093/nar/gkae1070"},
        ]

        async with ClientSession() as session:
            tasks = [
                doi_to_ref(session, key, doi)
                for summary, gene_refs in deepdive_results.values()
                for key, doi in gene_refs.items()
            ]
            refs.extend(await asyncio.gather(*tasks))

        log("References complete.")
        return refs

    async def stage_sections(self,  deepdive_results):
        introduction, discussion, references = await asyncio.gather(
            
            self.write_report_introduction(),
            self.write_report_discussion(deepdive_results),
            self.make_references(deepdive_results)
        )
        abstract = await self.write_report_abstract()

        title_prompt = '''Write a title for this gene set crossing report using the following structure:

        Crossing Gene Sets from [DATASETS] Reveals [ENTITY] as a Potential Modulator of [PROCESS] and [PROCESS]

        where:

        [DATASETS]: list the source resource names exactly as provided, separated by commas (if 3+) with "and" before the final resource.
        [ENTITY]: select a single intersecting gene, protein, or compound from the provided terms. Choose the entity that
        best represents the shared biological theme of the crossed gene sets. Select which intersecting entity would best serve as
        a marker, perturbagen, or target; do not select entities that are not in the crossing terms.
        [PROCESS]: summarize each crossed term into one or two concise biological concepts (such as a biological process, phenotype,
        disease, cellular function, or signaling pathway). Use broad, recognizable concepts rather than copying full term names.
        Include two or three processes when appropriate. Merge redundant concepts and avoid repeating synonymous ideas.

        Resource descriptions are provided to indicate what each gene set collection represents
         (for example, GlyGen — Glycosylated Proteins; GTEx — Tissue-Specific Aging Signatures).

        Use these descriptions only to infer the underlying biological themes when selecting the process keywords. In the
        generated title, include only the resource names exactly as provided and do not include the descriptive phrases. Derive
        the process keywords only from the crossed term names, using the resource descriptions only as
        supporting context when the biological theme is unclear. Processes should be derived from the dataset descriptions and ge set names.
        Do not include a process if the corresponding term is already selected as the entity (Example: a glycan is selected for an aging/
        exercise/glycosylation corssing, so only aging and exercise are selected as processes).
        

        The title should summarize the common biology implied by the intersecting gene sets rather than simply restating their names.
        Preserve acronyms exactly as written. Do not include statistics, gene set sizes, p-values, Jaccard values,
        parentheses, quotation marks, or other special symbols. Do not invent biological concepts that are unsupported 
        by the provided terms. Return only the title. Use title case, capitalizaing all keywords.'''

        resources = [
            dataset["dataset"]["resource"]
            for dataset in self.datasets
        ]
        dataset_descriptions = [
            dataset["dataset"]["name"]
            .replace(dataset["dataset"]["resource"], f'{dataset["dataset"]["resource"]} - ') 
            for dataset in self.datasets
        ]
        terms = [
            val for col,val in self.crossing.items() if (col.startswith("term") and "Length" not in col)
        ]
        title_data = f'Resources:{resources}\nDataset Descriptions:{dataset_descriptions}\nCrossing Terms:{terms}'

        REPORT_TITLE = await generate_section(self.client, self.model, self.system_prompt, title_prompt, title_data)
        self.title = REPORT_TITLE

        log("Sections complete.")
        return REPORT_TITLE,abstract,introduction,discussion,references

def construct_crossing_report(datasets, gmts, crossing, venn_diagram, intersecting_genes, enrichr):
    client = AsyncOpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    model = os.getenv("OPENAI_MODEL", "gpt-5-nano")

    enrichr_id = enrichr.pop("enrichr_id")["shortId"]
    enrichr_results = {
        library:pd.DataFrame(results["scored"]).head(10) for library,results in enrichr.items()
    }

    crossing_report = CrossingReport(client, model, datasets, gmts, crossing, venn_diagram, intersecting_genes["set"], enrichr_results)

    log("Fetching intersecting gene DeepDives...")
    deepdive_results = {}
    for gene in crossing_report.intersecting_genes:
        summary = get_gpt4o_gene_summary(gene)
        deepdive_results[gene] = format_deepdive(gene,summary)

    log("Generating sections...")
    methods = crossing_report.write_report_methods()
    results = crossing_report.write_report_results(deepdive_results)
    title,abstract,introduction,discussion,references = asyncio.run(crossing_report.stage_sections(deepdive_results))
    

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

    venn_diagram_pdf = make_venn_diagram(venn_diagram)
    venn_diagram_file =  {
        "file":venn_diagram_pdf,
        "caption":f"A venn diagram of the crossed gene sets. Sizes of each set and {'2-,3-, and 4-way' if len(datasets)==4 else '2- and 3-way'} intersection are labelled.",
    }
    enrichr_bars_pdf = make_enrichr_barplot(crossing_report.enrichr)
    enrichr_bars_file =  {
        "file":enrichr_bars_pdf,
        "caption":f"Enrichment analysis of the up-regulated gene signature using Enrichr~\\cite{{Enrichr}}. Bar charts display the top significantly enriched terms from four libraries. A. GO Biological Process 2023~\\cite{{GO}}, B. KEGG 2021 Human~\\cite{{KEGG}}, C. ChEA 2022~\\cite{{ChEA}}, and D. GWAS Catalog 2019~\\cite{{GWAS}}. Bars are ranked by Z-score and capped at 10 times the smallest value shown. Color indicates the source library. The full enrichment results are available to view at \\href{{https://maayanlab.cloud/enrichr/enrich?dataset={enrichr_id}}}{{Enrichr}}.",
    }
    figures = {
        "vennDiagram": venn_diagram_file,
        "enrichrBars": enrichr_bars_file,
    }
    log("Figures complete.")
    return {
        "title": crossing_report.title,
        "abstract": crossing_report.abstract,
        "crossing": crossing_report.crossing,
        "introduction": crossing_report.introduction,
        "methods": crossing_report.methods,
        "results": crossing_report.results,
        "discussion": crossing_report.discussion,
        "figures": figures,
        "supplement": {
            "enrichrId":enrichr_id,
            "datasets":datasets
        },
        "references": references,
        "model": crossing_report.model
    }


def construct_crossing_references(references: list[dict]) -> str:
    entries = []
    for ref in references:
        fields = [
            f"  title     = {{{ref['title']}}}" if ref.get('title') else None,
            f"  author    = {{{' and '.join(ref['authors'])}}}" if ref.get('authors') else None,
            f"  year      = {{{ref['year']}}}"    if ref.get('year')    else None,
            f"  journal   = {{{ref['journal']}}}" if ref.get('journal') else None,
            f"  volume    = {{{ref['volume']}}}"  if ref.get('volume')  else None,
            f"  pages     = {{{ref['pages']}}}"   if ref.get('pages')   else None,
            f"  doi       = {{{ref['doi']}}}"     if ref.get('doi')     else None,
            f"  url       = {{{ref['url']}}}"     if ref.get('url')     else None,
        ]
        fields = [f for f in fields if f is not None]
        entries.append(f"@{ref['type']}{{{ref['id']}, {', '.join(fields)} }}")

    return '\n\n'.join(entries)

def render_crossing_report(report):
    log('Preparing LaTeX bundle...')
    crossing = report["crossing"]
    n = (len(crossing)-5)//2
    datasets = [dataset["key"] for dataset in report["supplement"]["datasets"]]
    crossing_name = f'Crossing_{"".join(datasets)}_{crossing["rank"]}'
    with TemporaryDirectory() as texdir:
        public_url = os.getenv("PUBLIC_URL")
        urllib.request.urlretrieve(f"{public_url}/crossing-report/crossing_report_{n}.tex", f"{texdir}/crossing_report.tex")
        urllib.request.urlretrieve(f"{public_url}/crossing-report/ExcelAtFIT.cls", f"{texdir}/ExcelAtFIT.cls")
        urllib.request.urlretrieve(f"{public_url}/crossing-report/maayan-lab-logo.png", f"{texdir}/maayan-lab-logo.png")
        os.makedirs(f"{texdir}/figures", exist_ok=True)
        for figure,figurefile in report["figures"].items():
            with file_as_path(figurefile["file"]) as path:
                shutil.copy(path, f"{texdir}/figures/{figure}.pdf")
        os.makedirs(f"{texdir}/tables", exist_ok=True)
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(texdir),
            variable_start_string='(((',
            variable_end_string=')))',
        )
        template = env.get_template("crossing_report.tex").render(
            title=report['title'],
            abstract=report['abstract'],
            introduction = report['introduction'],
            methods = report['methods'],
            results = report['results'],
            discussion = report['discussion'],
            figures = report['figures'],
            model = report['model'],
        )
        log('Rendering LaTeX template...')
        with open(f"{texdir}/{crossing_name}.tex", "w") as texfile:
            texfile.write(template)
        with open(f"{texdir}/references.bib", "w") as bibfile:
            bibfile.write(construct_crossing_references(report["references"]))
        with upsert_file('.zip') as zipfile:
            base = zipfile.file.rsplit('.', 1)[0]
            shutil.make_archive(base, "zip", texdir)
        log('Bundle complete.')
        log('Beginning PDF compilation...')

        with subprocess.Popen(
            ['latexmk', '-pdf', '-lualatex', '-cd', '-f', '-interaction=nonstopmode', f"{texdir}/{crossing_name}.tex"],
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            text=True) as proc:
            for line in proc.stdout:
                if line.startswith("Run"):
                    log(line.strip())

        with upsert_file(".pdf") as pdffile:
            reader = pypdf.PdfReader(f'{texdir}/{crossing_name}.pdf')
            writer = pypdf.PdfWriter()
            for page in reader.pages:
                writer.add_page(page)
            writer.write(pdffile.file)

    zipfile["filename"] = f"{crossing_name}.zip"
    pdffile["filename"] = f"{crossing_name}.pdf"

    return {"bundle":zipfile, "pdf":pdffile}