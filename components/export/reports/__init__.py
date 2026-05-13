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
                "content": '''Review the paragraph above against the instructions given earlier in this conversation.
                    Check and correct each of the following:
                        1. All citations use the exact format [[[key]]] with no variations. You must use triple-brackets for formatting purposes.
                        2. No citation keys have been added that were not present in the source material or allowed reference lists.
                        3. The output is plain prose only with no markdown, bullet points, or headings.
                        4. No claims are made that go beyond the scope of the provided data or source material.
                    Return only the corrected paragraph. Do not explain your changes.''',
            },
        )
    )

    response = await client.responses.create(
        model=model,
        input=history

    )

    return format_citations(response.output_text.replace("%","\\%"))
    

async def write_report_abstract(client: AsyncOpenAI, model:str, geo_accession:str, introduction:dict[str,str], methods:str, results:str, discussion:dict[str,str]):
    log("Writing abstract...")
    prompt = '''Write a single plain text paragraph abstract of approximately 200 words for a re-analysis report.

    Cover the following points in order, briefly:
    1. The biological or clinical problem the original study addressed.
    2. The motivation for re-analyzing the data from that study.
    3. What the re-analysis workflow did (include tools used but not specific methods)
    4. One or two specific findings from the discussion that best represent the re-analysis results.

    Do not copy sentences verbatim from the provided sections.
    You MUST NOT include citations or reference keys of any kind.'''
    data = f'''GEO Accession:{geo_accession}
    Introduction:{introduction}
    Methods: {methods}
    Results: {results}
    Discussion: {discussion}'''

    abstract = await generate_section(client, model, prompt, data)
    abstract = re.sub(r"\[\[\[(.+?)\]\]\]", "", abstract)

    log("Abstract complete.")
    return abstract


async def write_report_introduction(client: AsyncOpenAI, model:str, geo_accession:str, pmc_articles:list, labelled_samples:pd.DataFrame):
    log("Writing introduction...")
    pmc_articles = str(pmc_articles)

    problem_prompt = '''Write a plain text paragraph introducing the biological or clinical problem  investigated in the GEO study.
    Use only the titles, abstracts, and introductions of the provided PMC articles as your source.
    Describe the problem being studied and why it is relevant. Do not describe the methods used in the original study.
    Do not mention the re-analysis. 
    Include inline citations for any claims you make using references that appear in the articles.'''

    background_prompt = '''Write a plain text paragraph providing biological background that contextualizes the problem introduced in the GEO study.
    Use only the titles, abstracts, and introductions of the provided PMC articles as your source.
    Focus on established biology, prior work, or relevant context that motivates the study.
    Expand on the previous problem paragraph, but avoid repeating information or abbreviations that have already been introduced.
    Do not describe any analysis methods from the original study or the re-analysis.
    Do not repeat the specific problem statement - assume that has already been introduced.
    Include inline citations for any claims you make using references that appear in the articles.'''

    signature_prompt = '''Create a short name for the signature being re-analyzed by identifying the condition(s)
    that are being used as perturbations in the labelled samples. Return only the name of the signature. For example:
    "metformin signature" or "melanoma cell line signature". Avoid acronyms and be as concise as possible.'''


    introduction_problem, SIGNATURE_NAME= await asyncio.gather(
        generate_section(client, model,problem_prompt, pmc_articles),
        generate_section(client, model, signature_prompt, labelled_samples)
    )
    introduction_background = await generate_section(client, model, background_prompt, f"Articles:{pmc_articles}\nIntroduction Problem:{introduction_problem}")

    introduction_motivation = f'''In order to further investigate any underlying mechanisms or regulatory activity, we performed a re-analysis of samples from {geo_accession} ~\\cite{{{geo_accession}}}
    to create a workflow analyzing the {SIGNATURE_NAME} signature utillizing bioinformatics tools. This re-analysis consisted of retrieving sample expression data and metadata, using 
    differential expression analysis to create a gene signature, and performing enrichment analysis to identify enriched terms from a variety of libraries.'''.replace('\n','').replace('\t','')
    log("Introduction complete.")

    return {
        "problem":introduction_problem, 
        "background":introduction_background,
        "motivation":introduction_motivation
    }


def write_report_methods(geo_accession:str, labelld_samples:pd.DataFrame):
    template = f'''The workflow starts with selecting {geo_accession} as the search term. GEO studies were identified matching {geo_accession} using 
    ARCHS4 ~\\cite{{ARCHS4}} term search. The GEO study accession was used to fetch the linked publication accession from PMC. Gene expression counts and sample 
    metadata for published samples were obtained from ARCHS4 ~\\cite{{ARCHS4}}. An AnnData file was prepared from the input data and metadata ~\\cite{{AnnData}}. Genes from the 
    anndata matrix were filtered to include protein-coding genes. The samples were then labeled as either control or perturbation to allow for 
    further analysis. The AnnData file was then visualized as a bar plot representing library sizes. Dimensionality reduction of the data was performed 
    using PCA with the normalization set to log-counts-per-million (logCPM). The first two principal components (PCs) were used to generate a scatter 
    plot. The AnnData file was then analyzed using differential expression by Limma-Voom ~\\cite{{limma, voom}} to create a gene signature using the selected conditions. 
    The data in the differential expression table was then visualized as a volcano plot. The up-regulated genes were extracted from the gene signature 
    computed by the Limma-Voom analysis from the file. The gene set containing significant up genes was extracted from the gene signatureand submitted to 
    Enrichr ~\\\{{Enrichr}}. The gene set was enriched against the GO Biological Process 2023 ~\\cite{{GO}}, KEGG 2021 Human ~\\cite{{KEGG}}, ChEA 2022 ~\\cite{{ChEA}}, and KOMP2 Mouse Phenotypes 
    2022 ~\\cite{{KOMP2}} libraries to identify statistically significant enriched biological processes, pathways, transcription factors and phenotypes. 
    The gene set containing significant down genes was extracted from the gene signatureand submitted to Enrichr ~\\cite{{Enrichr}}. The gene set was enriched against 
    the GO Biological Process 2023 ~\\cite{{GO}}, KEGG 2021 Human ~\\cite{{KEGG}}, ChEA 2022 ~\\cite{{ChEA}}, and KOMP2 Mouse Phenotypes 2022 ~\\cite{{KOMP2}} libraries to identify statistically 
    significant enriched biological processes, pathways, transcription factors and phenotypes. Significant genes were extracted from the gene signature 
    and submitted to Perturb-Seqr [Perturb-Seqr] to identify small molecules and single gene CRISPR KOs producing gene expression profiles similar or opposite to 
    the signature.'''

    return template.replace('\n','').replace('\t','')


def write_report_results(geo_accession:str, labelld_samples:pd.DataFrame, enrichr_results:dict[str,dict[str,pd.DataFrame]], perturbseqr_results:dict[str,pd.DataFrame]):
    text = f'''
    Library size analysis was also performed to document the number of reads included in each sample. Library sizes of each sample were plotted using a bar plot (Figure \\ref{{fig:Libraries}}).
    The first two principal components (PCs) were used to generate a scatter plot (Figure \\ref{{fig:PCA}}). The plot highlights the control and perturbation samples used and displays unused samples from the study.
    Following differential expression analysis with limma-voom, logFC and p-value scores were used to construct a volcano plot (Figure \\ref{{fig:Volcano}}).
    This plot shows the relative change in expression of up- and down-regulated genes, as well as highlighting specific genes with highly significant scores.
    '''

    return text.replace('\n','').replace('\t','')


async def write_report_discussion(client: AsyncOpenAI, model:str, geo_accession:str, pmc_articles:list, labelled_samples:pd.DataFrame, enrichr_results:dict[str,dict[str,pd.DataFrame]], perturbseqr_results:dict[str,pd.DataFrame]):
    log("Writing discussion...")
    pmc_articles = str(pmc_articles)

    enrichr_prompt = '''Write a plain text paragraph analyzing the Enrichr enrichment results provided for the up-regulated and down-regulated gene sets from this re-analysis.
    Focus on: terms that appear or are consistent across multiple libraries, and any complementary or contrasting patterns between the up and down gene sets.
    Only discuss terms that are present in the results provided.
    Do not introduce terms or biological processes not present in the data.
    Cite only the libraries whose results you are directly referencing, choosing from: 
    [[[GO]]], [[[KEGG]]], [[[ChEA]]], [[[KOMP2]]].'''

    perturbseqr_prompt = '''Write a plain text paragraph analyzing the Perturb-Seqr results provided for the gene signature from this re-analysis.
    Focus on: small molecules or genetic perturbations that appear in each mimicker and reverser result table and what those patterns suggest about the biology of the signature.
    Highlight similarities in perturbations within each result table.
    Only discuss entries that are present in the results provided.
    Do not introduce perturbations or mechanisms not present in the data.
    When you reference a mimicker or reverser signature, you MUST cite the associated resource (found in Dataset field), choosing from:
    [[[CMap]]], [[[CM4AI]]], [[[CREEDS]]], [[[DeepCoverMoA]]], [[[Ginkgo]]], [[[LINCS]]], [[[NIBR]]], 
    [[[PerturbAtlas]]], [[[RummaGEO]]], [[[Replogle]]], [[[SciPlex]]], [[[Tahoe]]].'''

    discussion_enrichr, discussion_perturbseqr = await asyncio.gather(
        generate_section(client, model, enrichr_prompt, str(enrichr_results)),
        generate_section(client, model, perturbseqr_prompt, str(perturbseqr_results)),
    )

    conclusion_prompt = '''Write a plain text concluding paragraph for a re-analysis report.
    Summarize the most notable findings from the Enrichr and Perturb-Seqr analyses provided.
    Focus on findings that are consistent with or directly relevant to the biology of the samples.
    Where a finding is unexpected relative to the sample context, note it briefly.
    Do not propose mechanisms unless they are directly supported by the enrichment results provided.
    Connect the findings back to the research question described in the Enrichr and Perturb-Seqr paragraphs.
    Do not include citations.'''

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
    req = requests.get(
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',
        params={
            "id":{",".join(pmc_set).replace("PMC","")},
            "db":"pmc"
        }
    )
    req.raise_for_status()
    return req.text


def extract_text_with_citations(elem):
    parts = []

    if elem.text:
        parts.append(elem.text)

    for child in elem:
        if child.tag == "xref" and child.attrib.get("ref-type") == "bibr":
            label = child.get("rid")
            parts.append(f"[{label}]")
        else:
            parts.append(extract_text_with_citations(child))

        if child.tail:
            parts.append(child.tail)

    return "".join(parts)


def parse_section(sec):
    sec_title_el = sec.find("title")
    sec_title = "".join(sec_title_el.itertext()).strip() if sec_title_el is not None else None

    content_blocks = []

    for p in sec.findall("./p"):
        text = extract_text_with_citations(p).strip()
        if text:
            content_blocks.append(text)

    for subsec in sec.findall("./sec"):
        content_blocks.extend(str(parse_section(subsec)))

    return [{
        "title": sec_title,
        "text": "".join(content_blocks)
    }]


def parse_pmc_xml(geo_accession:str, pmc_set: set[str]):
    xml_string = retrieve_pmc_articles(pmc_set)
    root = ET.fromstring(xml_string)
    articles = []

    for article in root:
        article_title = article.find('.//article-title').text
        front = article.find('.//front')
        abstract_el = article.find('.//abstract')
        abstract = " ".join(
            "".join(p.itertext()).strip()
            for p in abstract_el.findall(".//p")
        ) if abstract_el is not None else ""

        sections = []
        body = article.find(".//body")
        if body is not None:
            # Top-level bare paragraphs (not inside a sec)
            for child in body:
                if child.tag == "p":
                    if text := extract_text_with_citations(child).strip():
                        sections.append({"title": "", "text": text})

            for sec in body.findall("./sec"):
                sections.extend(parse_section(sec))

        references = []
        authors = []
        contrib_group = front.find(".//contrib-group")
        if contrib_group is not None:
            for name in contrib_group.findall('.//name'):
                surname = name.findtext('surname', '')
                given = name.findtext('given-names', '')
                if surname:
                    authors.append(f"{surname}, {given}")
        author_str = " and ".join(authors)

        doi = front.findtext(".//article-id[@pub-id-type='doi']", "")
        pmid = front.findtext(".//article-id[@pub-id-type='pmid']", "")
        fields = {
            "author": author_str,
            "title": article_title,
            "journal": front.findtext(".//journal-id[@journal-id-type='nlm-ta']", ""),
            "year": front.findtext('.//year', ''),
            "volume": front.findtext('.//volume', ''),
            "pages": f"{front.findtext('.//fpage', '')}--{front.findtext('.//lpage', '')}" if front.find('.//fpage') is not None else "",
            "doi": doi,
            "pmid": pmid,
            "url": f"https://doi.org/{doi}" if doi else ""
        }

        bib_lines = [f"@article{{{geo_accession},"]
        for key, value in fields.items():
            if value and value.strip() not in ("", "--"):
                clean_val = html.unescape(str(value)).replace('&', '\\&')
                bib_lines.append(f"  {key} = {{{clean_val}}},")
        bib_lines.append("}")
        references.append("\n".join(bib_lines))

        for ref in article.findall(".//ref"):
            ref_id = ref.attrib.get('id', 'ref_unknown')
            cit = ref.find('.//element-citation') or ref.find('.//mixed-citation')
            if cit is None:
                continue
            pub_type = cit.attrib.get('publication-type', 'article')
            bib_type = "article" if pub_type == "journal" else "misc"

            authors = []
            person_group = cit.find("person-group[@person-group-type='author']")
            if person_group is not None:
                for name in person_group.findall('name')+person_group.findall('string-name'):
                    surname = name.findtext('surname', '')
                    given = name.findtext('given-names', '')
                    if surname:
                        authors.append(f"{surname}, {given}")
            author_str = " and ".join(authors)

            doi = cit.findtext(".//pub-id[@pub-id-type='doi']", "")
            pmid = cit.findtext(".//pub-id[@pub-id-type='pmid']", "")
            fields = {
                "author": author_str,
                "title": cit.findtext('article-title', '').strip(),
                "journal": cit.findtext('source', ''),
                "year": cit.findtext('year', ''),
                "volume": cit.findtext('volume', ''),
                "pages": f"{cit.findtext('fpage', '')}--{cit.findtext('lpage', '')}" if cit.find('fpage') is not None else "",
                "doi": doi,
                "pmid": pmid,
                "url": f"https://doi.org/{doi}" if doi else ""
            }

            bib_lines = [f"@{bib_type}{{{ref_id},"]
            for key, value in fields.items():
                if value and value.strip() not in ("", "--"):
                    clean_val = html.unescape(str(value)).replace('&', '\\&')
                    bib_lines.append(f"  {key} = {{{clean_val}}},")
            bib_lines.append("}")
            references.append("\n".join(bib_lines))

        articles.append({
            "title": article_title,
            "abstract": abstract,
            "body": sections,
            "references": references
        })

    return articles


def extract_labelled_samples(anndata_file):
    meta_df = anndata_from_file(anndata_file).obs
    meta_df['Type: Control or Perturbation'] = meta_df['Type: Control or Perturbation'].astype(str).replace("nan", pd.NA).copy()
    return meta_df.dropna()


def extract_enrichr_results(up_enrichment, down_enrichment):
    up_result = {
        library:pd.DataFrame(results["scored"]).head(10) for library,results in up_enrichment.items()
    }
    down_result = {
        library:pd.DataFrame(results["scored"]).head(10) for library,results in down_enrichment.items()
    }
    return {"enrichr_up":up_result, "enrichr_down":down_result}


def extract_perturbseqr_results(perturbseqr_ids:dict[str,str], n:int=20):
    up_id = perturbseqr_ids['up_id']
    down_id = perturbseqr_ids['down_id']
    perturbseqr_url = 'https://perturbseqr.maayanlab.cloud/enrichpair/download'
    colnames = ["Dataset", "Perturbation", "Perturbation ID", "Cell Line", "Timepoint", "Concentration", "MoA", "FDA Approved", "Gene Set Size Up", "Gene Set Size Down", "n Overlap", "Odds Ratio", "p-value", "Adjusted p-value"]
    mimic_req = requests.get(
            perturbseqr_url,
            headers={'Accept': 'text/tab-separated-values'},
            params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_mimic","maxTotal":n}
        )
    mimic_req.raise_for_status()
    mimic_df = pd.read_csv(io.BytesIO(mimic_req.content), sep='\t').head(n).drop(columns=["signatureCount", "nReverseOverlap", "oddsRatioReverse", "pvalueReverse", "adjPvalueReverse"]).rename_axis("Rank", inplace=False)
    mimic_df["oddsRatioMimic"] = mimic_df["oddsRatioMimic"].astype(float).map(lambda x: np.format_float_positional(x, precision=5))
    mimic_df["pvalueMimic"] = mimic_df["pvalueMimic"].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
    mimic_df["adjPvalueMimic"] = mimic_df["adjPvalueMimic"].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
    mimic_df.index = mimic_df.index.astype(int) + 1
    mimic_df.columns=colnames
    reverse_req = requests.get(
            perturbseqr_url,
            headers={'Accept': 'text/tab-separated-values'},
            params={"datasetup":up_id,"datasetdown":down_id,"sort":"pvalue_reverse","maxTotal":n}
        )
    reverse_req.raise_for_status()
    reverse_df = pd.read_csv(io.BytesIO(reverse_req.content), sep='\t').head(n).drop(columns=["signatureCount", "nMimicOverlap", "oddsRatioMimic", "pvalueMimic", "adjPvalueMimic"]).rename_axis("Rank", inplace=False)
    reverse_df["oddsRatioReverse"] = reverse_df["oddsRatioReverse"].astype(float).map(lambda x: np.format_float_positional(x, precision=5))
    reverse_df["pvalueReverse"] = reverse_df["pvalueReverse"].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
    reverse_df["adjPvalueReverse"] = reverse_df["adjPvalueReverse"].astype(float).map(lambda x: np.format_float_scientific(x, precision=5))
    reverse_df.index = reverse_df.index.astype(int) + 1
    reverse_df.columns=colnames
    return {"mimic_signatures":mimic_df, "reverse_signatures":reverse_df}


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


    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)

    ax.scatter(x, y, c=mpl_colors, s=20, alpha=0.7, linewidths=0)

    ax.axvline(0, linestyle="--", linewidth=1, alpha=0.4)
    ax.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1, alpha=0.4)

    xpad = (x.max() - x.min()) * 0.05
    ypad = (y.max() - y.min()) * 0.08
    ax.set_xlim(x.min() - xpad, x.max() + xpad)
    ax.set_ylim(y.min() - ypad, y.max() + ypad)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

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


def make_figures(plots, enrichr_results:dict[str,dict[str,pd.DataFrame]]):
    return {
        "librarySizes": make_library_size_barplot(plots["library_sizes_plot"]),
        "PCAScatter": make_pca_scatter(plots["pca_plot"]),
        "volcanoScatter": make_volcano_scatter(plots["volcano_plot"]),
        "upEnrichrBars": make_enrichr_barplot(enrichr_results["enrichr_up"]),
        "downEnrichrBars": make_enrichr_barplot(enrichr_results["enrichr_down"])
    }


# Table utility functions
def make_table(table_data:pd.DataFrame):
    with upsert_file('.csv') as f:
        table_data.to_csv(f.file)
    return f


def make_tables(perturbseqr_results):
    return {
        "perturbseqrMimickers": make_table(perturbseqr_results["mimic_signatures"]),
        "perturbseqrReversers": make_table(perturbseqr_results["reverse_signatures"])
    }


def make_references(pmc_articles):
    refs = '''@article{AnnData, title={anndata: Annotated data}, url={https://doi.org/10.1101/2021.12.16.473007}, doi={10.1101/2021.12.16.473007}, abstractNote={<jats:title>Summary</jats:title><jats:p>anndata is a Python package for handling annotated data matrices in memory and on disk (<jats:ext-link xmlns:xlink="http://www.w3.org/1999/xlink" ext-link-type="uri" xlink:href="http://github.com/theislab/anndata">github.com/theislab/anndata</jats:ext-link>), positioned between pandas and xarray. anndata offers a broad range of computationally efficient features including, among others, sparse data support, lazy operations, and a PyTorch interface.</jats:p><jats:sec><jats:title>Statement of need</jats:title><jats:p>Generating insight from high-dimensional data matrices typically works through training models that annotate observations and variables via low-dimensional representations. In exploratory data analysis, this involves<jats:italic>iterative</jats:italic>training and analysis using original and learned annotations and task-associated representations. anndata offers a canonical data structure for book-keeping these, which is neither addressed by pandas (McKinney, 2010), nor xarray (Hoyer &amp; Hamman, 2017), nor commonly-used modeling packages like scikit-learn (Pedregosa et al., 2011).</jats:p></jats:sec>}, publisher={Cold Spring Harbor Laboratory}, author={Virshup, Isaac and Rybakov, Sergei and Theis, Fabian J. and Angerer, Philipp and Wolf, F. Alexander}, year={2021}, month=dec }

@article{ARCHS4, title={Massive mining of publicly available RNA-seq data from human and mouse}, volume={9}, url={https://doi.org/10.1038/s41467-018-03751-6}, doi={10.1038/s41467-018-03751-6}, abstractNote={<jats:title>Abstract</jats:title><jats:p>RNA sequencing (RNA-seq) is the leading technology for genome-wide transcript quantification. However, publicly available RNA-seq data is currently provided mostly in raw form, a significant barrier for global and integrative retrospective analyses. ARCHS4 is a web resource that makes the majority of published RNA-seq data from human and mouse available at the gene and transcript levels. For developing ARCHS4, available FASTQ files from RNA-seq experiments from the Gene Expression Omnibus (GEO) were aligned using a cloud-based infrastructure. In total 187,946 samples are accessible through ARCHS4 with 103,083 mouse and 84,863 human. Additionally, the ARCHS4 web interface provides intuitive exploration of the processed data through querying tools, interactive visualization, and gene pages that provide average expression across cell lines and tissues, top co-expressed genes for each gene, and predicted biological functions and protein–protein interactions for each gene based on prior knowledge combined with co-expression.</jats:p>}, number={1}, journal={Nature Communications}, publisher={Springer Science and Business Media LLC}, author={Lachmann, Alexander and Torre, Denis and Keenan, Alexandra B. and Jagodnik, Kathleen M. and Lee, Hoyjin J. and Wang, Lily and Silverstein, Moshe C. and Ma’ayan, Avi}, year={2018}, month=apr, language={en} }

@article{voom, title={voom: precision weights unlock linear model analysis tools for RNA-seq read counts}, volume={15}, url={https://doi.org/10.1186/gb-2014-15-2-r29}, doi={10.1186/gb-2014-15-2-r29}, number={2}, journal={Genome Biology}, publisher={Springer Science and Business Media LLC}, author={Law, Charity W and Chen, Yunshun and Shi, Wei and Smyth, Gordon K}, year={2014}, pages={R29}, language={en} }

@article{limma, title={limma powers differential expression analyses for RNA-sequencing and microarray studies}, volume={43}, url={https://doi.org/10.1093/nar/gkv007}, doi={10.1093/nar/gkv007}, number={7}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Ritchie, Matthew E. and Phipson, Belinda and Wu, Di and Hu, Yifang and Law, Charity W. and Shi, Wei and Smyth, Gordon K.}, year={2015}, month=jan, pages={e47–e47}, language={en} }

@article{DEA, title={Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data}, volume={14}, url={https://doi.org/10.1186/gb-2013-14-9-r95}, doi={10.1186/gb-2013-14-9-r95}, number={9}, journal={Genome Biology}, publisher={Springer Science and Business Media LLC}, author={Rapaport, Franck and Khanin, Raya and Liang, Yupu and Pirun, Mono and Krek, Azra and Zumbo, Paul and Mason, Christopher E and Socci, Nicholas D and Betel, Doron}, year={2013}, pages={R95}, language={en} }

@article{Enrichr, title={Gene Set Knowledge Discovery with Enrichr}, volume={1}, url={https://doi.org/10.1002/cpz1.90}, doi={10.1002/cpz1.90}, abstractNote={<jats:title>Abstract</jats:title><jats:p>Profiling samples from patients, tissues, and cells with genomics, transcriptomics, epigenomics, proteomics, and metabolomics ultimately produces lists of genes and proteins that need to be further analyzed and integrated in the context of known biology. Enrichr (Chen et al., 2013; Kuleshov et al., 2016) is a gene set search engine that enables the querying of hundreds of thousands of annotated gene sets. Enrichr uniquely integrates knowledge from many high‐profile projects to provide synthesized information about mammalian genes and gene sets. The platform provides various methods to compute gene set enrichment, and the results are visualized in several interactive ways. This protocol provides a summary of the key features of Enrichr, which include using Enrichr programmatically and embedding an Enrichr button on any website. © 2021 Wiley Periodicals LLC.</jats:p><jats:p><jats:bold>Basic Protocol 1</jats:bold>: Analyzing lists of differentially expressed genes from transcriptomics, proteomics and phosphoproteomics, GWAS studies, or other experimental studies</jats:p><jats:p><jats:bold>Basic Protocol 2</jats:bold>: Searching Enrichr by a single gene or key search term</jats:p><jats:p><jats:bold>Basic Protocol 3</jats:bold>: Preparing raw or processed RNA‐seq data through BioJupies in preparation for Enrichr analysis</jats:p><jats:p><jats:bold>Basic Protocol 4</jats:bold>: Analyzing gene sets for model organisms using modEnrichr</jats:p><jats:p><jats:bold>Basic Protocol 5</jats:bold>: Using Enrichr in Geneshot</jats:p><jats:p><jats:bold>Basic Protocol 6</jats:bold>: Using Enrichr in ARCHS4</jats:p><jats:p><jats:bold>Basic Protocol 7</jats:bold>: Using the enrichment analysis visualization Appyter to visualize Enrichr results</jats:p><jats:p><jats:bold>Basic Protocol 8</jats:bold>: Using the Enrichr API</jats:p><jats:p><jats:bold>Basic Protocol 9</jats:bold>: Adding an Enrichr button to a website</jats:p>}, number={3}, journal={Current Protocols}, publisher={Wiley}, author={Xie, Zhuorui and Bailey, Allison and Kuleshov, Maxim V. and Clarke, Daniel J. B. and Evangelista, John E. and Jenkins, Sherry L. and Lachmann, Alexander and Wojciechowicz, Megan L. and Kropiwnicki, Eryk and Jagodnik, Kathleen M. and Jeon, Minji and Ma'ayan, Avi}, year={2021}, month=mar, language={en} }

@article{GO, title={Gene Ontology: tool for the unification of biology}, volume={25}, url={https://doi.org/10.1038/75556}, doi={10.1038/75556}, number={1}, journal={Nature Genetics}, publisher={Springer Science and Business Media LLC}, author={Ashburner, Michael and Ball, Catherine A. and Blake, Judith A. and Botstein, David and Butler, Heather and Cherry, J. Michael and Davis, Allan P. and Dolinski, Kara and Dwight, Selina S. and Eppig, Janan T. and Harris, Midori A. and Hill, David P. and Issel-Tarver, Laurie and Kasarskis, Andrew and Lewis, Suzanna and Matese, John C. and Richardson, Joel E. and Ringwald, Martin and Rubin, Gerald M. and Sherlock, Gavin}, year={2000}, month=may, pages={25-29}, language={en} }

@article{KEGG, title={KEGG for taxonomy-based analysis of pathways and genomes}, volume={51}, url={https://doi.org/10.1093/nar/gkac963}, doi={10.1093/nar/gkac963}, abstractNote={<jats:title>Abstract</jats:title>
               <jats:p>KEGG (https://www.kegg.jp) is a manually curated database resource integrating various biological objects categorized into systems, genomic, chemical and health information. Each object (database entry) is identified by the KEGG identifier (kid), which generally takes the form of a prefix followed by a five-digit number, and can be retrieved by appending /entry/kid in the URL. The KEGG pathway map viewer, the Brite hierarchy viewer and the newly released KEGG genome browser can be launched by appending /pathway/kid, /brite/kid and /genome/kid, respectively, in the URL. Together with an improved annotation procedure for KO (KEGG Orthology) assignment, an increasing number of eukaryotic genomes have been included in KEGG for better representation of organisms in the taxonomic tree. Multiple taxonomy files are generated for classification of KEGG organisms and viruses, and the Brite hierarchy viewer is used for taxonomy mapping, a variant of Brite mapping in the new KEGG Mapper suite. The taxonomy mapping enables analysis of, for example, how functional links of genes in the pathway and physical links of genes on the chromosome are conserved among organism groups.</jats:p>}, number={D1}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Kanehisa, Minoru and Furumichi, Miho and Sato, Yoko and Kawashima, Masayuki and Ishiguro-Watanabe, Mari}, year={2022}, month=oct, pages={D587–D592}, language={en} }

@article{ChEA, title={ChEA3: transcription factor enrichment analysis by orthogonal omics integration}, volume={47}, url={https://doi.org/10.1093/nar/gkz446}, doi={10.1093/nar/gkz446}, abstractNote={<jats:title>Abstract</jats:title><jats:p>Identifying the transcription factors (TFs) responsible for observed changes in gene expression is an important step in understanding gene regulatory networks. ChIP-X Enrichment Analysis 3 (ChEA3) is a transcription factor enrichment analysis tool that ranks TFs associated with user-submitted gene sets. The ChEA3 background database contains a collection of gene set libraries generated from multiple sources including TF–gene co-expression from RNA-seq studies, TF-target associations from ChIP-seq experiments, and TF-gene co-occurrence computed from crowd-submitted gene lists. Enrichment results from these distinct sources are integrated to generate a composite rank that improves the prediction of the correct upstream TF compared to ranks produced by individual libraries. We compare ChEA3 with existing TF prediction tools and show that ChEA3 performs better. By integrating the ChEA3 libraries, we illuminate general transcription factor properties such as whether the TF behaves as an activator or a repressor. The ChEA3 web-server is available from https://amp.pharm.mssm.edu/ChEA3.</jats:p>}, number={W1}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Keenan, Alexandra B and Torre, Denis and Lachmann, Alexander and Leong, Ariel K and Wojciechowicz, Megan L and Utti, Vivian and Jagodnik, Kathleen M and Kropiwnicki, Eryk and Wang, Zichen and Ma'ayan, Avi}, year={2019}, month=may, pages={W212-W224}, language={en} }

@article{KOMP2, title={The International Mouse Phenotyping Consortium: comprehensive knockout phenotyping underpinning the study of human disease}, volume={51}, url={https://doi.org/10.1093/nar/gkac972}, doi={10.1093/nar/gkac972}, abstractNote={<jats:title>Abstract</jats:title>
               <jats:p>The International Mouse Phenotyping Consortium (IMPC; https://www.mousephenotype.org/) web portal makes available curated, integrated and analysed knockout mouse phenotyping data generated by the IMPC project consisting of 85M data points and over 95,000 statistically significant phenotype hits mapped to human diseases. The IMPC portal delivers a substantial reference dataset that supports the enrichment of various domain-specific projects and databases, as well as the wider research and clinical community, where the IMPC genotype–phenotype knowledge contributes to the molecular diagnosis of patients affected by rare disorders. Data from 9,000 mouse lines and 750 000 images provides vital resources enabling the interpretation of the ignorome, and advancing our knowledge on mammalian gene function and the mechanisms underlying phenotypes associated with human diseases. The resource is widely integrated and the lines have been used in over 4,600 publications indicating the value of the data and the materials.</jats:p>}, number={D1}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Groza, Tudor and Gomez, Federico Lopez and Mashhadi, Hamed Haseli and Muñoz-Fuentes, Violeta and Gunes, Osman and Wilson, Robert and Cacheiro, Pilar and Frost, Anthony and Keskivali-Bond, Piia and Vardal, Bora and McCoy, Aaron and Cheng, Tsz Kwan and Santos, Luis and Wells, Sara and Smedley, Damian and Mallon, Ann-Marie and Parkinson, Helen}, year={2022}, month=oct, pages={D1038–D1045}, language={en} }

@article{CMap, title={The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease}, volume={313}, ISSN={1095-9203}, url={https://doi.org/10.1126/science.1132939}, doi={10.1126/science.1132939}, number={5795}, journal={Science}, publisher={American Association for the Advancement of Science (AAAS)}, author={Lamb, Justin and Crawford, Emily D. and Peck, David and Modell, Joshua W. and Blat, Irene C. and Wrobel, Matthew J. and Lerner, Jim and Brunet, Jean-Philippe and Subramanian, Aravind and Ross, Kenneth N. and Reich, Michael and Hieronymus, Haley and Wei, Guo and Armstrong, Scott A. and Haggarty, Stephen J. and Clemons, Paul A. and Wei, Ru and Carr, Steven A. and Lander, Eric S. and Golub, Todd R.}, year={2006}, month=sep, pages={1929-1935} }

@article{CM4AI, title={A Perturbation Cell Atlas of Human Induced Pluripotent Stem Cells}, url={https://doi.org/10.1101/2024.11.03.621734}, doi={10.1101/2024.11.03.621734}, publisher={openRxiv}, author={Nourreddine, Sami and Doctor, Yesh and Dailamy, Amir and Forget, Antoine and Lee, Yi-Hung and Chinn, Becky and Khaliq, Hammza and Polacco, Benjamin and Muralidharan, Monita and Pan, Emily and Zhang, Yifan and Sigaeva, Alina and Hansen, Jan Niklas and Gao, Jiahao and Parker, Jillian A. and Obernier, Kirsten and Clark, Timothy and Chen, Jake Y. and Metallo, Christian and Lundberg, Emma and Ideker, Trey and Krogan, Nevan and Mali, Prashant}, year={2024}, month=nov }

@article{CREEDS, title={Extraction and analysis of signatures from the Gene Expression Omnibus by the crowd}, volume={7}, ISSN={2041-1723}, url={https://doi.org/10.1038/ncomms12846}, doi={10.1038/ncomms12846}, number={1}, journal={Nature Communications}, publisher={Springer Science and Business Media LLC}, author={Wang, Zichen and Monteiro, Caroline D. and Jagodnik, Kathleen M. and Fernandez, Nicolas F. and Gundersen, Gregory W. and Rouillard, Andrew D. and Jenkins, Sherry L. and Feldmann, Axel S. and Hu, Kevin S. and McDermott, Michael G. and Duan, Qiaonan and Clark, Neil R. and Jones, Matthew R. and Kou, Yan and Goff, Troy and Woodland, Holly and Amaral, Fabio M R. and Szeto, Gregory L. and Fuchs, Oliver and Schüssler-Fiorenza Rose, Sophia M. and Sharma, Shvetank and Schwartz, Uwe and Bausela, Xabier Bengoetxea and Szymkiewicz, Maciej and Maroulis, Vasileios and Salykin, Anton and Barra, Carolina M. and Kruth, Candice D. and Bongio, Nicholas J. and Mathur, Vaibhav and Todoric, Radmila D and Rubin, Udi E. and Malatras, Apostolos and Fulp, Carl T. and Galindo, John A. and Motiejunaite, Ruta and Jüschke, Christoph and Dishuck, Philip C. and Lahl, Katharina and Jafari, Mohieddin and Aibar, Sara and Zaravinos, Apostolos and Steenhuizen, Linda H. and Allison, Lindsey R. and Gamallo, Pablo and de Andres Segura, Fernando and Dae Devlin, Tyler and Pérez-García, Vicente and Ma’ayan, Avi}, year={2016}, month=sep }

@article{DeepCoverMOA, title={A proteome-wide atlas of drug mechanism of action}, volume={41}, ISSN={1546-1696}, url={https://doi.org/10.1038/s41587-022-01539-0}, doi={10.1038/s41587-022-01539-0}, number={6}, journal={Nature Biotechnology}, publisher={Springer Science and Business Media LLC}, author={Mitchell, Dylan C. and Kuljanin, Miljan and Li, Jiaming and Van Vranken, Jonathan G. and Bulloch, Nathan and Schweppe, Devin K. and Huttlin, Edward L. and Gygi, Steven P.}, year={2023}, month=jan, pages={845–857} }

@article{Ginkgo, title={Ginkgo Bioworks}, url={https://www.ginkgo.bio/} }

@article{LINCS, title={SigCom LINCS: data and metadata search engine for a million gene expression signatures}, volume={50}, ISSN={1362-4962}, url={https://doi.org/10.1093/nar/gkac328}, doi={10.1093/nar/gkac328}, number={W1}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Evangelista, John Erol and Clarke, Daniel J B and Xie, Zhuorui and Lachmann, Alexander and Jeon, Minji and Chen, Kerwin and Jagodnik, Kathleen M and Jenkins, Sherry L and Kuleshov, Maxim V and Wojciechowicz, Megan L and Schürer, Stephan C and Medvedovic, Mario and Ma’ayan, Avi}, year={2022}, month=may, pages={W697-W709} }

@article{NIBR, title={DRUG-seq Provides Unbiased Biological Activity Readouts for Neuroscience Drug Discovery}, volume={17}, ISSN={1554-8937}, url={https://doi.org/10.1021/acschembio.1c00920}, doi={10.1021/acschembio.1c00920}, number={6}, journal={ACS Chemical Biology}, publisher={American Chemical Society (ACS)}, author={Li, Jingyao and Ho, Daniel J. and Henault, Martin and Yang, Chian and Neri, Marilisa and Ge, Robin and Renner, Steffen and Mansur, Leandra and Lindeman, Alicia and Kelly, Brian and Tumkaya, Tayfun and Ke, Xiaoling and Soler-Llavina, Gilberto and Shanker, Gopi and Russ, Carsten and Hild, Marc and Gubser Keller, Caroline and Jenkins, Jeremy L. and Worringer, Kathleen A. and Sigoillot, Frederic D. and Ihry, Robert J.}, year={2022}, month=may, pages={1401–1414} }

@article{PerturbAtlas, title={PerturbAtlas: a comprehensive atlas of public genetic perturbation bulk RNA-seq datasets}, volume={53}, ISSN={1362-4962}, url={https://doi.org/10.1093/nar/gkae851}, doi={10.1093/nar/gkae851}, number={D1}, journal={Nucleic Acids Research}, publisher={Oxford University Press (OUP)}, author={Zhang, Yiming and Zhang, Ting and Yang, Gaoxia and Pan, Zhenzhong and Tang, Min and Wen, Yue and He, Ping and Wang, Yuan and Zhou, Ran}, year={2024}, month=oct, pages={D1112–D1119} }

@online{Perturb-Seqr, title={Perturb-Seqr}, year=2026, url={https://perturbseqr.maayanlab.cloud}} }

@article{Replogle, title={Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq}, volume={185}, ISSN={0092-8674}, url={https://doi.org/10.1016/j.cell.2022.05.013}, doi={10.1016/j.cell.2022.05.013}, number={14}, journal={Cell}, publisher={Elsevier BV}, author={Replogle, Joseph M. and Saunders, Reuben A. and Pogson, Angela N. and Hussmann, Jeffrey A. and Lenail, Alexander and Guna, Alina and Mascibroda, Lauren and Wagner, Eric J. and Adelman, Karen and Lithwick-Yanai, Gila and Iremadze, Nika and Oberstrass, Florian and Lipson, Doron and Bonnar, Jessica L. and Jost, Marco and Norman, Thomas M. and Weissman, Jonathan S.}, year={2022}, month=july, pages={2559-2575.e28} }

@article{RummaGEO, title={RummaGEO: Automatic mining of human and mouse gene sets from GEO}, volume={5}, ISSN={2666-3899}, url={https://doi.org/10.1016/j.patter.2024.101072}, doi={10.1016/j.patter.2024.101072}, number={10}, journal={Patterns}, publisher={Elsevier BV}, author={Marino, Giacomo B. and Clarke, Daniel J.B. and Lachmann, Alexander and Deng, Eden Z. and Ma'ayan, Avi}, year={2024}, month=oct, pages={101072} }

@article{SciPlex, title={Massively multiplex chemical transcriptomics at single-cell resolution}, volume={367}, ISSN={1095-9203}, url={https://doi.org/10.1126/science.aax6234}, doi={10.1126/science.aax6234}, number={6473}, journal={Science}, publisher={American Association for the Advancement of Science (AAAS)}, author={Srivatsan, Sanjay R. and McFaline-Figueroa, José L. and Ramani, Vijay and Saunders, Lauren and Cao, Junyue and Packer, Jonathan and Pliner, Hannah A. and Jackson, Dana L. and Daza, Riza M. and Christiansen, Lena and Zhang, Fan and Steemers, Frank and Shendure, Jay and Trapnell, Cole}, year={2020}, month=jan, pages={45-51} }

@article{Tahoe, title={Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling}, url={https://doi.org/10.1101/2025.02.20.639398}, doi={10.1101/2025.02.20.639398}, publisher={openRxiv}, author={Zhang, Jesse and Ubas, Airol A and de Borja, Richard and Svensson, Valentine and Thomas, Nicole and Thakar, Neha and Lai, Ian and Winters, Aidan and Khan, Umair and Jones, Matthew G. and Thompson, John D. and Tran, Vuong and Pangallo, Joseph and Papalexi, Efthymia and Sapre, Ajay and Nguyen, Hoai and Sanderson, Oliver and Nigos, Maria and Kaplan, Olivia and Schroeder, Sarah and Hariadi, Bryan and Marrujo, Simone and Salvino, Crina Curca Alec and Gallareta Olivares, Guillermo and Koehler, Ryan and Geiss, Gary and Rosenberg, Alexander and Roco, Charles and Merico, Daniele and Alidoust, Nima and Goodarzi, Hani and Yu, Johnny}, year={2025}, month=feb }'''.split('\n\n')

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


def construct_georeanalysis_report(geo_accession, pmc_set, labelled_samples_anndata, plots, enrichr_up, enrichr_down, perturbseqr):
    client = AsyncOpenAI(api_key=os.getenv("OPENAI_API_KEY"))
    model = os.getenv("OPENAI_MODEL", "gpt-5-nano")
    labelled_samples = extract_labelled_samples(labelled_samples_anndata).to_json()
    enrichr_results = extract_enrichr_results(enrichr_up,enrichr_down)
    perturbseqr_results = extract_perturbseqr_results(perturbseqr)
    pmc_articles = parse_pmc_xml(geo_accession,pmc_set)
    log("PMC articles retrieved.")
    
    log("Generating sections...")
    methods = write_report_methods(geo_accession, labelled_samples)
    results = write_report_results(geo_accession, labelled_samples, enrichr_results, perturbseqr_results)
    title,abstract,introduction,discussion = asyncio.run(stage_sections(client, model, geo_accession, pmc_articles, labelled_samples, enrichr_results, perturbseqr_results, methods, results))
    log("Sections complete.")

    log("Collecting references...")
    references = make_references(pmc_articles)
    log("References complete.")
    
    figures = make_figures(plots, enrichr_results)
    log("Figures complete.")
    tables = make_tables(perturbseqr_results)

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
        "references":references
    }


