import React from "react";
import { MetaNode } from "@/spec/metanode";
import {
  GeneInfo,
  GeneInfoFromGeneTerm,
} from "@/components/service/mygeneinfo";
import { z } from "zod";
import { glygen_icon } from "@/icons";
import { ProteinTerm, GeneTerm } from "@/components/core/term";
import { ProteinSet } from "@/components/core/set";
import {
  GlyGenProteinResponse,
  GlyGenProteinSetResponse,
  GlyGenProteinResponseArray,
  GlycosylationData,
  PhosphorylationData,
} from "./data_models";
import {
  get_single_protein_data,
  get_protein_set_data,
  extract_specific_data,
  glycosylation_check,
  phosphorylation_check,
} from "./sup_functions";
import { GlycosylationTable, PhosphorylationTable } from "./sup_components";

// --- DATA METANODES --- //

// Individual protein data metanode
export const GlyGenProteinResponseNode = MetaNode("GlyGenProteinResponse")
  .meta({
    label: "GlyGen Protein Products",
    description: "Protein product records in GlyGen",
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinResponse)
  .view((data) => {
    const glyGenLink = `http://www.glygen.org/protein/${data.uniprot.uniprot_canonical_ac}`;

    return (
      <div className="prose max-w-none">
        <div>
          Gene Name: <b>{data.gene.name}</b>
        </div>
        <div>
          <span>UniProtKB Accession: </span>
          <b>
            <a
              href={glyGenLink}
              target="_blank"
              rel="noopener noreferrer"
              style={{ color: "blue" }}
            >
              <u style={{ color: "blue" }}>
                {data.uniprot.uniprot_canonical_ac}
              </u>
            </a>
          </b>
        </div>
        <div>
          Gene location: Chromosome: {data.gene.locus.chromosome} (
          {data.gene.locus.start_pos} - {data.gene.locus.end_pos}, '
          {data.gene.locus.strand}' strand)
        </div>
        <div>UniProtKB ID: {data.uniprot.uniprot_id}</div>
        <div>
          Protein Length: <b>{data.uniprot.length}</b>
        </div>
        <div>UniProtKB Protein Name(s): {data.protein_names.name}</div>
        <div>
          Organism:{" "}
          <b>
            {data.species.name} ({data.species.common_name}; TaxID:{" "}
            {data.species.taxid})
          </b>
        </div>
        <div>
          Phosphoprotein:{" "}
          {data.phosphorylation.phosphorylation ? "True" : "False"}
        </div>
        <div>
          Glycoprotein: {data.glycoprotein.glycosylation ? "True" : "False"}
        </div>
        <br />
        <div>
          {data.glycoprotein.glycosylation &&
            data.glycoprotein.glycosylation_data &&
            data.glycoprotein.glycosylation_data.length > 0 && (
              <GlycosylationTable
                glycosylation_data={
                  data.glycoprotein.glycosylation_data.length > 5
                    ? data.glycoprotein.glycosylation_data.slice(0, 5)
                    : data.glycoprotein.glycosylation_data
                }
                is_preview={data.glycoprotein.glycosylation_data.length > 5}
              />
            )}
        </div>
        <div>
          {data.phosphorylation.phosphorylation &&
            data.phosphorylation.phosphorylation_data &&
            data.phosphorylation.phosphorylation_data.length > 0 && (
              <PhosphorylationTable
                phosphorylation_data={
                  data.phosphorylation.phosphorylation_data.length > 5
                    ? data.phosphorylation.phosphorylation_data.slice(0, 5)
                    : data.phosphorylation.phosphorylation_data
                }
                is_preview={
                  data.phosphorylation.phosphorylation_data.length > 5
                }
              />
            )}
        </div>
      </div>
    );
  })
  .build();

// Protein set data metanode
export const GlyGenProteinSetResponseNode = MetaNode("GlyGenProteinSetResponse")
  .meta({
    label: "GlyGen Protein Products",
    description: "Protein product records in GlyGen",
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinSetResponse)
  .view((data) => {
    return (
      <div className="prose max-w-none">
        <table>
          <thead>
            <tr>
              <th>Gene name</th>
              <th>Uniprot Accession</th>
              <th>Protein Name</th>
              <th>Species</th>
              <th>Glycosylation</th>
              <th>Phosphorylation</th>
              <th>SNV</th>
            </tr>
          </thead>
          <tbody>
            {data.map((entry, index) => (
              <tr key={index}>
                <td>{entry.gene.name}</td>
                <td>
                  <a
                    href={`http://www.glygen.org/protein/${entry.uniprot.uniprot_canonical_ac}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    style={{ color: "blue" }}
                  >
                    <u style={{ color: "blue" }}>
                      {entry.uniprot.uniprot_canonical_ac}
                    </u>
                  </a>
                </td>
                <td>{entry.protein_names.name}</td>
                <td>{`${entry.species.name} (TaxID: ${entry.species.taxid})`}</td>
                <td>
                  {entry.bools.total_n_glycosites +
                    entry.bools.total_o_glycosites >
                  0
                    ? "Yes"
                    : "No"}
                </td>
                <td>{entry.bools.reported_phosphosites > 0 ? "Yes" : "No"}</td>
                <td>{entry.bools.reported_snv > 0 ? "Yes" : "No"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    );
  })
  .build();

// Glycosylation data metanode to display the glycosylation table
export const GlycosylationViewResponseNode = MetaNode(
  "GlycosylationViewResponse",
)
  .meta({
    label: "Glycosylation Information for Glycoproteins",
    description: "Glycosylation product records in GlyGen",
    icon: [glygen_icon],
  })
  .codec(GlycosylationData)
  .view((data) => {
    const emptyGlytoucanAcCount =
      data.glycosylation_data?.filter((entry) => entry.glytoucan_ac === "")
        .length ?? 0;
    const glyGenLink = `http://www.glygen.org/protein/${data.protein_accession}`;

    if ((data.glycosylation_data?.length ?? 0) > 0) {
      return (
        <div className="prose max-w-none">
          <div>
            <span>UniProtKB Accession: </span>
            <b>
              <a
                href={glyGenLink}
                target="_blank"
                rel="noopener noreferrer"
                style={{ color: "blue" }}
              >
                {" "}
                <u style={{ color: "blue" }}>{data.protein_accession}</u>
              </a>
            </b>
          </div>
          <div>
            Total Records: <b>{data.glycosylation_data?.length ?? 0}</b>
          </div>
          <div>
            Records Without GlyTouCan Accessions: <b>{emptyGlytoucanAcCount}</b>
          </div>
          <br />
          <GlycosylationTable
            glycosylation_data={data.glycosylation_data ?? []}
          />
        </div>
      );
    } else {
      return (
        <div className="prose max-w-none">
          <div>
            UniProtKB Accession:
            <b>
              <a
                href={glyGenLink}
                target="_blank"
                rel="noopener noreferrer"
                style={{ color: "blue" }}
              >
                {" "}
                <u style={{ color: "blue" }}>{data.protein_accession}</u>
              </a>
            </b>
          </div>
          <div>
            <b>No Glycosylation Information to Display</b>
          </div>
        </div>
      );
    }
  })
  .build();

export const PhosphorylationViewResponseNode = MetaNode(
  "PhosphorylationViewResponseNode",
)
  .meta({
    label: "Phosphorylation Information for the Phosphoprotein",
    description: "Phosphorylation product records in GlyGen",
    icon: [glygen_icon],
  })
  .codec(PhosphorylationData)
  .view((data) => {
    if ((data.phosphorylation_data?.length ?? 0) > 0) {
      return (
        <div className="prose max-w-none">
          Total Records: <b>{data.phosphorylation_data?.length ?? 0}</b>
          <PhosphorylationTable
            phosphorylation_data={data.phosphorylation_data ?? []}
          />
        </div>
      );
    } else {
      return (
        <div>
          <b>No Phosphorylation Information to Display</b>
        </div>
      );
    }
  })
  .build();

export const SNVViewResponseNode = MetaNode("SNVViewResponseNode")
  .meta({
    label: "SNV Information",
    description: "SNV Information",
    icon: [glygen_icon],
  })
  .codec(GlyGenProteinResponseArray)
  .view((data) => {
    const get_keywords = (keywords: Array<String>) => {
      const valid_keywords = keywords.filter(
        (kw) =>
          kw.toLowerCase() === "somatic" || kw.toLowerCase() === "germline",
      );
      if (valid_keywords.length === 0) return "NA";
      return valid_keywords.join("; ");
    };
    return (
      <div className="prose max-w-none">
        <table>
          <thead>
            <tr>
              <th>UniProt Ac</th>
              <th>Genomic Locus</th>
              <th>Ref NT</th>
              <th>Alt AT</th>
              <th>Start Pos</th>
              <th>End Pos</th>
              <th>Sequence</th>
              <th>Disease</th>
              <th>Keywords</th>
            </tr>
          </thead>
          <tbody>
            {data.map((protein, index) =>
              protein.snv.map((snv, snvIndex) => (
                <tr key={`${index}-${snvIndex}`}>
                  <td>
                    <a
                      href={`https://www.glygen.org/protein/${protein.uniprot.uniprot_canonical_ac}`}
                      target="_blank"
                      rel="noopener noreferrer"
                      style={{ color: "blue" }}
                    >
                      <u style={{ color: "blue" }}>
                        {protein.uniprot.uniprot_canonical_ac}
                      </u>
                    </a>
                  </td>
                  <td>{`Chr${snv.chr_id}:${snv.chr_pos}`}</td>
                  <td>{snv.ref_nt}</td>
                  <td>{snv.alt_nt}</td>
                  <td>{snv.start_pos}</td>
                  <td>{snv.end_pos}</td>
                  <td>{`${snv.sequence_org} -> ${snv.sequence_mut}`}</td>
                  <td>
                    {snv.disease && snv.disease.length > 0
                      ? snv.disease.map((disease, idx) => (
                          <span key={idx}>
                            {` ${disease.recommended_name.name}`} (
                            <a
                              href={`https://disease-ontology.org/term/${disease.disease_id}`}
                              target="_blank"
                              rel="noopener noreferrer"
                              style={{ color: "blue" }}
                            >
                              {disease.disease_id}
                            </a>
                            )
                          </span>
                        ))
                      : "NA"}
                  </td>
                  <td>{get_keywords(snv.keywords)}</td>
                </tr>
              )),
            )}
          </tbody>
        </table>
      </div>
    );
  })
  .build();

// --- PROCESS METANODES --- //

// Individual protein process metanode
export const GlyGenProtein = MetaNode("GGP")
  .meta({
    label: "Search GlyGen by Protein Name for Protein Products",
    description: "Find protein product records in GlyGen for the gene",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ protein_uniprot_canonical_ac: ProteinTerm })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    const protein_response = await get_single_protein_data(
      props.inputs.protein_uniprot_canonical_ac,
    );
    return protein_response;
  })
  .story((props) => ({
    abstract: `Next, the GlyGen database\\ref{doi:10.1093/glycob/cwz080} was searched to identify a relevant set of proteins that originate from.`,
  }))
  .build();

// Protein set protein process metanode
export const GlyGenProteinSet = MetaNode("GGPS")
  .meta({
    label: "Search GlyGen by Protein Name for Protein Products",
    description: "Find protein product records in GlyGen.",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ protein_uniprot_canonical_acs: ProteinSet })
  .output(GlyGenProteinSetResponseNode)
  .resolve(async (props) => {
    console.log(
      "===> Got protein(s) input: ",
      props.inputs.protein_uniprot_canonical_acs.set,
    );
    const protein_response = await get_protein_set_data(
      props.inputs.protein_uniprot_canonical_acs.set,
    );
    return protein_response;
  })
  .story((props) => ({
    abstract: `Next, the GlyGen database\\ref{doi:10.1093/glycob/cwz080} was searched to identify a relevant set of proteins that originate from.`,
  }))
  .build();

// Glycosylation table process metanode from the protein product
export const GlycosylationInformation = MetaNode("GlycosylationInformation")
  .meta({
    label: "Get Glycosylation Data from GlyGen Protein Products",
    description: "Glycosylation information for Glycoproteins",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glygenProteinResponse: GlyGenProteinResponseNode })
  .output(GlycosylationViewResponseNode)
  .resolve(async (props) => {
    const data = extract_specific_data(
      props.inputs.glygenProteinResponse,
      "glycosylation",
    );
    if (!glycosylation_check(data)) {
      throw new Error("Expected glycosylation data but got something else.");
    }
    return data;
  })
  .story((props) => ({
    abstract: "The glycosylation data was extracted from the GlyGen protein response and prepared for presentation in the view metanode.",
  }))
  .build();

// Phosphoprotein table process metanode from the protein product
export const PhosphorylationInformation = MetaNode("PhosphorylationInformation")
  .meta({
    label: "Get Phosphorylation Data from the GlyGen Protein Products",
    description: "Phosphorylation information for Phosphoproteins",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glyGenProteinResponse: GlyGenProteinResponseNode })
  .output(PhosphorylationViewResponseNode)
  .resolve(async (props) => {
    const data = extract_specific_data(
      props.inputs.glyGenProteinResponse,
      "phosphorylation",
    );
    if (!phosphorylation_check(data)) {
      throw new Error("Expected phosphorylation data but got something else.");
    }
    return data;
  })
  .story((props) => ({
    abstract: "The phosphorylation data was extracted from the GlyGen protein response and prepared for presentation in the view metanode.",
  }))
  .build();

// SNV table process metanode from the protein set
export const SNVInformation = MetaNode("SNVInformation")
  .meta({
    label: "Get SNV Data for the Protein Set",
    description: "SNV information",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glyGenSetProteinResponse: GlyGenProteinSetResponseNode })
  .output(SNVViewResponseNode)
  .resolve(async (props) => {
    const results = await Promise.all(
      props.inputs.glyGenSetProteinResponse.map((protein) => {
        const curr_accession = protein.uniprot.uniprot_canonical_ac;
        return get_single_protein_data(curr_accession);
      }),
    );
    console.log(results);
    return results;
  })
  .story((props) => "The SNV data is parsed from the GlyGen Protein data and prepared for presentation in the data view metanode.")
  .build();

// Links the Protein metanode chains to the Gene process metanodes
export const ProteinLink = MetaNode("ProteinLinkMetanode")
  .meta({
    label: "Get Additional Gene Data",
    description: "Protein link",
    icon: [glygen_icon],
  })
  .inputs({ glyGenProteinResponse: GlyGenProteinResponseNode })
  .output(GeneTerm)
  .resolve(async (props) => {
    const geneName = props.inputs.glyGenProteinResponse.gene.name;
    return geneName;
  })
  .story((props) => "The gene name was extracted from the protein response data in order to further explore the gene data.")
  .build();
