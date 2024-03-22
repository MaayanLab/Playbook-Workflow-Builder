import React from "react";
import { MetaNode } from "@/spec/metanode";
import {
  GeneInfo,
  GeneInfoFromGeneTerm,
} from "@/components/service/mygeneinfo";
import { z } from "zod";
import { glygen_icon } from "@/icons";
import { ProteinTerm } from "@/components/core/input/term";
import { ProteinSet } from "@/components/core/input/set";
import {
  GlyGenProteinResponse,
  GlyGenProteinSetResponse,
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

// --- Data Metanodes --- //

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

// --- Process Metanodes --- //

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
  .story(
    (props) =>
      `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from.`,
  )
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
  .story(
    (props) =>
      `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from.`,
  )
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
  .story(
    (props) =>
      "The glycosylation data was extracted from the GlyGen protein response and prepared for presentation in the view metanode.",
  )
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
  .story(
    (props) =>
      "The phosphorylation data was extracted from the GlyGen protein response and prepared for presentation in the view metanode.",
  )
  .build();
