import React from "react";
import { MetaNode } from "@/spec/metanode";
import { z } from "zod";
import { glygen_icon } from "@/icons";
import { GlycanTerm } from "@/components/core/input/term";
import { GlycanSet } from "@/components/core/input/set";
import {
  EnzymeData,
  GlycanResponse,
  GlyGenGlycanSetResponse,
} from "./data_models";
import {
  GlycanClassification,
  GlycanCrossRef,
  GlycoenzymeTable,
} from "./sup_components";
import {
  glycan_search,
  glycan_set_search_query,
  extract_enzyme_data,
  glycoenzyme_check,
} from "./sup_functions";

// --- Data Metanodes --- //

// Single glycan data metanode
export const GlycanViewResponseNode = MetaNode("GlycanViewResponse")
  .meta({
    label: "Glycan information",
    description: "Glycan information from GlyGen",
  })
  .codec(GlycanResponse)
  .view((data) => {
    const glyGenLink = `http://www.glygen.org/glycan/${data.glytoucan.glytoucan_ac}`;

    return (
      <div className="prose">
        <div>
          <span>GlyTouCan Accession: </span>
          <b>
            <a
              href={glyGenLink}
              target="_blank"
              rel="noopener nonreferrer"
              style={{ color: "blue" }}
            >
              {" "}
              <u style={{ color: "blue" }}>{data.glytoucan.glytoucan_ac}</u>
            </a>
          </b>
        </div>
        <div>
          Monoisotopic Mass: <b>{data.mass} Da</b>
        </div>
        <div>
          Monoisotopic Mass-pMe (Da): <b>{data.mass_pme} Da</b>
        </div>
        <div>
          Glycan Type / Glycan Subtype:{" "}
          <b>
            <GlycanClassification classification={data.classification} />
          </b>
        </div>
        <div>
          <GlycanCrossRef crossref={data.crossref} />
        </div>
        <div>Glycoenzyme: {data.enzyme.enzyme ? "True" : "False"}</div>
        <div>
          Glycan Image:
          <img
            src={`https://api.glygen.org/glycan/image/${data.glytoucan.glytoucan_ac}/`}
            alt="Glycan Image"
          />
        </div>
        <br />
        <div>
          {data.enzyme.enzyme &&
            data.enzyme.enzyme_data &&
            data.enzyme.enzyme_data.length > 0 && (
              <GlycoenzymeTable
                enzyme_data={
                  data.enzyme.enzyme_data.length > 5
                    ? data.enzyme.enzyme_data.slice(0, 5)
                    : data.enzyme.enzyme_data
                }
                is_preview={data.enzyme.enzyme_data.length > 5}
              />
            )}
        </div>
      </div>
    );
  })
  .build();

// Glycan set data metanode
export const GlyGenGlycanSetResponseNode = MetaNode("GlyGenGlycanSetResponse")
  .meta({
    label: "GlyGen Glycans",
    description: "Protein product records in GlyGen",
    icon: [glygen_icon],
  })
  .codec(GlyGenGlycanSetResponse)
  .view((data) => {
    return (
      <div>
        <table>
          <thead>
            <tr>
              <th>Glycan ID</th>
              <th>Glycan Image</th>
              <th>Hit Score</th>
              <th>Monoisotopic Mass</th>
              <th>Monoisotopic Mass-pMe (Da)</th>
              <th>No of Sugars</th>
              <th>No of Glycoproteins</th>
            </tr>
          </thead>
          <tbody>
            {data.map((entry, index) => (
              <tr key={index}>
                <td>{entry.glytoucan.glytoucan_ac}</td>
                <td>
                  <img
                    src={`https://api.glygen.org/glycan/image/${entry.glytoucan.glytoucan_ac}/`}
                    alt="Glycan Image"
                  />
                </td>
                <td>{entry.hit_score}</td>
                <td>{`${entry.mass}`}</td>
                <td>{entry.mass_pme}</td>
                <td>{entry.sugar_count}</td>
                <td>{entry.glycoprotein_count}</td>
                <td>{entry.associated_enzymes}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    );
  })
  .build();

// Glycoenzyme data metanode
export const GlycoenzymeViewResponseNode = MetaNode(
  "GlycoenzymeViewResponseNode",
)
  .meta({
    label: "GlyGen Glycoenzyme",
    description: "Glycoenzyme product records in GlyGen",
    icons: [glygen_icon],
  })
  .codec(EnzymeData)
  .view((data) => {
    if ((data.enzyme_data?.length ?? 0) > 0) {
      return (
        <div className="prose">
          Total Records: <b>{data.enzyme_data?.length ?? 0}</b>
          <GlycoenzymeTable enzyme_data={data.enzyme_data ?? []} />
        </div>
      );
    } else {
      return (
        <div>
          <b>No Enzyme Information to Display</b>
        </div>
      );
    }
  })
  .build();

// --- Process Metanodes --- //

// Single glycan process metanode
export const GlycanInformation = MetaNode("GlycanInformation")
  .meta({
    label: "Search GlyGen by GlyTouCan Accession",
    description: "Search for Glycan information",
    // icon: []
    pagerank: 2,
  })
  .inputs({ glycan: GlycanTerm })
  .output(GlycanViewResponseNode)
  .resolve(async (props) => {
    console.log(props.inputs.glycan);
    const glycan_data = await glycan_search(props.inputs.glycan);
    return glycan_data;
  })
  .story(
    (props) =>
      `The GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a information about ${props.inputs ? props.inputs.glycan : "the glycan"}.`,
  )
  .build();

// Glycan set process metanode
export const GlyGenGlycanSet = MetaNode("GGGS")
  .meta({
    label: "Search GlyGen by GlyTouCan Accession for Glycans",
    description: "Find glycan records in GlyGen.",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glycan_glytoucan_acc_set: GlycanSet })
  .output(GlyGenGlycanSetResponseNode)
  .resolve(async (props) => {
    console.log(
      "===> Got glycan set input: ",
      props.inputs.glycan_glytoucan_acc_set.set,
    );
    const glycan_response = await glycan_set_search_query(
      props.inputs.glycan_glytoucan_acc_set.set,
    );
    return glycan_response;
  })
  .story(
    (props) =>
      // TODO: re-write story sentence to make sense with protein term input (previous gene value removed to prevent `npm run build` error)
      `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of glycans that originate from.`,
  )
  .build();

// Glycoenzyme table process metanode from the glycan input
export const GlycoenzymeInformation = MetaNode("GlycoenzymeInformation")
  .meta({
    label: "Get Glycoenzyme Data from GlyGen Glycan Product",
    description: "Glycoenzyme information for Glycoenzymes",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ glyGenGlycanResponse: GlycanViewResponseNode })
  .output(GlycoenzymeViewResponseNode)
  .resolve(async (props) => {
    const data = extract_enzyme_data(props.inputs.glyGenGlycanResponse);
    if (!glycoenzyme_check(data)) {
      throw new Error("Expected glycoenzyme data but got something else.");
    }
    return data;
  })
  .story(
    (props) =>
      "The glycoenzyme data was extracted from the GlyGen glycan response and prepared for presentation in the view metanode.",
  )
  .build();
