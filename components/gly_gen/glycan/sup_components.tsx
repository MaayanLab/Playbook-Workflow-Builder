import React from "react";
import { glygen_icon } from "@/icons";
import {
  GlycanResponse,
  GlyGenGlycanSetResponse,
  ClassificationArray,
  CrossRefArray,
  EnzymeArray,
} from "./data_models";
import { z } from "zod";

// Infer the Typescript type from the Zod schemas
type GlycanResponseType = z.infer<typeof GlycanResponse>;
type GlyGenGlycanSetResponse = z.infer<typeof GlyGenGlycanSetResponse>;
export type ClassificationArrayType = z.infer<typeof ClassificationArray>;
type CrossRefArrayType = z.infer<typeof CrossRefArray>;
type EnzymeArrayType = z.infer<typeof EnzymeArray>;

export function GlycanClassification({
  classification,
}: {
  classification: ClassificationArrayType;
}): String {
  /* Takes the glycan classification array and formats it into "Glycan Type / GlyCan Subtype"
   * delimited by "|".
   *
   * classification: the classification array in the format of ClassificationArray.
   *
   * Returns: a string formatted like "Glycan Type / Glycan Subtype" delmited by pipes ("|").
   */
  return classification
    .map((entry) => `${entry.type.name} / ${entry.subtype.name}`)
    .join(" | ");
}

export function GlycanCrossRef({ crossref }: { crossref: CrossRefArrayType }) {
  /* Takes the full cross reference array and pulls the links for just the PubChem Coumpound
   * and the PubChem Substance.
   *
   * crossref: the cross reference array in the format of CrossRefArray.
   *
   * Returns: formatted tsx tags in the format of "<database>: <id>" where the id is linked using
   * the corresponding url.
   */
  const pubMedCrossRefs = crossref.filter(
    (entry) =>
      ["PubChem Compound", "PubChem Substance"].includes(entry.database) &&
      entry.url,
  );

  if (pubMedCrossRefs.length === 0) {
    return null;
  }

  return (
    <div className="prose">
      {pubMedCrossRefs.map((entry) => (
        <div key={entry.id}>
          <b>{entry.database}: </b>
          <a
            href={entry.url}
            target="_blank"
            rel="noopener noreferrer"
            style={{ color: "blue" }}
          >
            <u style={{ color: "blue" }}>{entry.id}</u>
          </a>
        </div>
      ))}
    </div>
  );
}

export function GlycoenzymeTable({
  enzyme_data,
  is_preview = false,
}: {
  enzyme_data: EnzymeArrayType;
  is_preview?: boolean;
}) {
  /* Creates a table fo the glycoenzyme information.
   *
   * enzyme_data: An array of the glycoenzyme data in the format of EnzymeArray.
   * is_preview: Boolean that determines whether to create the full table or just the first 5 entries.
   *
   * Returns: A table of the extracted glycoenzyme data.
   */
  return (
    <div className="prose">
      <table>
        <thead>
          <tr>
            <th>UniProtKB Accession</th>
            <th>Gene Name</th>
            <th>Protein Name</th>
            <th>Organism</th>
          </tr>
        </thead>
        <tbody>
          {enzyme_data.map((entry, index) => (
            <tr key={index}>
              <td>
                <a
                  href={`http://www.glygen.org/protein/${entry.uniprot_canonical_ac}`}
                  target="_blank"
                  rel="noopener noreferrer"
                  style={{ color: "blue" }}
                >
                  <u style={{ color: "blue" }}>{entry.uniprot_canonical_ac}</u>
                </a>
              </td>
              <td>{entry.gene}</td>
              <td>{entry.protein_name}</td>
              <td>
                {entry.tax_name} ({entry.tax_common_name})
              </td>
            </tr>
          ))}
        </tbody>
      </table>
      {is_preview && (
        <div style={{ marginTop: "10px", fontSynthesis: "italic" }}>
          *This is a preview of the full glycoenzyme data table.
        </div>
      )}
    </div>
  );
}
