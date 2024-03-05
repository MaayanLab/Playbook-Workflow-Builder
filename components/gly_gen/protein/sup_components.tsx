import React from "react";
import { glygen_icon } from "@/icons";
import { GlycosylationArray, PhosphorylationArray } from "./data_models";
import { z } from "zod";

// Infer the Typescript type from the Zod schemas
type GlycosylationArrayType = z.infer<typeof GlycosylationArray>;
type PhosphorylationArrayType = z.infer<typeof PhosphorylationArray>;

export function GlycosylationTable({
  glycosylation_data,
  is_preview = false,
}: {
  glycosylation_data: GlycosylationArrayType;
  is_preview?: boolean;
}) {
  /* Creates a table of the glycosylation information.
   *
   * glycosylation_data: an array of the glycosylation data in the format of GlycosylationArray.
   * is_preview: boolean that determines whether to create the full table or just the first 5 entries.
   *
   * Returns: A table of the extracted glycosylation data.
   */
  return (
    <div className="prose">
      <table>
        <thead>
          <tr>
            <th>Residue</th>
            <th>Glycosylation Site Category</th>
            <th>Glycosylation Type</th>
            <th>GlyTouCan Ac</th>
          </tr>
        </thead>
        <tbody>
          {glycosylation_data.map((entry, index) => (
            <tr key={index}>
              <td>{entry.site_lbl}</td>
              <td>{entry.site_category}</td>
              <td>{entry.type}</td>
              <td>{entry.glytoucan_ac}</td>
            </tr>
          ))}
        </tbody>
      </table>
      {is_preview && (
        <div style={{ marginTop: "10px", fontSynthesis: "italic" }}>
          *This is a preview of the full glycosylation data data.
        </div>
      )}
    </div>
  );
}

export function PhosphorylationTable({
  phosphorylation_data,
}: {
  phosphorylation_data: PhosphorylationArrayType;
}) {
  /* Creates a table of the phosphorylation information.
   *
   * phosphorylation_data: an array of the phosphorylation data in the format of PhosphorylationArray.
   *
   * Returns: A table of the extracted phosphorylation data.
   */
  return (
    <div>
      <table>
        <thead>
          <tr>
            <th>Residue</th>
            <th>Kinase Protein</th>
            <th>Kinase Gene</th>
            <th>Notes</th>
          </tr>
        </thead>
        <tbody>
          {phosphorylation_data.map((entry, index) => (
            <tr key={index}>
              <td>{entry.residue}</td>
              <td>{entry.kinase_uniprot_canonical_ac}</td>
              <td>{entry.kinase_gene_name}</td>
              <td>{entry.comment}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
