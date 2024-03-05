import React from "react";
import { glygen_icon } from "@/icons";
import {
  GlycanResponse,
  GlyGenGlycanSetResponse,
  ClassificationArray,
  CrossRefArray
} from "./data_models";
import { z } from "zod";

// Infer the Typescript type from the Zod schemas
type GlycanResponseType = z.infer<typeof GlycanResponse>;
type GlyGenGlycanSetResponse = z.infer<typeof GlyGenGlycanSetResponse>;
export type ClassificationArrayType = z.infer<typeof ClassificationArray>;
type CrossRefArrayType = z.infer<typeof CrossRefArray>;

export function GlycanClassification({classification}: {classification: ClassificationArrayType}): String {
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

export function GlycanCrossRef( {crossref}: {crossref: CrossRefArrayType} ) {
  /* Takes the full cross reference array and pulls the links for just the PubChem Coumpound
   * and the PubChem Substance. 
   *
   * crossref: the cross reference array in the format of CrossRefArray.
   *
   * Returns: formatted tsx tags in the format of "<database>: <id>" where the id is linked using
   * the corresponding url.
   */
  const pubMedCrossRefs = crossref.filter(entry => 
    ['PubChem Compound', 'PubChem Substance'].includes(entry.database) && entry.url
  )

  if (pubMedCrossRefs.length === 0) {
    return null
  }

  return (
    <div className="prose">
      {pubMedCrossRefs.map(entry => (
        <div key={entry.id}>
          <b>{entry.database}: </b>
          <a href={entry.url} target='_blank' rel='noopener noreferrer' style={{color: 'blue'}}>
            <u style={{color: 'blue'}}>{entry.id}</u>
          </a>
        </div>
      ))}
    </div>
  ) 

}
