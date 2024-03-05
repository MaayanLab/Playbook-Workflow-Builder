import React from "react";
import { z } from "zod";

// --- Supplementary Data Models for the Main Data Models --- //

const ClassificationEntry = z.object({
  type: z.object({
    name: z.string(),
  }),
  subtype: z.object({
    name: z.string(),
  }),
});
export const ClassificationArray = z.array(ClassificationEntry);

const CrossRefEntry = z.object({
  id: z.string(),
  url: z.string().optional(),
  database: z.string(),
});
export const CrossRefArray = z.array(CrossRefEntry);

// --- Main Glycan Data Models --- //

export const GlycanResponse = z.object({
  glytoucan: z.object({
    glytoucan_ac: z.string(),
  }),
  mass: z.number(),
  mass_pme: z.number(),
  classification: ClassificationArray,
  crossref: CrossRefArray
});

export const GlyGenGlycanSetResponse = z.array(
  z.object({
    glytoucan: z.object({
      glytoucan_ac: z.string(),
    }),
    hit_score: z.number(),
    mass: z.number(),
    mass_pme: z.number(),
    sugar_count: z.number(),
    glycoprotein_count: z.number(),
    associated_enzymes: z.number(),
  }),
);
