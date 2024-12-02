import React from "react";
import { z } from "zod";

// --- Supplementary Data Models for the Main Data Models --- //

// Glycosylation entry data models
const GlycosylationEntry = z.object({
  site_lbl: z.string(),
  site_category: z.string(),
  type: z.string(),
  glytoucan_ac: z.string(),
});

export const GlycosylationArray = z.array(GlycosylationEntry);
export const GlycosylationData = z.object({
  glycosylation: z.boolean().optional(),
  glycosylation_data: GlycosylationArray.optional(),
  protein_accession: z.string(),
});

// Phosphorylation entry data models
const PhosphorylationEntry = z.object({
  start_pos: z.number(),
  end_pos: z.number(),
  kinase_uniprot_canonical_ac: z.string().optional(),
  kinase_gene_name: z.string().optional(),
  residue: z.string(),
  comment: z.string().optional(),
});

export const PhosphorylationArray = z.array(PhosphorylationEntry);
export const PhosphorylationData = z.object({
  phosphorylation: z.boolean(),
  phosphorylation_data: PhosphorylationArray.optional(),
});

// Section stats entry data models
const SectionStatsEntry = z.object({
  table_id: z.string(),
  table_stats: z.array(
    z.object({
      field: z.string(),
      count: z.number(),
    }),
  ),
  sort_fields: z.array(z.string()).optional(),
});
export const SectionStatsArray = z.array(SectionStatsEntry);

// SNV entry data models
const SNVEntry = z.object({
  evidence: z.array(
    z.object({
      id: z.string(),
      database: z.string(),
      url: z.string(),
    }),
  ),
  glycoeffect: z.array(z.string()).optional(),
  start_pos: z.number(),
  end_pos: z.number(),
  sequence_org: z.string(),
  sequence_mut: z.string(),
  comment: z.string(),
  chr_id: z.string(),
  chr_pos: z.string(),
  ref_nt: z.string(),
  alt_nt: z.string(),
  site_lbl: z.string(),
  disease: z
    .array(
      z.object({
        disease_id: z.string(),
        recommended_name: z.object({
          id: z.string(),
          resource: z.string(),
          url: z.string(),
          name: z.string(),
          description: z.string().optional(),
        }),
      }),
    )
    .optional(),
  keywords: z.array(z.string()),
});
const SNVArray = z.array(SNVEntry);

// --- Main Protein Data Models --- //

// Formatted GlyGen protein response data model
export const GlyGenProteinResponse = z
  .object({
    snv: SNVArray.optional(),
    gene: z.object({
      name: z.string(),
      locus: z.object({
        chromosome: z.string(),
        start_pos: z.number(),
        end_pos: z.number(),
        strand: z.string(),
      }),
    }),
    uniprot: z.object({
      uniprot_id: z.string(),
      uniprot_canonical_ac: z.string(),
      length: z.number(),
    }),
    protein_names: z.object({
      name: z.string(),
    }),
    species: z.object({
      name: z.string(),
      common_name: z.string(),
      taxid: z.string(),
    }),
    glycoprotein: GlycosylationData,
    phosphorylation: PhosphorylationData,
    section_stats: SectionStatsArray,
  })
  .strip();
export const GlyGenProteinResponseArray = z.array(GlyGenProteinResponse)

// Formatted GlyGen protein set response data model
export const GlyGenProteinSetResponse = z.array(
  z.object({
    gene: z.object({
      name: z.string(),
    }),
    uniprot: z.object({
      uniprot_canonical_ac: z.string(),
    }),
    protein_names: z.object({
      name: z.string(),
    }),
    species: z.object({
      name: z.string(),
      taxid: z.string(),
    }),
    bools: z.object({
      total_n_glycosites: z.number(),
      total_o_glycosites: z.number(),
      reported_phosphosites: z.number(),
      reported_snv: z.number(),
    }),
  }),
);

// --- Supplementary Data Models for Handling API Responses in the Supplementary Functions -- //

// Name model for GlyGen Protein detail endpoint
export const ProteinNameAPIResponse = z.object({
  name: z.string(),
  resource: z.string().optional(),
  type: z.string(),
  id: z.string().optional(),
});

// Glycosylation model for GlyGen Protein detail endpoint
export const GlycosylationAPIResponse = z
  .object({
    site_lbl: z.string().optional(),
    site_category: z.string().optional(),
    type: z.string().optional(),
    glytoucan_ac: z.string().optional(),
  })
  .partial();

// Phosphorylation model for GlyGen Protein detail endpoint
export const PhosphorylationAPIResponse = z
  .object({
    start_pos: z.number().optional(),
    end_pos: z.number().optional(),
    kinase_uniprot_canonical_ac: z.string().optional(),
    kinase_gene_name: z.string().optional(),
    residue: z.string().optional(),
    comment: z.string().optional(),
  })
  .partial();
