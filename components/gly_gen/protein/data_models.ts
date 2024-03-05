import React from 'react'
import { z } from 'zod';

// --- Supplementary Data Models for the Main Data Models --- //

// Glycosylation entry data models
const GlycosylationEntry = z.object({
  site_lbl: z.string(),
  site_category: z.string(),
  type: z.string(),
  glytoucan_ac: z.string()
});

export const GlycosylationArray = z.array(GlycosylationEntry)
export const GlycosylationData = z.object({
  glycosylation: z.boolean(),
  glycosylation_data: GlycosylationArray.optional(),
  protein_accession: z.string()
})

// Phosphorylation entry data models
const PhosphorylationEntry = z.object({
  start_pos: z.number(),
  end_pos: z.number(),
  kinase_uniprot_canonical_ac: z.string().optional(),
  kinase_gene_name: z.string().optional(),
  residue: z.string(),
  comment: z.string().optional()
})

export const PhosphorylationArray = z.array(PhosphorylationEntry)
export const PhosphorylationData = z.object({
  phosphorylation: z.boolean(),
  phosphorylation_data: PhosphorylationArray.optional()
})

// --- Main Protein Data Models --- //

// Formatted GlyGen protein response data model
export const GlyGenProteinResponse = z.object({
  gene: z.object({
    name: z.string(),
    locus: z.object({
      chromosome: z.string(),
      start_pos: z.number(),
      end_pos: z.number(),
      strand: z.string()
    })
  }),
  uniprot: z.object({
    uniprot_id: z.string(),
    uniprot_canonical_ac: z.string(),
    length: z.number()
  }),
  protein_names: z.object({
    name: z.string()
  }),
  species: z.object({
    name: z.string(),
    common_name: z.string(),
    taxid: z.string(),
  }),
  glycoprotein: GlycosylationData,
  phosphorylation: PhosphorylationData
}).strip();

// Formatted GlyGen protein set response data model 
export const GlyGenProteinSetResponse = z.array(
  z.object({
    gene: z.object({
      name: z.string()
    }),
    uniprot: z.object({
      uniprot_canonical_ac: z.string()
    }),
    protein_names: z.object({
      name: z.string()
    }),
    species: z.object({
      name: z.string(),
      taxid: z.string()
    }),
    bools: z.object({
      total_n_glycosites: z.number(),
      total_o_glycosites: z.number(),
      reported_phosphosites: z.number(),
      reported_snv: z.number()
    })
  })
)

// --- Supplementary Data Models for Handling API Responses -- //

// Name model for GlyGen Protein detail endpoint
export const ProteinNameAPIResponse = z.object({
    name: z.string(),
    resource: z.string().optional(),
    type: z.string(),
    id: z.string().optional()
})

// Glycosylation model for GlyGen Protein detail endpoint
export const GlycosylationAPIResponse = z.object({
  site_lbl: z.string().optional(),
  site_category: z.string().optional(),
  type: z.string().optional(),
  glytoucan_ac: z.string().optional()
}).partial();

// Phosphorylation model for GlyGen Protein detail endpoint
export const PhosphorylationAPIResponse = z.object({
  start_pos: z.number().optional(),
  end_pos: z.number().optional(),
  kinase_uniprot_canonical_ac: z.string().optional(),
  kinase_gene_name: z.string().optional(),
  residue: z.string().optional(),
  comment: z.string().optional()
}).partial();


