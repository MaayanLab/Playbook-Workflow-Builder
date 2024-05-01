import { z } from 'zod'

const xQTL_EvidenceEntContent = z.object({
  GTExIri: z.string(),
  score: z.number(),
  sourceDescription: z.string(),
  eQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  sQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  esQTL: z.object({
    nes: z.string(),
    sig: z.string()
  }).optional(),
  type: z.string().optional()
});

export const xQTL_EvidenceArray = z.array(
  z.object({
    ldhId: z.string(),
    xQTLEntContent: xQTL_EvidenceEntContent
  })
);

export const GitHubVariantInfoC =  z.object({
    data: z.object({
      entId: z.string(),
      entType: z.string(),
      ld: z.object({
        AlleleSpecificEvidence: z.array(
          z.object({
            ldhId: z.string(),
            ldhIri: z.string(),
            entContent: z.object({
              sourceDescription: z.string(),
              AlleleSpecificity:z.any().optional()
            })
          })
        ),
        xqtlEvidence:  z.array(z.object({
          entId: z.string(),
          ldhId: z.string(),
          ldhIri: z.string(),
          entContent: xQTL_EvidenceEntContent
          })
        )
      }),
      ldFor: z.object({
        RegulatoryElement: z.array(
            z.object({ entId: z.string() })
          )
      })
    })
  })
 export type GitHubVariantInfo = z.infer<typeof GitHubVariantInfoC>

export async function getGitDataHubVariantInfo(variantId: string): Promise<GitHubVariantInfo> {
  const res = await fetch(`https://ldh.genome.network/cfde/ldh/Variant/id/${encodeURIComponent(variantId)}`)
  return await res.json()
}
