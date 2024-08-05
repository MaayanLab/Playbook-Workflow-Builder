import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

//General RE info

export const RegulatoryElementInfoC = z.object({
    data: z.object({
      entId: z.string(),
      entType: z.string(),
      coordinates: z.object({
        "chromosome": z.string(),
        "start": z.number(),
        "end": z.number()
      }),
      ld: z.object({
        ENCODERegulatoryElementEvidence: z.array(z.object({
          ldhId: z.string()
        })),
        Variant: z.array(z.object({
          entId: z.string(),
          entIri: z.string(),
          ldhId: z.string(),
        }))
      }),
      ldFor: z.object({
        Gene: z.array(z.object({
          entId: z.string(),
          ldhId: z.string()
        }))
      })
    })
  })
  
  export type RegulatoryElementInfo = z.infer<typeof RegulatoryElementInfoC>

export async function myRegElemInfo_query(regElemId: string): Promise<RegulatoryElementInfo> {
    const res = await fetch(`https://ldh.genome.network/cfde/ldh/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
    return await res.json()
}


//Position data

export const RE_PositionalDataC = z.object({
    data: z.object({
      cCREQuery: z.array(z.object({
        coordinates: z.object({
          chromosome: z.string(),
          start: z.number(),
          end: z.number()
        })
      }))
    })
  });
type RE_PositionalData = z.infer<typeof RE_PositionalDataC>

export async function getRegElemPositionData(regElemId: string): Promise<RE_PositionalData> {
    let bodyString = '{\"query\":\"query CCRE{\\n  cCREQuery(assembly: \\"GRCh38\\", accession:\\"'+encodeURIComponent(regElemId)+'\\") {\\n  coordinates {\\n  chromosome\\n  start\\n  end\\n  }\\n  }\\n}"}';
    const res = await fetch(`https://ga.staging.wenglab.org/graphql`, {
        method: 'POST',
        headers: {
        'Accept-Encoding': 'gzip, deflate, br',
        'Content-Type': 'application/json',
        'Accept': 'application/json',
        'Connection': 'keep-alive',
        'DNT': '1',
        'Origin': 'https://ga.staging.wenglab.org'
        },
        body: bodyString
    })
    return await res.json()
}

//Unique Regions

export const RE_UniqueRegionC = z.array(
    z.object({
      '@id': z.string(),
      communityStandardTitle: z.string(),
      id: z.string(),
      type:  z.string()
    })
  );
export type RE_UniqueRegion = z.infer<typeof RE_UniqueRegionC>;

export const RE_UniqueRegionSetC = z.array(
    z.object({
      reId: z.string(),
      reference: RE_UniqueRegionC
    })
  );
export type RE_UniqueRegionSet = z.infer<typeof RE_UniqueRegionSetC>;
  
export async function getUniqueGenomicRegions(genomicRegion: string):  Promise<RE_UniqueRegion> {
    const res = await fetch(`https://reg.genome.network/reg/loc/desc/GRCh38%20(${encodeURIComponent(genomicRegion)})`)
    return await res.json()
}


//Unique Regions Cross Reference

export async function getUniqueGenomicRegionsCrossReference(genomicRegion: string): Promise<string>{
    const res = await fetch(`https://reg.clinicalgenome.org/coordinateTransform?glDesc=GRCh37%20(${encodeURIComponent(genomicRegion)})`)
    return await res.text();
}

export function reformatUniqueGenomicRegions(response: string): string[]{
    let genomes = ['GRCh38','GRCh37','NCBI36'];
    let referencesArray = [];
    for(let i in genomes){
      let g = genomes[i];
     let startIndx = response.indexOf(g);
     let endIndex = response.indexOf(")", startIndx);
     let reference = response.substring(startIndx,endIndex+1);
     if(reference != null && reference != ''){
      referencesArray.push(reference);
     }
    }
    return referencesArray;
  }
