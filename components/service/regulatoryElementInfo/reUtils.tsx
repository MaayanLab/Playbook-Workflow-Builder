import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'

//General RE info
export const RegulatoryElementInfoC = z.object({
    data: z.object({
      entId: z.string(),
      entType: z.string(),
      entContent: z.object({
        coordinates: z.object({
          chromosome: z.string(),
          start: z.number(),
          end: z.number()
        }),
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

export async function myRegElemInfoFromLDH(regElemId: string): Promise<RegulatoryElementInfo> {
    const res = await fetch(`https://ldh.genome.network/cfde/ldh/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
    return await res.json()
}

//positional data from Weng Lab
export const RE_PositionalDataWengLabC = z.object({
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
type RE_PositionalDataWengLab = z.infer<typeof RE_PositionalDataWengLabC>

export async function getRegElemPositionDataWengLab(regElemId: string): Promise<RE_PositionalDataWengLab> {
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

//positional data from LDH
export const RegulatoryElementPositionC = z.object({
    entId: z.string(),
    coordinates: z.object({
      "chromosome": z.string(),
      "start": z.number(),
      "end": z.number()
    }) 
})
export type RegulatoryElementPosition = z.infer<typeof RegulatoryElementPositionC>

export const RegulatoryElementSetPositionC = z.array( RegulatoryElementPositionC )
export type RegulatoryElementSetPosition = z.infer<typeof RegulatoryElementSetPositionC>

export async function getRegulatoryElementPosition(reIdentifier: string): Promise<RegulatoryElementPosition>{
  const ldhResponse = await myRegElemInfoFromLDH(reIdentifier);
  let coordinatesObject: any = null;
  if(ldhResponse != null && ldhResponse.data != null && ldhResponse.data.entContent != null){
    coordinatesObject = ldhResponse.data.entContent.coordinates
  }else{
    let wengLabReposne = await getRegElemPositionDataWengLab(reIdentifier);
    if(wengLabReposne != null && wengLabReposne.data.cCREQuery[0].coordinates != null){
      coordinatesObject = wengLabReposne.data.cCREQuery[0].coordinates;
    }
  }
  let reCoordinatesObj = {
    entId: ldhResponse.data.entId,
    coordinates: coordinatesObject
  };

  return reCoordinatesObj;
}

export async function getRegulatoryElementsSetPosition(regElemeIdsSet: string[]): Promise<RegulatoryElementSetPosition>{
  let reSetCoordinates: RegulatoryElementSetPosition = [];
  for(let i in regElemeIdsSet){
      let reId = regElemeIdsSet[i];
      let coordinatesObjFromLDH = await getRegulatoryElementPosition(reId);
      if(coordinatesObjFromLDH == null || coordinatesObjFromLDH.coordinates == null){
        continue;
      }
      reSetCoordinates.push(coordinatesObjFromLDH);
  }
  return reSetCoordinates;
}

const MyRegulatoryElementC = z.object({
  entId: z.string(),
  ldhId: z.string(),
  entContent: z.object({
    coordinates: z.object({    
      chromosome: z.string(),
      end: z.any(),
      start: z.any()
    })
  })
})
export type MyRegulatoryElement = z.infer<typeof MyRegulatoryElementC>

export const MyRegulatoryElementSetInfoC = z.array(
  MyRegulatoryElementC
)
export type MyRegulatoryElementSetInfo = z.infer<typeof MyRegulatoryElementSetInfoC>

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
