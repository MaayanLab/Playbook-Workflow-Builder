import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { RegulatoryElementTerm } from '@/components/core/term'
import { GeneSet, RegulatoryElementSet, VariantSet } from '@/components/core/set'
import { z } from 'zod'
import { linkeddatahub_icon, datafile_icon } from '@/icons'
import { Table, Cell, Column} from '@/app/components/Table'
import { downloadBlob } from '@/utils/download'

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

export const RegulatoryElementPosition = MetaNode('RegulatoryElementPosition')
  .meta({
    label: 'Regulatory Element position',
    description: 'Regulatory Element position',
    icon: [linkeddatahub_icon],
  })
  .codec(RegulatoryElementInfoC)
  .view(regElem => (
    <div className="prose max-w-none"> {regElem.data.entId} Regulatory Element<br></br>Position: {regElem.data.coordinates.chromosome}: {regElem.data.coordinates.start}-{regElem.data.coordinates.end} (GRCh38)</div>
  ))
  .build()

export async function myRegElemInfo_query(regElemId: string): Promise<RegulatoryElementInfo> {
  const res = await fetch(`https://ldh.genome.network/cfde/ldh/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
  return await res.json()
}

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

export const GetRegulatoryElementPosition = MetaNode('GetRegulatoryElementPosition')
  .meta({
    label: 'Regulatory Element Position Info',
    description: 'Regulatory Element Position Info',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm })
  .output(RegulatoryElementPosition)
  .resolve(async (props) => {
    const rePositionData = await getRegElemPositionData(props.inputs.regulatoryElement);
    
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(response == null || response.data == null){
      throw new Error("Unable to get data from Git Data Hub API, please try again or wait a few minutes before the next atempt!");
    }

    if(rePositionData != null && rePositionData.data.cCREQuery[0].coordinates != null){
      response.data.coordinates = rePositionData.data.cCREQuery[0].coordinates;
    }else{
      response.data.coordinates = {
          chromosome: "",
          start: 0,
          end: 0
      };
    }
    return response;
  })
  .story(props => ({ abstract: `Regulatory element genomic position.` }))
  .build()

export const GetGenesForRegulatoryElementInfo = MetaNode('GetGenesForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Genes for Regulatory Elements Info',
    description: 'Get Linked Genes For Regulatory Element Info',
  })
  .inputs({  regulatoryElement: RegulatoryElementTerm   })
  .output(GeneSet)
  .resolve(async (props) => {
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(response == null || response.data == null){
      throw new Error("Unable to get data from LInked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }

    let geneNames =  response.data.ldFor.Gene.map(({ entId }) => entId);
    let geneSet = {
      description: 'Gene set for regulatory element '+response.data.entId,
      set: geneNames
    };
    return geneSet;
  })
  .story(props => ({ abstract: `Genes linked to the regulatory element${props.inputs ? ` ${props.inputs.regulatoryElement}` : ''} were resolved.` }))
  .build()

export const GetVariantsForRegulatoryElementInfo = MetaNode('GetVariantListForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Variants for Regulatory Elements Info',
    description: 'Get Linked Variants For Regulatory Element Info',
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm  })
  .output(VariantSet)
  .resolve(async (props) => {
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(response == null || response.data == null){
      throw new Error("Unable to get data from LInked Data Hub API, please try again or wait a few minutes before the next atempt!");
    }

    let variantNames =  response.data.ld.Variant.map(({ entId }) => entId);
    let variantSet = {
      description: 'Variant set for regulatory element '+response.data.entId,
      set: variantNames
    };
    return variantSet;
  })
  .story(props => ({ abstract: `Variants linked to the regulatory element${props.inputs ? ` ${props.inputs.regulatoryElement}` : ''} were resolved.` }))
  .build()

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

export const RegulatoryElementSetInfo = MetaNode('RegulatoryElementSetInfo')
  .meta({
    label: 'Regulatory Element Set Info',
    description: '',
    icon: [datafile_icon]
  })
  .codec(MyRegulatoryElementSetInfoC)
  .view(regulatoryElementSet => {
    return( 
      <>
        <Table
          height={500}
          cellRendererDependencies={[regulatoryElementSet]}
          numRows={regulatoryElementSet.length}
          enableGhostCells
          enableFocusedCell
          downloads={{
            JSON: () => downloadBlob(new Blob([JSON.stringify(regulatoryElementSet)], { type: 'application/json;charset=utf-8' }), 'data.json')
          }}>
          <Column
            name="Entity id"
            cellRenderer={row => <Cell key={row+''}>{regulatoryElementSet[row].entId}</Cell>}
          />
          <Column
            name="Chromosome"
            cellRenderer={row => <Cell key={row+''}>{regulatoryElementSet[row].entContent.coordinates.chromosome}</Cell>}
          />
          <Column
            name="Start Pos."
            cellRenderer={row => <Cell key={row+''}>{regulatoryElementSet[row].entContent.coordinates.start}</Cell>}
          />
          <Column
            name="End Pos."
            cellRenderer={row => <Cell key={row+''}>{regulatoryElementSet[row].entContent.coordinates.end}</Cell>}
          />
        </Table>
      </>
    )
  })
  .build()

export const RegElementSetInfoFromRegElementTerm = MetaNode('RegElementSetInfoFromRegElementTerm')
  .meta({
    label: 'Resolve Regulatory Element Set Info from Term',
    description: 'Resolve Regulatory Element Set Info from term (id).',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElementSet: RegulatoryElementSet })
  .output(RegulatoryElementSetInfo)
  .resolve(async (props) => {
    let regElemeIdsSet = props.inputs.regulatoryElementSet.set;
    let response: MyRegulatoryElementSetInfo = [];
    for(let i in regElemeIdsSet){
      let rgId = regElemeIdsSet[i];
      const rePositionData = await getRegElemPositionData(rgId);

      let coordinates = null;
      if(rePositionData != null && rePositionData.data.cCREQuery[0] != null && rePositionData.data.cCREQuery[0].coordinates != null){
        coordinates = rePositionData.data.cCREQuery[0].coordinates;
      }else{
        coordinates = {
            chromosome: "",
            start: 0,
            end: 0
        };
      }
      let rgObj = {
        entId: rgId,
        ldhId: rgId,
        entContent: {
          coordinates: coordinates
        }
      }
      response.push(rgObj);
    }

    return response;
  })
  .story(props => ({ abstract: `Additional information about the regulatory elements was resolved.` }))
  .build()
