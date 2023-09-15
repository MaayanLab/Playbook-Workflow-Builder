import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { GeneSet, VariantSet } from '@/components/core/input/set'
import { z } from 'zod'
import { linkeddatahub_icon, datafile_icon } from '@/icons'
import { Table, Cell, Column} from '@/app/components/Table'

export const MyRegulatoryElementC = z.object({
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

export type MyRegulatoryElement = z.infer<typeof MyRegulatoryElementC>

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

export const RegulatoryElementInfo = MetaNode('RegulatoryElementInfo')
  .meta({
    label: 'Regulatory Element',
    description: 'Regulatory Element resolver',
    icon: [linkeddatahub_icon],
  })
  .codec(MyRegulatoryElementC)
  .view(regElem => (
    <div> {regElem.data.entId} {regElem.data.entType}<br></br>Position: {regElem.data.coordinates.chromosome}: {regElem.data.coordinates.start}-{regElem.data.coordinates.end} (GRCh38)</div>
  ))
  .build()

export async function myRegElemInfo_query(regElemId: string): Promise<MyRegulatoryElement> {
  const res = await fetch(`https://genboree.org/cfde-gene-dev/RegulatoryElement/id/${encodeURIComponent(regElemId)}`)
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

export const RegElementInfoFromRegElementTerm = MetaNode('RegElementInfoFromRegElementTerm')
  .meta({
    label: 'Resolve Regulatory Element Info from Term',
    description: 'Resolve Regulatory Element info from variant term',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm })
  .output(RegulatoryElementInfo)
  .resolve(async (props) => {
    const rePositionData = await getRegElemPositionData(props.inputs.regulatoryElement);
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(rePositionData.data.cCREQuery[0].coordinates != null){
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
  .story(props => `Additional information about the regulatory element${props.inputs ? ` ${props.inputs.regulatoryElement}` : ''} was resolved.`)
  .build()

export const GetGenesForRegulatoryElementInfo = MetaNode('GetGenesForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Genes for Regulatory Elements Info',
    description: 'Get Linked Genes For Regulatory Element Info',
  })
  .inputs({ regElemInfo: RegulatoryElementInfo  })
  .output(GeneSet)
  .resolve(async (props) => {
    let geneNames =  props.inputs.regElemInfo.data.ldFor.Gene.map(({ entId }) => entId);
    let geneSet = {
      description: 'Gene set for regulatory element '+props.inputs.regElemInfo.data.entId,
      set: geneNames
    };
    return geneSet;
  })
  .story(props => `Genes linked to the regulatory element${props.inputs ? ` ${props.inputs.regElemInfo.data.entId}` : ''} were resolved.`)
  .build()

export const GetVariantsForRegulatoryElementInfo = MetaNode('GetVariantListForRegulatoryElementInfo')
  .meta({
    label: 'Resolve Variants for Regulatory Elements Info',
    description: 'Get Linked Variants For Regulatory Element Info',
  })
  .inputs({ regElemInfo: RegulatoryElementInfo  })
  .output(VariantSet)
  .resolve(async (props) => {
    let variantNames =  props.inputs.regElemInfo.data.ld.Variant.map(({ entId }) => entId);
    let variantSet = {
      description: 'Variant set for regulatory element '+props.inputs.regElemInfo.data.entId,
      set: variantNames
    };
    return variantSet;
  })
  .story(props => `Variants linked to the regulatory element${props.inputs ? ` ${props.inputs.regElemInfo.data.entId}` : ''} were resolved.`)
  .build()

const MyRegulatoryElement = z.object({
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

export const MyRegulatoryElementSetInfoC = z.array(
  MyRegulatoryElement
)
export type MyRegulatoryElementSetInfo = z.infer<typeof MyRegulatoryElementSetInfoC>

export const RegulatoryElementSetInfo = MetaNode('RegulatoryElementSetInfo')
  .meta({
    label: 'RegulatoryElementSetInfo',
    description: '',
    icon: [datafile_icon]
  })
  .codec(MyRegulatoryElementSetInfoC)
  .view(regulatoryElementSet => {
    return( 
      <Table
        height={500}
        cellRendererDependencies={[regulatoryElementSet]}
        numRows={regulatoryElementSet.length}
        enableGhostCells
        enableFocusedCell>
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
    )
  })
  .build()
