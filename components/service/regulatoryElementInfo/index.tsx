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
      throw new Error("Unable to get data from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
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
    label: 'Identify Genes In Vicinity',
    description: 'Identify genes in 10kbps distance of regulatory element.',
  })
  .inputs({  regulatoryElement: RegulatoryElementTerm   })
  .output(GeneSet)
  .resolve(async (props) => {
    const response = await myRegElemInfo_query(props.inputs.regulatoryElement);
    if(response == null || response.data == null){
      throw new Error("Unable to get data from Linked Data Hub API, please try again or wait a few minutes before the next atempt!");
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
    label: 'Identify Variants Within Regulatory Element',
    description: 'Retrieve registered variants in the region of given regulatory element from Allele Registry',
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
        <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
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

  const RE_UniqueRegionC = z.array(
      z.object({
        '@id': z.string(),
        communityStandardTitle: z.string(),
        id: z.string(),
        type:  z.string()
      })
    );
  type RE_UniqueRegion = z.infer<typeof RE_UniqueRegionC>;

  const RE_UniqueRegionSetC = z.array(
    z.object({
      reId: z.string(),
      reference: RE_UniqueRegionC
    })
  );
  type RE_UniqueRegionSet = z.infer<typeof RE_UniqueRegionSetC>;
  
  async function getUniqueGenomicRegions(genomicRegion: string):  Promise<RE_UniqueRegion> {
    const res = await fetch(`https://reg.test.genome.network/reg/loc/desc/GRCh38%20(${encodeURIComponent(genomicRegion)})`)
    return await res.json()
  }

  export const UniqueGenomicRegion = MetaNode('UniqueGenomicRegion')
  .meta({
    label: 'Unique Genomic Region',
    description: '',
    icon: [datafile_icon]
  })
  .codec(RE_UniqueRegionC)
  .view(uniqueRegions => {
    if(uniqueRegions == null || uniqueRegions.length == 0){
      return (
        <p>Nothing found on the queried position!</p>
      )
    }else{
      return( 
        <>
          <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
          <Table
            height={500}
            cellRendererDependencies={[uniqueRegions]}
            numRows={uniqueRegions.length}
            enableGhostCells
            enableFocusedCell
            downloads={{
              JSON: () => downloadBlob(new Blob([JSON.stringify(uniqueRegions)], { type: 'application/json;charset=utf-8' }), 'data.json')
            }}>
            <Column
              name="@id"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegions[row]['@id']}</Cell>}
            />
            <Column
              name="Title"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegions[row].communityStandardTitle}</Cell>}
            />
            <Column
              name="Id"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegions[row].id}</Cell>}
            />
            <Column
              name="Type"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegions[row].type}</Cell>}
            />
          </Table>
        </>
      )
    }
  })
  .build()

  export const GetUniqueNameForGenomicRegions = MetaNode('GetUniqueNameForGenomicRegions')
  .meta({
    label: 'Get unique name for genomic region(s)',
    description: 'Get unique name for genomic region(s).',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm })
  .output(UniqueGenomicRegion)
  .resolve(async (props) => {
    const rePositionData = await getRegElemPositionData(props.inputs.regulatoryElement);
    let positionString = rePositionData.data.cCREQuery[0].coordinates.chromosome+":"+rePositionData.data.cCREQuery[0].coordinates.start+"-"+rePositionData.data.cCREQuery[0].coordinates.end;
    let response = await getUniqueGenomicRegions(positionString);
    if(response == null){
      throw new Error("Unable to find any data for the queried position!");
    }
    return response;
  })
  .story(props => ({ abstract: `Get unique name for genomic region(s).` }))
  .build()

  export const UniqueGenomicRegionRESet = MetaNode('UniqueGenomicRegionRESet')
  .meta({
    label: 'Unique Genomic Region For Regulatory Elem. Set',
    description: '',
    icon: [datafile_icon]
  })
  .codec(RE_UniqueRegionSetC)
  .view(uniqueRegions => {
    if(uniqueRegions == null || uniqueRegions.length == 0){
      return (
        <p>Nothing found on the queried position!</p>
      )
    }else{
      return( 
        <>
          <p style={{fontSize: '14px'}}><b>Note:</b> In order to view all data, if avaliable, please expand the table rows!</p>
          <Table
            height={500}
            cellRendererDependencies={[uniqueRegions]}
            numRows={uniqueRegions.length}
            enableGhostCells
            enableFocusedCell
            downloads={{
              JSON: () => downloadBlob(new Blob([JSON.stringify(uniqueRegions)], { type: 'application/json;charset=utf-8' }), 'data.json')
            }}>
            <Column
              name="Regulatory Elem. id"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegions[row].reId}</Cell>}
            />
            <Column
              name="Title"
              cellRenderer={row => <Cell key={row+''}>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {uniqueRegions[row].reference.map(reference =>
                        <tr><td>{ reference.communityStandardTitle }</td></tr>
                    )}
                  </table>
              </Cell>}
            />
            <Column
              name="Id"
              cellRenderer={row => <Cell key={row+''}>                  
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {uniqueRegions[row].reference.map(reference =>
                      <tr><td>{ reference.id }</td></tr>
                    )}
                  </table>
                </Cell>}
            />
            <Column
              name="Type"
              cellRenderer={row => <Cell key={row+''}>
                  <table style={{borderCollapse: 'collapse', width:'100%'}}>
                    {uniqueRegions[row].reference.map(reference =>
                      <tr><td>{ reference.type }</td></tr>
                    )}
                  </table>
              </Cell>}
            />
          </Table>
        </>
      )
    }
  })
  .build()

  export const GetUniqueNameForGenomicRegionsRESet = MetaNode('GetUniqueNameForGenomicRegionsRESet')
  .meta({
    label: 'Get unique name for genomic region(s)',
    description: 'Get unique name for genomic region(s)',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElementSet: RegulatoryElementSet })
  .output(UniqueGenomicRegionRESet)
  .resolve(async (props) => {
    let regElemeIdsSet = props.inputs.regulatoryElementSet.set;
    let referenceArray: RE_UniqueRegionSet = [];
    for(let i in regElemeIdsSet){
      let rgId = regElemeIdsSet[i];
      const rePositionData = await getRegElemPositionData(rgId);

      let coordinates = null;
      if(rePositionData != null && rePositionData.data.cCREQuery[0] != null && rePositionData.data.cCREQuery[0].coordinates != null){
        coordinates = rePositionData.data.cCREQuery[0].coordinates;
      }else{
        continue;
      }

      if(coordinates != null){
        let positionString = coordinates.chromosome+":"+coordinates.start+"-"+coordinates.end;
        let response = await getUniqueGenomicRegions(positionString);
        if(response != null && response.length > 0){
          let tempObj = {
            reId: rgId,
            reference: response
          }    
          referenceArray.push(tempObj);  
        }
      }
    }
    return referenceArray;
  })
  .story(props => ({ abstract: `Get unique name for genomic region(s)` }))
  .build()


  async function getUniqueGenomicRegionsCrossReference(genomicRegion: string): Promise<string>{
    const res = await fetch(`https://reg.clinicalgenome.org/coordinateTransform?glDesc=GRCh37%20(${encodeURIComponent(genomicRegion)})`)
    return await res.text();
  }

  export const UniqueGenomicRegionCrossReference = MetaNode('UniqueGenomicRegionCrossReference')
  .meta({
    label: 'Unique Genomic Region Cross Reference',
    description: '',
    icon: [datafile_icon]
  })
  .codec(z.string())
  .view(uniqueRegionsCrossReference => { 
      if(uniqueRegionsCrossReference == null || uniqueRegionsCrossReference == ''){
        return(<p>Unable to find any unique region(s) for the queried position!</p>)
      }else{
        return(<p>{uniqueRegionsCrossReference}</p>)
      }
  })
  .build()

  export const GenomicRegionCoordinateTransformationAcrossReferences = MetaNode('GenomicRegionCoordinateTransformationAcrossReferences')
  .meta({
    label: 'Genomic region coordinate transformation across references',
    description: 'Get unique name for genomic region(s) cross refenrce: GRCh38, GRCh37 and NCBI36.',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElement: RegulatoryElementTerm })
  .output(UniqueGenomicRegionCrossReference)
  .resolve(async (props) => {
    const rePositionData = await getRegElemPositionData(props.inputs.regulatoryElement);
    let positionString = rePositionData.data.cCREQuery[0].coordinates.chromosome+":"+rePositionData.data.cCREQuery[0].coordinates.start+"-"+rePositionData.data.cCREQuery[0].coordinates.end;
    return await getUniqueGenomicRegionsCrossReference(positionString);
  })
  .story(props => ({ abstract: `Get unique name for genomic region(s) cross refenrce: GRCh38, GRCh37 and NCBI36.` }))
  .build()


  export const UniqueGenomicRegionCrossReferenceRESet = MetaNode('UniqueGenomicRegionCrossReferenceRESet')
  .meta({
    label: 'Unique Genomic Region Cross Reference',
    description: '',
    icon: [datafile_icon]
  })
  .codec(z.array(
    z.object({
      reId: z.string().optional(),
      reference: z.string().optional()
    }).optional()
  ).nullable().optional())
  .view(uniqueRegionsCrossReferenceList => { 
      if(uniqueRegionsCrossReferenceList == null || uniqueRegionsCrossReferenceList.length == 0){
        return(<p>Unable to find any unique region(s) for the queried positions!</p>)
      }else{
        return(
          <Table
            height={500}
            cellRendererDependencies={[uniqueRegionsCrossReferenceList]}
            numRows={uniqueRegionsCrossReferenceList.length}
            enableGhostCells
            enableFocusedCell
            downloads={{
              JSON: () => downloadBlob(new Blob([JSON.stringify(uniqueRegionsCrossReferenceList)], { type: 'application/json;charset=utf-8' }), 'data.json')
            }}>
            <Column
              name="Regulatory Element id"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegionsCrossReferenceList[row]?.reId}</Cell>}
            />
            <Column
              name="Cross Reference(s)"
              cellRenderer={row => <Cell key={row+''}>{uniqueRegionsCrossReferenceList[row]?.reference}</Cell>}
            />
          </Table>

        )
      }
  })
  .build()

  export const GenomicRegionCoordinateTransformationAcrossReferencesRE_Set = MetaNode('GenomicRegionCoordinateTransformationAcrossReferencesRE_Set')
  .meta({
    label: 'Genomic region coordinate transformation across references',
    description: 'Get unique name for genomic region(s) cross refenrce: GRCh38, GRCh37 and NCBI36.',
    icon: [linkeddatahub_icon],
  })
  .inputs({ regulatoryElementSet: RegulatoryElementSet })
  .output(UniqueGenomicRegionCrossReferenceRESet)
  .resolve(async (props) => {
    let regElemeIdsSet = props.inputs.regulatoryElementSet.set;
    let referenceArray = [];
    for(let i in regElemeIdsSet){
      let rgId = regElemeIdsSet[i];
      const rePositionData = await getRegElemPositionData(rgId);

      let coordinates = null;
      if(rePositionData != null && rePositionData.data.cCREQuery[0] != null && rePositionData.data.cCREQuery[0].coordinates != null){
        coordinates = rePositionData.data.cCREQuery[0].coordinates;
      }else{
        continue;
      }

      if(coordinates != null){
        let positionString = coordinates.chromosome+":"+coordinates.start+"-"+coordinates.end;
        let response = await getUniqueGenomicRegionsCrossReference(positionString); 
        let tempObj = {
          reId: rgId,
          reference: response
        }    
        referenceArray.push(tempObj);  
      }
    }
    return referenceArray;
  })
  .story(props => ({ abstract: `Get unique name for genomic region(s) cross refenrce: GRCh38, GRCh37 and NCBI36.` }))
  .build()
