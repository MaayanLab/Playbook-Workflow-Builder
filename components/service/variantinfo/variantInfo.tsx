import { MetaNode } from '@/spec/metanode'
import { VariantTerm } from '@/components/core/input/term'
import { RegulatoryElementTerm } from '@/components/core/input/term'
import { VariantSet } from '@/components/core/input/set'
import { Table, Cell, Column} from '@/app/components/Table'
import { z } from 'zod'
import { linkeddatahub_icon } from '@/icons'
import { getRegElemPositionData } from '@/components/service/regulatoryElementInfo'
import { downloadBlob } from '@/utils/download'
import { resolveVarinatCaID, variantIdResolveErrorMessage, alleleRegRespErrorMessage, gitDataHubErroMessage } from './variantUtils'
import { xQTL_EvidenceEntContent } from './xQTLEvidences'
import { AlleleRegistryExternalRecordsTable, VarinatSetExternalRecordsInfo } from './externalRecords'
import { AlleleSpecificEvidenceInfoC } from './alleleSpecificEvidences'

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

const AlleleRegistryVariantInfoC = z.object({
    '@id': z.string(),
    entId: z.string(),
    externalRecords: z.any({}).optional()
  })
export type AlleleRegistryVariantInfo = z.infer<typeof AlleleRegistryVariantInfoC>

export async function getAlleleRegistryVariantInfo(variantId: string): Promise<AlleleRegistryVariantInfo> {
  const res = await fetch(`https://reg.genome.network/allele/${encodeURIComponent(variantId)}`);
  return await res.json()
}

const AlleleRegistryVariantSetInfoC = z.array(
  AlleleRegistryVariantInfoC
);


export const VariantInfo = MetaNode('VariantInfo')
  .meta({
    label: 'Variant Information',
    description: 'A Variant resolved with reg.clinicalgenome.org',
    icon: [linkeddatahub_icon],
  })
  .codec(AlleleRegistryVariantInfoC)
  .view(variantinfo => (
    <div className="prose max-w-none">
      <a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantinfo.entId}`}>{variantinfo.entId}</a> (variant)
    </div>
  ))
  .build()

export const VariantInfoFromVariantTermAlleleReg = MetaNode('VariantInfoFromVariantTermAlleleReg')
  .meta({
    label: 'Resolve Variant Info from Term (Allele Registry)',
    description: 'Resolve variant info from variant term using the Allele registry API.',
  })
  .inputs({ variant: VariantTerm })
  .output(VariantInfo)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const response = await getAlleleRegistryVariantInfo(varCaId);
    if(response == null){
      throw new Error(alleleRegRespErrorMessage);
    }

    response.entId = varCaId;
    return response;
  })
  .story(props => ``)
  .build()
  
  export const VariantSetInfo = MetaNode('VariantSetInfo')
  .meta({
    label: 'Variant Set Info',
    description: ''
  })
  .codec(AlleleRegistryVariantSetInfoC)
  .view(variantInfoSet => {
    return (
      <Table
        height={500}
        cellRendererDependencies={[variantInfoSet]}
        numRows={variantInfoSet.length}
        enableGhostCells
        enableFocusedCell
      >
        <Column
          name="Variant id (link)"
          cellRenderer={row => <Cell key={row+''}><a target="_blank" href={`https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_canonicalid?canonicalid=${variantInfoSet[row].entId}`}>{variantInfoSet[row].entId}</a></Cell>}
        />
      </Table>
    )
  }).build()

  export const VariantInfoFromVariantSet = MetaNode('VariantInfoFromVariantSet')
  .meta({
    label: `Variant Info From Variant Set (Allele reg.)`,
    description: "Get Variant Info from Allele Registry for a set of Variant CaId's."
  })
  .inputs({ variantset: VariantSet })
  .output(VariantSetInfo)
  .resolve(async (props) => {
    var varinatSetInpt = props.inputs.variantset.set;
    
    var varInfoSet = [];
    var varAlleleRegResponse = null;
    for(let indx in varinatSetInpt){
      varAlleleRegResponse = await getAlleleRegistryVariantInfo(varinatSetInpt[indx].trim());
      if(varAlleleRegResponse == null){
        continue;
      }
      varAlleleRegResponse.entId = varinatSetInpt[indx].trim();
      varInfoSet.push(varAlleleRegResponse);
    }

    if(varInfoSet.length == 0){
      throw new Error(alleleRegRespErrorMessage);
    }

    return varInfoSet;
  }).story(props =>
    `Get Allele Registry data for Variant Set.`
  ).build()

export const GetRegulatoryElementsForThisVariant = MetaNode('GetRegulatoryElementForThisVariant')
  .meta({
    label: 'Identify Variant And Regulatory Element Association',
    description: 'Get Regulatory Elements For This Variant ID',
  })
  .inputs({ variant: VariantTerm })
  .output(RegulatoryElementTerm)
  .resolve(async (props) => {
    var varCaId = await resolveVarinatCaID(props.inputs.variant);
    if(varCaId == null || varCaId == ''){
      throw new Error(variantIdResolveErrorMessage);
    }

    const response = await getGitDataHubVariantInfo(varCaId);
    if(response == null || response.data == null){
      throw new Error(gitDataHubErroMessage);
    }

    if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
      return response.data.ldFor.RegulatoryElement[0].entId;
    }
    return "N/A";
  })
  .story(props => ``)
  .build()

  export const REforVariantSetInfo = MetaNode('REforVariantSetInfo')
  .meta({
    label: 'Regulatory Elements For Variant Set Info',
    description: ''
  })
  .codec(
    z.array(z.object({
      caid: z.string(),
      reEntId: z.string(),
      rePosition: z.string()
    }))
  )
  .view( varAndReIdArray => {
    return ( 
      <Table
      height={500}
      cellRendererDependencies={[varAndReIdArray]}
      numRows={varAndReIdArray.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(varAndReIdArray)], { type: 'application/json;charset=utf-8' }), 'data.json')
      }}
      >
        <Column
          name="Variant CaId"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].caid }</Cell>}
        />
        <Column
          name="Regulatory Element EntId"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].reEntId }</Cell>}
        />
        <Column
          name="Regulatory Element Position"
          cellRenderer={row => <Cell key={row+''}>{ varAndReIdArray[row].rePosition }</Cell>}
        />
      </Table>  
    )
  })
  .build()

  export const GetRegulatoryElementsForVariantSet = MetaNode('GetRegulatoryElementsForVariantSet')
  .meta({
    label: `Regulatory Elements For Variant Set`,
    description: "Get Regulatory Elements for the Variant Set"
  })
  .inputs({ variantset: VariantSet })
  .output(REforVariantSetInfo)
  .resolve(async (props) => {
    var varinatSetInpt = props.inputs.variantset.set;
    
    let varAndRegulatoryElem = [];

    for(let indx in varinatSetInpt){
      let varCaID = varinatSetInpt[indx];
      const response = await getGitDataHubVariantInfo(varCaID);
      if(response == null || response.data == null){
        continue;
      }

      if(response.data.ldFor.RegulatoryElement != null && response.data.ldFor.RegulatoryElement.length == 1){
        let reEntId = response.data.ldFor.RegulatoryElement[0].entId;
        let tempObj = {
          'caid':varCaID,
          'reEntId': reEntId,
          'rePosition':""
        }

        const rePositionData = await getRegElemPositionData(reEntId);
        if(rePositionData != null && rePositionData.data.cCREQuery[0].coordinates != null){
          tempObj.rePosition = rePositionData.data.cCREQuery[0].coordinates.chromosome+": "+rePositionData.data.cCREQuery[0].coordinates.start+"-"+rePositionData.data.cCREQuery[0].coordinates.end;
        }
        varAndRegulatoryElem.push(tempObj);
      }
    }

    if(varAndRegulatoryElem.length == 0){
      throw new Error(gitDataHubErroMessage);
    }

    return varAndRegulatoryElem;
  }).story(props =>
    `Get Regulatory Elements id and posotion for Variant Set.`
  ).build()
